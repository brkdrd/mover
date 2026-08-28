// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Collision/Shape/ConvexHullShape.h>
#include <Jolt/Physics/Collision/Shape/ScaleHelpers.h>
#include <Jolt/Physics/Collision/RayCast.h>
#include <Jolt/Physics/Collision/CastResult.h>
#include <Jolt/Physics/Collision/CollidePointResult.h>
#include <Jolt/Physics/Collision/TransformedShape.h>
#include <Jolt/Physics/Collision/CollideSoftBodyVertexIterator.h>
#include <Jolt/Geometry/ConvexHullBuilder.h>
#include <Jolt/Geometry/ClosestPoint.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>
#include <Jolt/Core/StringTools.h>
#include <Jolt/Core/StreamIn.h>
#include <Jolt/Core/StreamOut.h>
#include <Jolt/Core/UnorderedMap.h>

JPH_NAMESPACE_BEGIN

JPH_IMPLEMENT_SERIALIZABLE_VIRTUAL(ConvexHullShapeSettings)
{
	JPH_ADD_BASE_CLASS(ConvexHullShapeSettings, ConvexShapeSettings)

	JPH_ADD_ATTRIBUTE(ConvexHullShapeSettings, mPoints)
	JPH_ADD_ATTRIBUTE(ConvexHullShapeSettings, mMaxConvexRadius)
	JPH_ADD_ATTRIBUTE(ConvexHullShapeSettings, mMaxErrorConvexRadius)
	JPH_ADD_ATTRIBUTE(ConvexHullShapeSettings, mHullTolerance)
}

ShapeSettings::ShapeResult ConvexHullShapeSettings::Create() const
{
	if (mCachedResult.IsEmpty())
		Ref<Shape> shape = new ConvexHullShape(*this, mCachedResult);
	return mCachedResult;
}

ConvexHullShape::ConvexHullShape(const ConvexHullShapeSettings &inSettings, ShapeResult &outResult) :
	ConvexShape(EShapeSubType::ConvexHull, inSettings, outResult),
	mConvexRadius(inSettings.mMaxConvexRadius)
{
	using BuilderCell = ConvexHullBuilder::Cell;
	using Cells = ConvexHullBuilder::Cells;

	// Check convex radius
	if (mConvexRadius < 0.0f)
	{
		outResult.SetError("Invalid convex radius");
		return;
	}

	// Build convex hull (points are already 4D)
	const char *error = nullptr;
	ConvexHullBuilder builder(inSettings.mPoints);
	ConvexHullBuilder::EResult result = builder.Initialize(cMaxPointsInHull, inSettings.mHullTolerance, error);
	if (result != ConvexHullBuilder::EResult::Success && result != ConvexHullBuilder::EResult::MaxVerticesReached)
	{
		outResult.SetError(error);
		return;
	}
	const Cells &builder_cells = builder.GetCells();

	// Check the consistency of the resulting hull if we fully built it
	if (result == ConvexHullBuilder::EResult::Success)
	{
		ConvexHullBuilder::Cell *max_error_cell;
		float max_error_distance, coplanar_distance;
		int max_error_idx;
		builder.DetermineMaxError(max_error_cell, max_error_distance, max_error_idx, coplanar_distance);
		if (max_error_distance > 4.0f * max(coplanar_distance, inSettings.mHullTolerance)) // Coplanar distance could be bigger than the allowed tolerance if the points are far apart
		{
			outResult.SetError(StringFormat("Hull building failed, point %d had an error of %g (relative to tolerance: %g)", max_error_idx, (double)max_error_distance, double(max_error_distance / inSettings.mHullTolerance)));
			return;
		}
	}

	// Calculate center of mass and volume
	Vec4 com4;
	builder.GetCenterOfMassAndVolume(com4, mVolume);
	mCenterOfMass = com4;

	// Calculate the second-moment tensor M_ij = integral(x_i * x_j dV) of the solid 4-polytope.
	// Each boundary cell of a 4D hull is a tetrahedron; coning it to the center of mass (the origin
	// after centering) forms a pentachoron (4-simplex). We accumulate the covariance of each such
	// simplex using the canonical-simplex-transform trick: cov(S) = det(A) * A * C4 * A^T, where
	// A = [v1 v2 v3 v4] and C4 is the second-moment tensor of the standard 4-simplex.
	// For the standard n-simplex: integral(x_i^2) = V*2/((n+1)(n+2)), integral(x_i x_j) = V/((n+1)(n+2)),
	// with V = 1/n!. For n=4: V=1/24, diagonal = 1/360, off-diagonal = 1/720.
	Mat44 covariance_canonical(
		Vec4(1.0f / 360.0f, 1.0f / 720.0f, 1.0f / 720.0f, 1.0f / 720.0f),
		Vec4(1.0f / 720.0f, 1.0f / 360.0f, 1.0f / 720.0f, 1.0f / 720.0f),
		Vec4(1.0f / 720.0f, 1.0f / 720.0f, 1.0f / 360.0f, 1.0f / 720.0f),
		Vec4(1.0f / 720.0f, 1.0f / 720.0f, 1.0f / 720.0f, 1.0f / 360.0f));
	Mat44 covariance_matrix = Mat44::sZero();
	for (BuilderCell *cell : builder_cells)
	{
		// Each cell is a tetrahedron (4 vertices); the 5th pentachoron vertex is the center of mass (origin).
		Vec4 v1 = inSettings.mPoints[cell->mIdx[0]] - mCenterOfMass;
		Vec4 v2 = inSettings.mPoints[cell->mIdx[1]] - mCenterOfMass;
		Vec4 v3 = inSettings.mPoints[cell->mIdx[2]] - mCenterOfMass;
		Vec4 v4 = inSettings.mPoints[cell->mIdx[3]] - mCenterOfMass;

		Mat44 a(v1, v2, v3, v4);
		float det_a = v1.Dot(Vec4::sCross(v2, v3, v4)); // 4x4 determinant = generalized scalar quad product
		Mat44 c = det_a * (a * covariance_canonical * a.Transposed());
		covariance_matrix += c;
	}

	// mInertia stores the raw second-moment tensor (density 1), matching BoxShape/CapsuleShape.
	// TODO(4D): validate this pentachoron decomposition against a brute-force numerical integral.
	mInertia = covariance_matrix;

	// Convert cells from the builder to our internal representation
	using VtxMap = UnorderedMap<int, uint8>;
	VtxMap vertex_map;
	vertex_map.reserve(VtxMap::size_type(inSettings.mPoints.size()));
	for (BuilderCell *builder_cell : builder_cells)
	{
		JPH_ASSERT(mVertexIdx.size() <= 0xFFFF);
		uint16 first_vertex = (uint16)mVertexIdx.size();
		uint16 num_vertices = 0;

		// Loop over the 4 vertices in the cell
		for (int vi = 0; vi < 4; ++vi)
		{
			uint8 new_idx;
			int original_idx = builder_cell->mIdx[vi];
			VtxMap::iterator m = vertex_map.find(original_idx);
			if (m != vertex_map.end())
			{
				new_idx = m->second;
			}
			else
			{
				Vec4 p = inSettings.mPoints[original_idx] - mCenterOfMass;
				mLocalBounds.Encapsulate(p);
				JPH_ASSERT(mPoints.size() <= 0xff);
				new_idx = (uint8)mPoints.size();
				mPoints.push_back({ p });
				vertex_map[original_idx] = new_idx;
			}

			JPH_ASSERT(mVertexIdx.size() < 0xffff);
			mVertexIdx.push_back(new_idx);
			num_vertices++;
		}

		// Add face (4 vertices per cell)
		mFaces.push_back({ first_vertex, num_vertices });

		// Add plane
		Vec4 cell_centroid = builder_cell->mCentroid - com4;
		Vec4 cell_normal_normalized = builder_cell->mNormal;
		float nlen = cell_normal_normalized.Length();
		if (nlen > 0.0f) cell_normal_normalized /= nlen;
		Plane plane = Plane::sFromPointAndNormal(cell_centroid, cell_normal_normalized);
		mPlanes.push_back(plane);
	}

	// Test if GetSupportFunction can support this many points
	if (mPoints.size() > cMaxPointsInHull)
	{
		outResult.SetError(StringFormat("Internal error: Too many points in hull (%u), max allowed %d", (uint)mPoints.size(), cMaxPointsInHull));
		return;
	}

	for (int p = 0; p < (int)mPoints.size(); ++p)
	{
		// For each point, find faces that use the point
		Array<int> faces;
		for (int f = 0; f < (int)mFaces.size(); ++f)
		{
			const Face &face = mFaces[f];
			for (int v = 0; v < face.mNumVertices; ++v)
				if (mVertexIdx[face.mFirstVertex + v] == p)
				{
					faces.push_back(f);
					break;
				}
		}

		if (faces.size() < 2)
		{
			outResult.SetError("A point must be connected to 2 or more faces!");
			return;
		}

		// Record up to 3 adjacent faces for this point.
		// TODO(4D): in 3D these were chosen to form the largest tetrahedron of normals (used by the
		// convex-radius vertex shift). In 4D a vertex meets >=4 hyperplanes and the shift needs a 4x4
		// solve; since the convex radius is currently deferred (forced to 0 below), we just record the
		// first few adjacent faces so the data structure is populated.
		Point &point = mPoints[p];
		point.mNumFaces = (int)min<size_t>(faces.size(), 3);
		for (int i = 0; i < point.mNumFaces; ++i)
			point.mFaces[i] = faces[i];
	}

	// TODO(4D): the convex radius (a GJK margin optimization) is deferred. Shrinking the hull inward by
	// the convex radius requires intersecting the >=4 hyperplanes that meet at each 4D vertex (a 4x4
	// solve), plus a 4D thickness/sharp-edge analysis. For v1 we ship a zero-margin ("sharp") hull,
	// which is fully correct, just numerically less smooth in GJK. See GetSupportFunction.
	mConvexRadius = 0.0f;

	// Calculate the inner radius by getting the minimum distance from the origin to the planes of the hull
	mInnerRadius = FLT_MAX;
	for (const Plane &p : mPlanes)
		mInnerRadius = min(mInnerRadius, -p.GetConstant());
	mInnerRadius = max(0.0f, mInnerRadius); // Clamp against zero, this should do nothing as the shape is centered around the center of mass but for flat convex hulls there may be numerical round off issues

	outResult.Set(this);
}

MassProperties ConvexHullShape::GetMassProperties() const
{
	MassProperties p;

	float density = GetDensity();

	// Calculate mass
	p.mMass = density * mVolume;

	// Calculate the second-moment tensor (full 4x4, no homogeneous padding element)
	p.mInertia = density * mInertia;

	return p;
}

Vec4 ConvexHullShape::GetSurfaceNormal(const SubShapeID &inSubShapeID, Vec4Arg inLocalSurfacePosition) const
{
	JPH_ASSERT(inSubShapeID.IsEmpty(), "Invalid subshape ID");

	const Plane &first_plane = mPlanes[0];
	Vec4 best_normal = first_plane.GetNormal();
	float best_dist = abs(first_plane.SignedDistance(inLocalSurfacePosition));

	// Find the face that has the shortest distance to the surface point
	for (Array<Face>::size_type i = 1; i < mFaces.size(); ++i)
	{
		const Plane &plane = mPlanes[i];
		Vec4 plane_normal = plane.GetNormal();
		float dist = abs(plane.SignedDistance(inLocalSurfacePosition));
		if (dist < best_dist)
		{
			best_dist = dist;
			best_normal = plane_normal;
		}
	}

	return best_normal;
}

class ConvexHullShape::HullNoConvex final : public Support
{
public:
	explicit				HullNoConvex(float inConvexRadius) :
		mConvexRadius(inConvexRadius)
	{
		static_assert(sizeof(HullNoConvex) <= sizeof(SupportBuffer), "Buffer size too small");
		JPH_ASSERT(IsAligned(this, alignof(HullNoConvex)));
	}

	virtual Vec4			GetSupport(Vec4Arg inDirection) const override
	{
		// Find the point with the highest projection on inDirection
		float best_dot = -FLT_MAX;
		Vec4 best_point = Vec4::sZero();

		for (Vec4 point : mPoints)
		{
			// Check if its support is bigger than the current max
			float dot = point.Dot(inDirection);
			if (dot > best_dot)
			{
				best_dot = dot;
				best_point = point;
			}
		}

		return best_point;
	}

	virtual float			GetConvexRadius() const override
	{
		return mConvexRadius;
	}

	using PointsArray = StaticArray<Vec4, cMaxPointsInHull>;

	inline PointsArray &	GetPoints()
	{
		return mPoints;
	}

	const PointsArray &		GetPoints() const
	{
		return mPoints;
	}

private:
	float					mConvexRadius;
	PointsArray				mPoints;
};

class ConvexHullShape::HullWithConvex final : public Support
{
public:
	explicit				HullWithConvex(const ConvexHullShape *inShape) :
		mShape(inShape)
	{
		static_assert(sizeof(HullWithConvex) <= sizeof(SupportBuffer), "Buffer size too small");
		JPH_ASSERT(IsAligned(this, alignof(HullWithConvex)));
	}

	virtual Vec4			GetSupport(Vec4Arg inDirection) const override
	{
		// Find the point with the highest projection on inDirection
		float best_dot = -FLT_MAX;
		Vec4 best_point = Vec4::sZero();

		for (const Point &point : mShape->mPoints)
		{
			// Check if its support is bigger than the current max
			float dot = point.mPosition.Dot(inDirection);
			if (dot > best_dot)
			{
				best_dot = dot;
				best_point = point.mPosition;
			}
		}

		return best_point;
	}

	virtual float			GetConvexRadius() const override
	{
		return 0.0f;
	}

private:
	const ConvexHullShape *	mShape;
};

class ConvexHullShape::HullWithConvexScaled final : public Support
{
public:
							HullWithConvexScaled(const ConvexHullShape *inShape, Vec4Arg inScale) :
		mShape(inShape),
		mScale(inScale)
	{
		static_assert(sizeof(HullWithConvexScaled) <= sizeof(SupportBuffer), "Buffer size too small");
		JPH_ASSERT(IsAligned(this, alignof(HullWithConvexScaled)));
	}

	virtual Vec4			GetSupport(Vec4Arg inDirection) const override
	{
		// Find the point with the highest projection on inDirection
		float best_dot = -FLT_MAX;
		Vec4 best_point = Vec4::sZero();

		for (const Point &point : mShape->mPoints)
		{
			// Calculate scaled position
			Vec4 pos = mScale * point.mPosition;

			// Check if its support is bigger than the current max
			float dot = pos.Dot(inDirection);
			if (dot > best_dot)
			{
				best_dot = dot;
				best_point = pos;
			}
		}

		return best_point;
	}

	virtual float			GetConvexRadius() const override
	{
		return 0.0f;
	}

private:
	const ConvexHullShape *	mShape;
	Vec4					mScale;
};

const ConvexShape::Support *ConvexHullShape::GetSupportFunction([[maybe_unused]] ESupportMode inMode, SupportBuffer &inBuffer, Vec4Arg inScale) const
{
	// TODO(4D): the convex radius is currently forced to 0 at construction (shrinking the hull inward
	// needs a 4-hyperplane vertex solve; deferred). So the support is the raw (zero-margin) hull for
	// every support mode, and HullNoConvex (the vertex-shrink support) is unused.
	JPH_ASSERT(mConvexRadius == 0.0f);
	if (ScaleHelpers::IsNotScaled(inScale))
		return new (&inBuffer) HullWithConvex(this);
	else
		return new (&inBuffer) HullWithConvexScaled(this, inScale);
}

void ConvexHullShape::GetSupportingFace(const SubShapeID &inSubShapeID, Vec4Arg inDirection, Vec4Arg inScale, RMat44Arg inCenterOfMassTransform, SupportingFace &outVertices) const
{
	JPH_ASSERT(inSubShapeID.IsEmpty(), "Invalid subshape ID");

	Vec4 inv_scale = inScale.Reciprocal();

	// Need to transform the plane normals using inScale
	// Transforming a direction with matrix M is done through multiplying by (M^-1)^T
	// In this case M is a diagonal matrix with the scale vector, so we need to multiply our normal by 1 / scale and renormalize afterwards
	Vec4 plane0_normal = inv_scale * mPlanes[0].GetNormal();
	float best_dot = plane0_normal.Dot(inDirection) / plane0_normal.Length();
	int best_face_idx = 0;

	for (Array<Plane>::size_type i = 1; i < mPlanes.size(); ++i)
	{
		Vec4 plane_normal = inv_scale * mPlanes[i].GetNormal();
		float dot = plane_normal.Dot(inDirection) / plane_normal.Length();
		if (dot < best_dot)
		{
			best_dot = dot;
			best_face_idx = (int)i;
		}
	}

	// Get vertices
	const Face &best_face = mFaces[best_face_idx];
	const uint8 *first_vtx = mVertexIdx.data() + best_face.mFirstVertex;
	const uint8 *end_vtx = first_vtx + best_face.mNumVertices;

	// If we have more than 1/2 the capacity of outVertices worth of vertices, we start skipping vertices (note we can't fill the buffer completely since extra edges will be generated by clipping).
	// TODO: This really needs a better algorithm to determine which vertices are important!
	int max_vertices_to_return = outVertices.capacity() / 2;
	int delta_vtx = (int(best_face.mNumVertices) + max_vertices_to_return) / max_vertices_to_return;

	// Calculate transform with scale
	RMat44 transform = inCenterOfMassTransform.PreScaled(inScale);

	// TODO(4D): a 4D hull facet is a 3D polytope; the winding-flip below is the 3D heuristic and the
	// contact-manifold face pipeline (conceptual gap 1) still needs a real 4D treatment.
	if (ScaleHelpers::IsInsideOut(inScale))
	{
		// Flip winding of supporting face
		for (const uint8 *v = end_vtx - 1; v >= first_vtx; v -= delta_vtx)
			outVertices.push_back(Vec4(transform * mPoints[*v].mPosition));
	}
	else
	{
		// Normal winding of supporting face
		for (const uint8 *v = first_vtx; v < end_vtx; v += delta_vtx)
			outVertices.push_back(Vec4(transform * mPoints[*v].mPosition));
	}
}

void ConvexHullShape::GetSubmergedVolume(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const Plane &inSurface, float &outTotalVolume, float &outSubmergedVolume, Vec4 &outCenterOfBuoyancy, [[maybe_unused]] RVec4Arg inBaseOffset) const
{
	// Trivially calculate total 4-volume (scaled by all 4 axes)
	Vec4 abs_scale = inScale.Abs();
	outTotalVolume = mVolume * abs_scale.GetX() * abs_scale.GetY() * abs_scale.GetZ() * abs_scale.GetW();

	// TODO(4D): PolyhedronSubmergedVolumeCalculator is a 3D polyhedron clip. Computing the submerged
	// 4-volume of a 4-polytope against a hyperplane needs a 4D clip (part of the deferred buoyancy /
	// contact-manifold work). For now report nothing submerged.
	(void)inCenterOfMassTransform; (void)inSurface;
	outSubmergedVolume = 0.0f;
	outCenterOfBuoyancy = Vec4::sZero();
}

#ifdef JPH_DEBUG_RENDERER
void ConvexHullShape::Draw(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inUseMaterialColors, bool inDrawWireframe) const
{
	if (mGeometry == nullptr)
	{
		Array<DebugRenderer::Triangle> triangles;
		for (const Face &f : mFaces)
		{
			const uint8 *first_vtx = mVertexIdx.data() + f.mFirstVertex;
			const uint8 *end_vtx = first_vtx + f.mNumVertices;

			// Draw first triangle of polygon
			Vec3 v0 = mPoints[first_vtx[0]].mPosition;
			Vec3 v1 = mPoints[first_vtx[1]].mPosition;
			Vec3 v2 = mPoints[first_vtx[2]].mPosition;
			Vec3 uv_direction = (v1 - v0).Normalized();
			triangles.push_back({ v0, v1, v2, Color::sWhite, v0, uv_direction });

			// Draw any other triangles in this polygon
			for (const uint8 *v = first_vtx + 3; v < end_vtx; ++v)
				triangles.push_back({ v0, mPoints[*(v - 1)].mPosition, mPoints[*v].mPosition, Color::sWhite, v0, uv_direction });
		}
		mGeometry = new DebugRenderer::Geometry(inRenderer->CreateTriangleBatch(triangles), GetLocalBounds());
	}

	// Test if the shape is scaled inside out
	DebugRenderer::ECullMode cull_mode = ScaleHelpers::IsInsideOut(inScale)? DebugRenderer::ECullMode::CullFrontFace : DebugRenderer::ECullMode::CullBackFace;

	// Determine the draw mode
	DebugRenderer::EDrawMode draw_mode = inDrawWireframe? DebugRenderer::EDrawMode::Wireframe : DebugRenderer::EDrawMode::Solid;

	// Draw the geometry
	Color color = inUseMaterialColors? GetMaterial()->GetDebugColor() : inColor;
	RMat44 transform = inCenterOfMassTransform.PreScaled(inScale);
	inRenderer->DrawGeometry(transform, color, mGeometry, cull_mode, DebugRenderer::ECastShadow::On, draw_mode);

	// Draw the outline if requested
	if (sDrawFaceOutlines)
		for (const Face &f : mFaces)
		{
			const uint8 *first_vtx = mVertexIdx.data() + f.mFirstVertex;
			const uint8 *end_vtx = first_vtx + f.mNumVertices;

			// Draw edges of face
			inRenderer->DrawLine(transform * mPoints[*(end_vtx - 1)].mPosition, transform * mPoints[*first_vtx].mPosition, Color::sGrey);
			for (const uint8 *v = first_vtx + 1; v < end_vtx; ++v)
				inRenderer->DrawLine(transform * mPoints[*(v - 1)].mPosition, transform * mPoints[*v].mPosition, Color::sGrey);
		}
}

void ConvexHullShape::DrawShrunkShape(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale) const
{
	// Get the shrunk points
	SupportBuffer buffer;
	const HullNoConvex *support = mConvexRadius > 0.0f? static_cast<const HullNoConvex *>(GetSupportFunction(ESupportMode::ExcludeConvexRadius, buffer, inScale)) : nullptr;

	RMat44 transform = inCenterOfMassTransform * Mat44::sScale(inScale);

	for (int p = 0; p < (int)mPoints.size(); ++p)
	{
		const Point &point = mPoints[p];
		RVec3 position = transform * point.mPosition;
		RVec3 shrunk_point = support != nullptr? transform * support->GetPoints()[p] : position;

		// Draw difference between shrunk position and position
		inRenderer->DrawLine(position, shrunk_point, Color::sGreen);

		// Draw face normals that are contributing
		for (int i = 0; i < point.mNumFaces; ++i)
			inRenderer->DrawLine(position, position + 0.1f * mPlanes[point.mFaces[i]].GetNormal(), Color::sYellow);

		// Draw point index
		inRenderer->DrawText3D(position, ConvertToString(p), Color::sWhite, 0.1f);
	}
}
#endif // JPH_DEBUG_RENDERER

bool ConvexHullShape::CastRayHelper(const RayCast &inRay, float &outMinFraction, float &outMaxFraction) const
{
	// Clip the ray against all bounding hyperplanes (dimension-agnostic).
	// TODO(4D): 3D Jolt had a special case for flat hulls (mFaces.size() == 2) that tested polygon
	// edges. In 4D a "flat" (degenerate) hull lives in a 3-flat and would need an in-cell point test;
	// that edge case is deferred - the halfspace clip below handles all non-degenerate hulls.
	int fractions_set = 0;
	bool all_inside = true;
	float min_fraction = 0.0f, max_fraction = 1.0f + FLT_EPSILON;
	for (const Plane &p : mPlanes)
	{
		// Check if the ray origin is behind this plane
		Vec4 plane_normal = p.GetNormal();
		float distance_to_plane = inRay.mOrigin.Dot(plane_normal) + p.GetConstant();
		bool is_outside = distance_to_plane > 0.0f;
		all_inside &= !is_outside;

		// Check if plane is parallel to ray
		float direction_projection = inRay.mDirection.Dot(plane_normal);
		if (abs(direction_projection) >= 1.0e-12f)
		{
			// Get intersection fraction between ray and plane
			float fraction = -distance_to_plane / direction_projection;

			// Update interval of ray that is inside the hull
			if (direction_projection < 0.0f)
			{
				min_fraction = max(fraction, min_fraction);
				fractions_set |= 1;
			}
			else
			{
				max_fraction = min(fraction, max_fraction);
				fractions_set |= 2;
			}
		}
		else if (is_outside)
			return false; // Outside the plane and parallel, no hit!
	}

	// Test if both min and max have been set
	if (fractions_set == 3)
	{
		// Output fractions
		outMinFraction = min_fraction;
		outMaxFraction = max_fraction;

		// Test if the infinite ray intersects with the hull (the length will be checked later)
		return min_fraction <= max_fraction && max_fraction >= 0.0f;
	}
	else
	{
		// Degenerate case, either the ray is parallel to all planes or the ray has zero length
		outMinFraction = 0.0f;
		outMaxFraction = 1.0f + FLT_EPSILON;

		// Return if the origin is inside the hull
		return all_inside;
	}
}

bool ConvexHullShape::CastRay(const RayCast &inRay, const SubShapeIDCreator &inSubShapeIDCreator, RayCastResult &ioHit) const
{
	// Determine if ray hits the shape
	float min_fraction, max_fraction;
	if (CastRayHelper(inRay, min_fraction, max_fraction)
		&& min_fraction < ioHit.mFraction) // Check if this is a closer hit
	{
		// Better hit than the current hit
		ioHit.mFraction = min_fraction;
		ioHit.mSubShapeID2 = inSubShapeIDCreator.GetID();
		return true;
	}
	return false;
}

void ConvexHullShape::CastRay(const RayCast &inRay, const RayCastSettings &inRayCastSettings, const SubShapeIDCreator &inSubShapeIDCreator, CastRayCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	// Determine if ray hits the shape
	float min_fraction, max_fraction;
	if (CastRayHelper(inRay, min_fraction, max_fraction)
		&& min_fraction < ioCollector.GetEarlyOutFraction()) // Check if this is closer than the early out fraction
	{
		// Better hit than the current hit
		RayCastResult hit;
		hit.mBodyID = TransformedShape::sGetBodyID(ioCollector.GetContext());
		hit.mSubShapeID2 = inSubShapeIDCreator.GetID();

		// Check front side hit
		if (inRayCastSettings.mTreatConvexAsSolid || min_fraction > 0.0f)
		{
			hit.mFraction = min_fraction;
			ioCollector.AddHit(hit);
		}

		// Check back side hit
		if (inRayCastSettings.mBackFaceModeConvex == EBackFaceMode::CollideWithBackFaces
			&& max_fraction < ioCollector.GetEarlyOutFraction())
		{
			hit.mFraction = max_fraction;
			ioCollector.AddHit(hit);
		}
	}
}

void ConvexHullShape::CollidePoint(Vec4Arg inPoint, const SubShapeIDCreator &inSubShapeIDCreator, CollidePointCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	// Check if point is behind all planes
	for (const Plane &p : mPlanes)
		if (p.SignedDistance(inPoint) > 0.0f)
			return;

	// Point is inside
	ioCollector.AddHit({ TransformedShape::sGetBodyID(ioCollector.GetContext()), inSubShapeIDCreator.GetID() });
}

void ConvexHullShape::CollideSoftBodyVertices(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const CollideSoftBodyVertexIterator &inVertices, uint inNumVertices, int inCollidingShapeIndex) const
{
	// TODO(4D): CollideSoftBodyVertexIterator is still 3D. Stubbed until the soft-body subsystem is
	// migrated to 4D (SoftBody is cut for v1).
	(void)inCenterOfMassTransform; (void)inScale; (void)inVertices; (void)inNumVertices; (void)inCollidingShapeIndex;
}

#if 0 // TODO(4D): original 3D soft-body-vertex collision, kept for reference during the 4D port
void ConvexHullShape::CollideSoftBodyVerticesUnported(Mat44Arg inCenterOfMassTransform, Vec3Arg inScale, const CollideSoftBodyVertexIterator &inVertices, uint inNumVertices, int inCollidingShapeIndex) const
{
	Mat44 inverse_transform = inCenterOfMassTransform.InversedRotationTranslation();

	Vec3 inv_scale = inScale.Reciprocal();
	bool is_not_scaled = ScaleHelpers::IsNotScaled(inScale);
	float scale_flip = ScaleHelpers::IsInsideOut(inScale)? -1.0f : 1.0f;

	for (CollideSoftBodyVertexIterator v = inVertices, sbv_end = inVertices + inNumVertices; v != sbv_end; ++v)
		if (v.GetInvMass() > 0.0f)
		{
			Vec3 local_pos = inverse_transform * v.GetPosition();

			// Find most facing plane
			float max_distance = -FLT_MAX;
			Vec3 max_plane_normal = Vec3::sZero();
			uint max_plane_idx = 0;
			if (is_not_scaled)
			{
				// Without scale, it is trivial to calculate the distance to the hull
				for (const Plane &p : mPlanes)
				{
					float distance = p.SignedDistance(local_pos);
					if (distance > max_distance)
					{
						max_distance = distance;
						max_plane_normal = p.GetNormal();
						max_plane_idx = uint(&p - mPlanes.data());
					}
				}
			}
			else
			{
				// When there's scale we need to calculate the planes first
				for (uint i = 0; i < (uint)mPlanes.size(); ++i)
				{
					// Calculate plane normal and point by scaling the original plane
					Vec3 plane_normal = (inv_scale * mPlanes[i].GetNormal()).Normalized();
					Vec3 plane_point = inScale * mPoints[mVertexIdx[mFaces[i].mFirstVertex]].mPosition;

					float distance = plane_normal.Dot(local_pos - plane_point);
					if (distance > max_distance)
					{
						max_distance = distance;
						max_plane_normal = plane_normal;
						max_plane_idx = i;
					}
				}
			}

			// Project point onto that plane, in local space to the vertex
			Vec3 closest_point = -max_distance * max_plane_normal;

			// Check edges if we're outside the hull (when inside we know the closest face is also the closest point to the surface)
			bool is_outside = max_distance > 0.0f;
			if (is_outside)
			{
				// Loop over edges
				float closest_point_dist_sq = FLT_MAX;
				const Face &face = mFaces[max_plane_idx];
				for (const uint8 *v_start = &mVertexIdx[face.mFirstVertex], *v1 = v_start, *v_end = v_start + face.mNumVertices; v1 < v_end; ++v1)
				{
					// Find second point
					const uint8 *v2 = v1 + 1;
					if (v2 == v_end)
						v2 = v_start;

					// Get edge points
					Vec3 p1 = inScale * mPoints[*v1].mPosition;
					Vec3 p2 = inScale * mPoints[*v2].mPosition;

					// Check if the position is outside the edge (if not, the face will be closer)
					Vec3 edge_normal = (p2 - p1).Cross(max_plane_normal);
					if (scale_flip * edge_normal.Dot(local_pos - p1) > 0.0f)
					{
						// Get closest point on edge
						uint32 set;
						Vec3 closest = ClosestPoint::GetClosestPointOnLine(p1 - local_pos, p2 - local_pos, set);
						float distance_sq = closest.LengthSq();
						if (distance_sq < closest_point_dist_sq)
						{
							closest_point_dist_sq = distance_sq;
							closest_point = closest;
						}
					}
				}
			}

			// Check if this is the largest penetration
			Vec3 normal = -closest_point;
			float normal_length = normal.Length();
			float penetration = normal_length;
			if (is_outside)
				penetration = -penetration;
			else
				normal = -normal;
			if (v.UpdatePenetration(penetration))
			{
				// Calculate contact plane
				normal = normal_length > 1.0e-12f? normal / normal_length : max_plane_normal;
				Plane plane = Plane::sFromPointAndNormal(local_pos + closest_point, normal);

				// Store collision
				v.SetCollision(plane.GetTransformed(inCenterOfMassTransform), inCollidingShapeIndex);
			}
		}
}
#endif // 0

void ConvexHullShape::GetTrianglesStart(GetTrianglesContext &ioContext, const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale) const
{
	// TODO(4D): the GetTriangles pipeline is still 3D. Stubbed until the 4D (tetrahedra-based) boundary
	// tessellation is defined.
	(void)ioContext; (void)inBox; (void)inPositionCOM; (void)inRotation; (void)inScale;
}

int ConvexHullShape::GetTrianglesNext(GetTrianglesContext &ioContext, int inMaxTrianglesRequested, Float4 *outTriangleVertices, const PhysicsMaterial **outMaterials) const
{
	// TODO(4D): See GetTrianglesStart. Returns 0 to indicate no tessellation data is currently produced.
	(void)ioContext; (void)inMaxTrianglesRequested; (void)outTriangleVertices; (void)outMaterials;
	return 0;
}

void ConvexHullShape::SaveBinaryState(StreamOut &inStream) const
{
	ConvexShape::SaveBinaryState(inStream);

	inStream.Write(mCenterOfMass);
	inStream.Write(mInertia);
	inStream.Write(mLocalBounds.mMin);
	inStream.Write(mLocalBounds.mMax);
	inStream.Write(mPoints);
	inStream.Write(mFaces);
	inStream.Write(mPlanes);
	inStream.Write(mVertexIdx);
	inStream.Write(mConvexRadius);
	inStream.Write(mVolume);
	inStream.Write(mInnerRadius);
}

void ConvexHullShape::RestoreBinaryState(StreamIn &inStream)
{
	ConvexShape::RestoreBinaryState(inStream);

	inStream.Read(mCenterOfMass);
	inStream.Read(mInertia);
	inStream.Read(mLocalBounds.mMin);
	inStream.Read(mLocalBounds.mMax);
	inStream.Read(mPoints);
	inStream.Read(mFaces);
	inStream.Read(mPlanes);
	inStream.Read(mVertexIdx);
	inStream.Read(mConvexRadius);
	inStream.Read(mVolume);
	inStream.Read(mInnerRadius);
}

Shape::Stats ConvexHullShape::GetStats() const
{
	// Count number of triangles
	uint triangle_count = 0;
	for (const Face &f : mFaces)
		triangle_count += f.mNumVertices - 2;

	return Stats(
		sizeof(*this)
			+ mPoints.size() * sizeof(Point)
			+ mFaces.size() * sizeof(Face)
			+ mPlanes.size() * sizeof(Plane)
			+ mVertexIdx.size() * sizeof(uint8),
		triangle_count);
}

void ConvexHullShape::sRegister()
{
	ShapeFunctions &f = ShapeFunctions::sGet(EShapeSubType::ConvexHull);
	f.mConstruct = []() -> Shape * { return new ConvexHullShape; };
	f.mColor = Color::sGreen;
}

JPH_NAMESPACE_END
