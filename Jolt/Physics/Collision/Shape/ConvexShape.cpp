// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Collision/Shape/ConvexShape.h>
#include <Jolt/Physics/Collision/RayCast.h>
#include <Jolt/Physics/Collision/ShapeCast.h>
#include <Jolt/Physics/Collision/CollideShape.h>
#include <Jolt/Physics/Collision/CastResult.h>
#include <Jolt/Physics/Collision/CollidePointResult.h>
#include <Jolt/Physics/Collision/Shape/ScaleHelpers.h>
#include <Jolt/Physics/Collision/TransformedShape.h>
#include <Jolt/Physics/Collision/CollisionDispatch.h>
#include <Jolt/Physics/Collision/NarrowPhaseStats.h>
#include <Jolt/Physics/PhysicsSettings.h>
#include <Jolt/Core/StreamIn.h>
#include <Jolt/Core/StreamOut.h>
#include <Jolt/Geometry/EPAPenetrationDepth.h>
#include <Jolt/Geometry/OrientedBox.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>

JPH_NAMESPACE_BEGIN

JPH_IMPLEMENT_SERIALIZABLE_ABSTRACT(ConvexShapeSettings)
{
	JPH_ADD_BASE_CLASS(ConvexShapeSettings, ShapeSettings)

	JPH_ADD_ATTRIBUTE(ConvexShapeSettings, mDensity)
	JPH_ADD_ATTRIBUTE(ConvexShapeSettings, mMaterial)
}

// Narrow an RVec4 (world/relative position) to Vec4. The center-of-mass transforms passed to the
// narrow phase are relative to a base offset, so the values are always in single-precision range.
static JPH_INLINE Vec4 sNarrow(RVec4Arg inV)
{
	return Vec4(float(inV.GetX()), float(inV.GetY()), float(inV.GetZ()), float(inV.GetW()));
}

void ConvexShape::sCollideConvexVsConvex(const Shape *inShape1, const Shape *inShape2, Vec4Arg inScale1, Vec4Arg inScale2, RMat44Arg inCenterOfMassTransform1, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, [[maybe_unused]] const ShapeFilter &inShapeFilter)
{
	JPH_PROFILE_FUNCTION();

	// Get the shapes
	JPH_ASSERT(inShape1->GetType() == EShapeType::Convex);
	JPH_ASSERT(inShape2->GetType() == EShapeType::Convex);
	const ConvexShape *shape1 = static_cast<const ConvexShape *>(inShape1);
	const ConvexShape *shape2 = static_cast<const ConvexShape *>(inShape2);

	// Get transforms. The collision runs in shape 1's local space, so the relative transform's
	// rotation (Mat44) and translation (Vec4) drive the local-space support functions.
	RMat44 inverse_transform1 = inCenterOfMassTransform1.InversedRotationTranslation();
	RMat44 transform_2_to_1 = inverse_transform1 * inCenterOfMassTransform2;
	Mat44 rotation_2_to_1 = transform_2_to_1.GetRotation();
	Vec4 translation_2_to_1 = sNarrow(transform_2_to_1.GetTranslation());

	// Get bounding boxes
	float max_separation_distance = inCollideShapeSettings.mMaxSeparationDistance;
	AABox shape1_bbox = shape1->GetLocalBounds().Scaled(inScale1);
	shape1_bbox.ExpandBy(Vec4::sReplicate(max_separation_distance));
	AABox shape2_bbox = shape2->GetLocalBounds().Scaled(inScale2);

	// This is only a rejection filter; GJK below performs the actual shape test.
	// A full 4D OBB SAT tests 56 axes even for overlapping pairs. Use conservative
	// transformed bounds here instead, retaining the previous filter's tolerance.
	shape1_bbox.ExpandBy(Vec4::sReplicate(1.0e-6f));
	if (!shape2_bbox.Transformed(rotation_2_to_1, translation_2_to_1).Overlaps(shape1_bbox))
		return;

	// Note: As we don't remember the penetration axis from the last iteration, and it is likely that shape2 is pushed out of
	// collision relative to shape1 by comparing their COM's, we use that as an initial penetration axis: shape2.com - shape1.com
	// This has been seen to improve performance by approx. 1% over using a fixed axis like (1, 0, 0).
	Vec4 penetration_axis = translation_2_to_1;

	// Ensure that we do not pass in a near zero penetration axis
	if (penetration_axis.IsNearZero())
		penetration_axis = Vec4::sAxisX();

	Vec4 point1, point2;
	EPAPenetrationDepth pen_depth;
	EPAPenetrationDepth::EStatus status;

	// Scope to limit lifetime of SupportBuffer
	{
		// Create support function
		SupportBuffer buffer1_excl_cvx_radius, buffer2_excl_cvx_radius;
		const Support *shape1_excl_cvx_radius = shape1->GetSupportFunction(ConvexShape::ESupportMode::ExcludeConvexRadius, buffer1_excl_cvx_radius, inScale1);
		const Support *shape2_excl_cvx_radius = shape2->GetSupportFunction(ConvexShape::ESupportMode::ExcludeConvexRadius, buffer2_excl_cvx_radius, inScale2);

		// Transform shape 2 in the space of shape 1
		TransformedConvexObject transformed2_excl_cvx_radius(rotation_2_to_1, translation_2_to_1, *shape2_excl_cvx_radius);

		// Perform GJK step
		status = pen_depth.GetPenetrationDepthStepGJK(*shape1_excl_cvx_radius, shape1_excl_cvx_radius->GetConvexRadius() + max_separation_distance, transformed2_excl_cvx_radius, shape2_excl_cvx_radius->GetConvexRadius(), inCollideShapeSettings.mCollisionTolerance, penetration_axis, point1, point2);
	}

	// Check result of collision detection
	switch (status)
	{
	case EPAPenetrationDepth::EStatus::Colliding:
		break;

	case EPAPenetrationDepth::EStatus::NotColliding:
		return;

	case EPAPenetrationDepth::EStatus::Indeterminate:
		{
			// Need to run expensive EPA algorithm

			// TODO(4D): upstream inflates shape 1 by the max separation distance before running EPA
			// and subtracts it again afterwards (see the correction below). The margin guards against
			// EPA reporting "separated" for a pair GJK already called overlapping. In 4D it costs far
			// more than it does in 3D: the inflation is a sphere sweep, so it rounds every edge and
			// corner, and EPA's inscribed-polytope error on a *curved* surface falls off as V^(-2/3)
			// on S^3 against V^-1 on S^2. A tesseract resting on a tesseract then does not converge:
			// measured over a 1296-body pyramid, EPA averaged 62 support points and 330 hull cells
			// per call and ran into its 96-point cap on 57% of all calls - at which point it returns
			// the best facet it happened to find rather than the answer. Against the un-inflated pair
			// (two polytopes, which is what EPA is good at) the same scene averages 9.8 points and
			// 25.4 cells and never hits the cap. End to end the pyramid runs 10-25% faster depending
			// on phase and thread count, and with visibly less frame-time variance.
			//
			// This is safe to skip here: reaching EPA means the shapes *excluding* convex radius
			// already overlap, so the real shapes overlap too and the penetration is positive - the
			// margin was only ever a numerical guard, never part of the result. Speculative contacts
			// are unaffected because they come out of the GJK branch above. The margin is kept as a
			// fallback for the case it was written for: if EPA does report no collision, retry the
			// upstream way. (Measured 6.5% of EPA calls take the fallback; a 441-tesseract pyramid
			// run for 600 steps keeps every layer within 3.5 cm of its ideal height either way.)
			max_separation_distance = min(max_separation_distance, 1.0f);
			float epa_separation_distance = 0.0f;

			// Create support function
			SupportBuffer buffer1_incl_cvx_radius, buffer2_incl_cvx_radius;
			const Support *shape1_incl_cvx_radius = shape1->GetSupportFunction(ConvexShape::ESupportMode::IncludeConvexRadius, buffer1_incl_cvx_radius, inScale1);
			const Support *shape2_incl_cvx_radius = shape2->GetSupportFunction(ConvexShape::ESupportMode::IncludeConvexRadius, buffer2_incl_cvx_radius, inScale2);

			// Transform shape 2 in the space of shape 1
			TransformedConvexObject transformed2_incl_cvx_radius(rotation_2_to_1, translation_2_to_1, *shape2_incl_cvx_radius);

			// Perform EPA step against the shapes as they are
			if (!pen_depth.GetPenetrationDepthStepEPA(*shape1_incl_cvx_radius, transformed2_incl_cvx_radius, inCollideShapeSettings.mPenetrationTolerance, penetration_axis, point1, point2))
			{
				// Fall back to upstream's behaviour: inflate shape 1 by the separation distance so a
				// pair GJK called overlapping cannot come out of EPA separated
				if (max_separation_distance <= 0.0f)
					return;

				AddConvexRadius shape1_add_max_separation_distance(*shape1_incl_cvx_radius, max_separation_distance);

				// EPA leaves the axis in whatever state it failed in, so restart it the way the
				// first attempt was seeded
				penetration_axis = translation_2_to_1;
				if (penetration_axis.IsNearZero())
					penetration_axis = Vec4::sAxisX();

				if (!pen_depth.GetPenetrationDepthStepEPA(shape1_add_max_separation_distance, transformed2_incl_cvx_radius, inCollideShapeSettings.mPenetrationTolerance, penetration_axis, point1, point2))
					return;

				epa_separation_distance = max_separation_distance;
			}

			// The corrections below must undo whichever inflation was actually used
			max_separation_distance = epa_separation_distance;
			break;
		}
	}

	// Check if the penetration is bigger than the early out fraction
	float penetration_depth = (point2 - point1).Length() - max_separation_distance;
	if (-penetration_depth >= ioCollector.GetEarlyOutFraction())
		return;

	// Correct point1 for the added separation distance
	float penetration_axis_len = penetration_axis.Length();
	if (penetration_axis_len > 0.0f)
		point1 -= penetration_axis * (max_separation_distance / penetration_axis_len);

	// Convert to world space
	Vec4 point1_world = sNarrow(inCenterOfMassTransform1 * point1);
	Vec4 point2_world = sNarrow(inCenterOfMassTransform1 * point2);
	Vec4 penetration_axis_world = inCenterOfMassTransform1.Multiply3x3(penetration_axis);

	// Create collision result
	CollideShapeResult result(point1_world, point2_world, penetration_axis_world, penetration_depth, inSubShapeIDCreator1.GetID(), inSubShapeIDCreator2.GetID(), TransformedShape::sGetBodyID(ioCollector.GetContext()));

	// Gather faces
	if (inCollideShapeSettings.mCollectFacesMode == ECollectFacesMode::CollectFaces)
	{
		// Get supporting face of shape 1
		shape1->GetSupportingFace(SubShapeID(), -penetration_axis, inScale1, inCenterOfMassTransform1, result.mShape1Face);

		// Get supporting face of shape 2 (transform the axis into shape 2's local frame: R^T * axis)
		shape2->GetSupportingFace(SubShapeID(), Vec4(rotation_2_to_1.Transposed() * Lane4(penetration_axis)), inScale2, inCenterOfMassTransform2, result.mShape2Face);
	}

	// Notify the collector
	JPH_IF_TRACK_NARROWPHASE_STATS(TrackNarrowPhaseCollector track;)
	ioCollector.AddHit(result);
}

bool ConvexShape::CastRay(const RayCast &inRay, const SubShapeIDCreator &inSubShapeIDCreator, RayCastResult &ioHit) const
{
	// Note: This is a fallback routine, most convex shapes should implement a more performant version!

	JPH_PROFILE_FUNCTION();

	// Create support function
	SupportBuffer buffer;
	const Support *support = GetSupportFunction(ConvexShape::ESupportMode::IncludeConvexRadius, buffer, Vec4::sOne());

	// Cast ray
	GJKClosestPoint gjk;
	if (gjk.CastRay(inRay.mOrigin, inRay.mDirection, cDefaultCollisionTolerance, *support, ioHit.mFraction))
	{
		ioHit.mSubShapeID2 = inSubShapeIDCreator.GetID();
		return true;
	}

	return false;
}

void ConvexShape::CastRay(const RayCast &inRay, const RayCastSettings &inRayCastSettings, const SubShapeIDCreator &inSubShapeIDCreator, CastRayCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Note: This is a fallback routine, most convex shapes should implement a more performant version!

	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	// First do a normal raycast, limited to the early out fraction
	RayCastResult hit;
	hit.mFraction = ioCollector.GetEarlyOutFraction();
	if (CastRay(inRay, inSubShapeIDCreator, hit))
	{
		// Check front side
		if (inRayCastSettings.mTreatConvexAsSolid || hit.mFraction > 0.0f)
		{
			hit.mBodyID = TransformedShape::sGetBodyID(ioCollector.GetContext());
			ioCollector.AddHit(hit);
		}

		// Check if we want back facing hits and the collector still accepts additional hits
		if (inRayCastSettings.mBackFaceModeConvex == EBackFaceMode::CollideWithBackFaces && !ioCollector.ShouldEarlyOut())
		{
			// Invert the ray, going from the early out fraction back to the fraction where we found our forward hit
			float start_fraction = min(1.0f, ioCollector.GetEarlyOutFraction());
			float delta_fraction = hit.mFraction - start_fraction;
			if (delta_fraction < 0.0f)
			{
				RayCast inverted_ray { inRay.mOrigin + start_fraction * inRay.mDirection, delta_fraction * inRay.mDirection };

				// Cast another ray
				RayCastResult inverted_hit;
				inverted_hit.mFraction = 1.0f;
				if (CastRay(inverted_ray, inSubShapeIDCreator, inverted_hit)
					&& inverted_hit.mFraction > 0.0f) // Ignore hits with fraction 0, this means the ray ends inside the object and we don't want to report it as a back facing hit
				{
					// Invert fraction and rescale it to the fraction of the original ray
					inverted_hit.mFraction = hit.mFraction + (inverted_hit.mFraction - 1.0f) * delta_fraction;
					inverted_hit.mBodyID = TransformedShape::sGetBodyID(ioCollector.GetContext());
					ioCollector.AddHit(inverted_hit);
				}
			}
		}
	}
}

void ConvexShape::CollidePoint(Vec4Arg inPoint, const SubShapeIDCreator &inSubShapeIDCreator, CollidePointCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	// First test bounding box
	if (GetLocalBounds().Contains(inPoint))
	{
		// Create support function
		SupportBuffer buffer;
		const Support *support = GetSupportFunction(ConvexShape::ESupportMode::IncludeConvexRadius, buffer, Vec4::sOne());

		// Create support function for point
		PointConvexSupport point { inPoint };

		// Test intersection
		GJKClosestPoint gjk;
		Vec4 v = inPoint;
		if (gjk.Intersects(*support, point, cDefaultCollisionTolerance, v))
			ioCollector.AddHit({ TransformedShape::sGetBodyID(ioCollector.GetContext()), inSubShapeIDCreator.GetID() });
	}
}

void ConvexShape::sCastConvexVsConvex(const ShapeCast &inShapeCast, const ShapeCastSettings &inShapeCastSettings, const Shape *inShape, Vec4Arg inScale, [[maybe_unused]] const ShapeFilter &inShapeFilter, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, CastShapeCollector &ioCollector)
{
	JPH_PROFILE_FUNCTION();

	// Only supported for convex shapes
	JPH_ASSERT(inShapeCast.mShape->GetType() == EShapeType::Convex);
	const ConvexShape *cast_shape = static_cast<const ConvexShape *>(inShapeCast.mShape);

	JPH_ASSERT(inShape->GetType() == EShapeType::Convex);
	const ConvexShape *shape = static_cast<const ConvexShape *>(inShape);

	// Determine if we want to use the actual shape or a shrunken shape with convex radius
	ConvexShape::ESupportMode support_mode = inShapeCastSettings.mUseShrunkenShapeAndConvexRadius? ConvexShape::ESupportMode::ExcludeConvexRadius : ConvexShape::ESupportMode::Default;

	// Create support function for shape to cast
	SupportBuffer cast_buffer;
	const Support *cast_support = cast_shape->GetSupportFunction(support_mode, cast_buffer, inShapeCast.mScale);

	// Create support function for target shape
	SupportBuffer target_buffer;
	const Support *target_support = shape->GetSupportFunction(support_mode, target_buffer, inScale);

	// Do a raycast against the result. The cast runs in shape 2's local space, so the cast start
	// transform's rotation (Mat44) and translation (Vec4) are passed to the geometry separately.
	EPAPenetrationDepth epa;
	float fraction = ioCollector.GetEarlyOutFraction();
	Vec4 contact_point_a, contact_point_b, contact_normal;
	if (epa.CastShape(inShapeCast.mCenterOfMassStart.GetRotation(), sNarrow(inShapeCast.mCenterOfMassStart.GetTranslation()), inShapeCast.mDirection, inShapeCastSettings.mCollisionTolerance, inShapeCastSettings.mPenetrationTolerance, *cast_support, *target_support, cast_support->GetConvexRadius(), target_support->GetConvexRadius(), inShapeCastSettings.mReturnDeepestPoint, fraction, contact_point_a, contact_point_b, contact_normal)
		&& (inShapeCastSettings.mBackFaceModeConvex == EBackFaceMode::CollideWithBackFaces
			|| contact_normal.Dot(inShapeCast.mDirection) > 0.0f)) // Test if backfacing
	{
		// Convert to world space
		Vec4 contact_point_a_world = sNarrow(inCenterOfMassTransform2 * contact_point_a);
		Vec4 contact_point_b_world = sNarrow(inCenterOfMassTransform2 * contact_point_b);
		Vec4 contact_normal_world = inCenterOfMassTransform2.Multiply3x3(contact_normal);

		ShapeCastResult result(fraction, contact_point_a_world, contact_point_b_world, contact_normal_world, false, inSubShapeIDCreator1.GetID(), inSubShapeIDCreator2.GetID(), TransformedShape::sGetBodyID(ioCollector.GetContext()));

		// Early out if this hit is deeper than the collector's early out value
		if (fraction == 0.0f && -result.mPenetrationDepth >= ioCollector.GetEarlyOutFraction())
			return;

		// Gather faces
		if (inShapeCastSettings.mCollectFacesMode == ECollectFacesMode::CollectFaces)
		{
			// Get supporting face of shape 1 (at the point of contact along the sweep)
			RMat44 transform_1_to_2 = inShapeCast.mCenterOfMassStart;
			transform_1_to_2.SetTranslation(transform_1_to_2.GetTranslation() + RVec4(fraction * inShapeCast.mDirection));
			cast_shape->GetSupportingFace(SubShapeID(), transform_1_to_2.Multiply3x3Transposed(-contact_normal), inShapeCast.mScale, inCenterOfMassTransform2 * transform_1_to_2, result.mShape1Face);

			// Get supporting face of shape 2
			shape->GetSupportingFace(SubShapeID(), contact_normal, inScale, inCenterOfMassTransform2, result.mShape2Face);
		}

		JPH_IF_TRACK_NARROWPHASE_STATS(TrackNarrowPhaseCollector track;)
		ioCollector.AddHit(result);
	}
}

void ConvexShape::GetTrianglesStart(GetTrianglesContext &ioContext, const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale) const
{
	// TODO(4D): the convex boundary tessellation returned sUnitSphereTriangles (a 3D triangle
	// mesh) via GetTrianglesContextVertexList. In 4D the boundary is a 3-manifold tessellated by
	// tetrahedra; the GetTrianglesContext pipeline still assumes triangles. Stubbed until that is
	// defined.
	(void)ioContext; (void)inBox; (void)inPositionCOM; (void)inRotation; (void)inScale;
}

int ConvexShape::GetTrianglesNext(GetTrianglesContext &ioContext, int inMaxTrianglesRequested, Float4 *outTriangleVertices, const PhysicsMaterial **outMaterials) const
{
	// TODO(4D): see GetTrianglesStart.
	(void)ioContext; (void)inMaxTrianglesRequested; (void)outTriangleVertices; (void)outMaterials;
	return 0;
}

void ConvexShape::GetSubmergedVolume(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const Plane &inSurface, float &outTotalVolume, float &outSubmergedVolume, Vec4 &outCenterOfBuoyancy, [[maybe_unused]] RVec4Arg inBaseOffset) const
{
	// 4D hypervolume of the bounding box (2*extent per axis)
	Vec4 extent = GetLocalBounds().GetExtent() * inScale.Abs();
	outTotalVolume = 16.0f * extent.GetX() * extent.GetY() * extent.GetZ() * extent.GetW();

	// TODO(4D): the submerged-hypervolume of a convex 4-polytope cut by a hyperplane needs a 4D
	// polyhedron clip (PolyhedronSubmergedVolumeCalculator is the 3D box version). Approximate for
	// now with the fully-above / fully-below cases and leave the intersecting case at zero.
	(void)inSurface;
	outSubmergedVolume = 0.0f;
	outCenterOfBuoyancy = sNarrow(inCenterOfMassTransform.GetTranslation());
}

#ifdef JPH_DEBUG_RENDERER
void ConvexShape::DrawGetSupportFunction(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec3Arg inScale, ColorArg inColor, bool inDrawSupportDirection) const
{
	// Get the support function with convex radius
	SupportBuffer buffer;
	const Support *support = GetSupportFunction(ESupportMode::ExcludeConvexRadius, buffer, inScale);
	AddConvexRadius add_convex(*support, support->GetConvexRadius());

	// Draw the shape
	DebugRenderer::GeometryRef geometry = inRenderer->CreateTriangleGeometryForConvex([&add_convex](Vec3Arg inDirection) { return add_convex.GetSupport(inDirection); });
	AABox bounds = geometry->mBounds.Transformed(inCenterOfMassTransform);
	float lod_scale_sq = geometry->mBounds.GetExtent().LengthSq();
	inRenderer->DrawGeometry(inCenterOfMassTransform, bounds, lod_scale_sq, inColor, geometry);

	if (inDrawSupportDirection)
	{
		// Iterate on all directions and draw the support point and an arrow in the direction that was sampled to test if the support points make sense
		for (Vec3 v : Vec3::sUnitSphere)
		{
			Vec3 direction = 0.05f * v;
			Vec3 pos = add_convex.GetSupport(direction);
			RVec3 from = inCenterOfMassTransform * pos;
			RVec3 to = inCenterOfMassTransform * (pos + direction);
			inRenderer->DrawMarker(from, Color::sWhite, 0.001f);
			inRenderer->DrawArrow(from, to, Color::sWhite, 0.001f);
		}
	}
}

void ConvexShape::DrawGetSupportingFace(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec3Arg inScale) const
{
	// Sample directions and map which faces belong to which directions
	using FaceToDirection = UnorderedMap<SupportingFace, Array<Vec3>>;
	FaceToDirection faces;
	for (Vec3 v : Vec3::sUnitSphere)
	{
		Vec3 direction = 0.05f * v;

		SupportingFace face;
		GetSupportingFace(SubShapeID(), direction, inScale, Mat44::sIdentity(), face);

		if (!face.empty())
		{
			JPH_ASSERT(face.size() >= 2, "The GetSupportingFace function should either return nothing or at least an edge");
			faces[face].push_back(direction);
		}
	}

	// Draw each face in a unique color and draw corresponding directions
	int color_it = 0;
	for (FaceToDirection::value_type &ftd : faces)
	{
		Color color = Color::sGetDistinctColor(color_it++);

		// Create copy of face (key in map is read only)
		SupportingFace face = ftd.first;

		// Displace the face a little bit forward so it is easier to see
		Vec3 normal = face.size() >= 3? (face[2] - face[1]).Cross(face[0] - face[1]).NormalizedOr(Vec3::sZero()) : Vec3::sZero();
		Vec3 displacement = 0.001f * normal;

		// Transform face to world space and calculate center of mass
		Vec3 com_ls = Vec3::sZero();
		for (Vec3 &v : face)
		{
			v = inCenterOfMassTransform.Multiply3x3(v + displacement);
			com_ls += v;
		}
		RVec3 com = inCenterOfMassTransform.GetTranslation() + com_ls / (float)face.size();

		// Draw the polygon and directions
		inRenderer->DrawWirePolygon(RMat44::sTranslation(inCenterOfMassTransform.GetTranslation()), face, color, face.size() >= 3? 0.001f : 0.0f);
		if (face.size() >= 3)
			inRenderer->DrawArrow(com, com + inCenterOfMassTransform.Multiply3x3(normal), color, 0.01f);
		for (Vec3 &v : ftd.second)
			inRenderer->DrawArrow(com, com + inCenterOfMassTransform.Multiply3x3(-v), color, 0.001f);
	}
}
#endif // JPH_DEBUG_RENDERER

void ConvexShape::SaveBinaryState(StreamOut &inStream) const
{
	Shape::SaveBinaryState(inStream);

	inStream.Write(mDensity);
}

void ConvexShape::RestoreBinaryState(StreamIn &inStream)
{
	Shape::RestoreBinaryState(inStream);

	inStream.Read(mDensity);
}

void ConvexShape::SaveMaterialState(PhysicsMaterialList &outMaterials) const
{
	outMaterials.clear();
	outMaterials.push_back(mMaterial);
}

void ConvexShape::RestoreMaterialState(const PhysicsMaterialRefC *inMaterials, uint inNumMaterials)
{
	JPH_ASSERT(inNumMaterials == 1);
	mMaterial = inMaterials[0];
}

void ConvexShape::sRegister()
{
	for (EShapeSubType s1 : sConvexSubShapeTypes)
		for (EShapeSubType s2 : sConvexSubShapeTypes)
		{
			CollisionDispatch::sRegisterCollideShape(s1, s2, sCollideConvexVsConvex);
			CollisionDispatch::sRegisterCastShape(s1, s2, sCastConvexVsConvex);
		}
}

JPH_NAMESPACE_END
