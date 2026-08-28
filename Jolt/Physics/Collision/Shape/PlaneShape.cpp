// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2024 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Collision/Shape/PlaneShape.h>
#include <Jolt/Physics/Collision/Shape/ConvexShape.h>
#include <Jolt/Physics/Collision/Shape/ScaleHelpers.h>
#include <Jolt/Physics/Collision/RayCast.h>
#include <Jolt/Physics/Collision/ShapeCast.h>
#include <Jolt/Physics/Collision/ShapeFilter.h>
#include <Jolt/Physics/Collision/CastResult.h>
#include <Jolt/Physics/Collision/CollisionDispatch.h>
#include <Jolt/Physics/Collision/TransformedShape.h>
#include <Jolt/Physics/Collision/CollidePointResult.h>
#include <Jolt/Physics/Collision/CollideSoftBodyVertexIterator.h>
#include <Jolt/Core/Profiler.h>
#include <Jolt/Core/StreamIn.h>
#include <Jolt/Core/StreamOut.h>
#include <Jolt/Geometry/Plane.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_DEBUG_RENDERER

JPH_NAMESPACE_BEGIN

JPH_IMPLEMENT_SERIALIZABLE_VIRTUAL(PlaneShapeSettings)
{
	JPH_ADD_BASE_CLASS(PlaneShapeSettings, ShapeSettings)

	JPH_ADD_ATTRIBUTE(PlaneShapeSettings, mPlane)
	JPH_ADD_ATTRIBUTE(PlaneShapeSettings, mMaterial)
	JPH_ADD_ATTRIBUTE(PlaneShapeSettings, mHalfExtent)
}

ShapeSettings::ShapeResult PlaneShapeSettings::Create() const
{
	if (mCachedResult.IsEmpty())
		Ref<Shape> shape = new PlaneShape(*this, mCachedResult);
	return mCachedResult;
}

// Build a right-handed orthonormal basis {outPerp1, outPerp2, outPerp3} spanning the hyperplane
// perpendicular to inNormal (in 4D the plane's tangent space is 3-dimensional).
// TODO(4D): the orthonormal 3-frame of a hyperplane normal now lives on Vec4 itself
// (Vec4::GetNormalizedPerpendicularBasis) because contact friction needs the same 3 tangents.

void PlaneShape::GetVertices(Vec4 *outVertices) const
{
	// Create orthonormal basis for the plane's 3D tangent space
	Vec4 normal = mPlane.GetNormal();
	Vec4 perp1, perp2, perp3;
	normal.GetNormalizedPerpendicularBasis(perp1, perp2, perp3);

	// Scale basis
	perp1 *= mHalfExtent;
	perp2 *= mHalfExtent;
	perp3 *= mHalfExtent;

	// Calculate the 8 corners of the 3-cube patch
	Vec4 point = -normal * mPlane.GetConstant();
	for (int i = 0; i < 8; ++i)
		outVertices[i] = point + ((i & 1)? perp1 : -perp1) + ((i & 2)? perp2 : -perp2) + ((i & 4)? perp3 : -perp3);
}

void PlaneShape::CalculateLocalBounds()
{
	// Get the vertices of the plane
	Vec4 vertices[8];
	GetVertices(vertices);

	// Encapsulate the vertices and a point mHalfExtent behind the plane
	mLocalBounds = AABox();
	Vec4 normal = mPlane.GetNormal();
	for (const Vec4 &v : vertices)
	{
		mLocalBounds.Encapsulate(v);
		mLocalBounds.Encapsulate(v - mHalfExtent * normal);
	}
}

PlaneShape::PlaneShape(const PlaneShapeSettings &inSettings, ShapeResult &outResult) :
	Shape(EShapeType::Plane, EShapeSubType::Plane, inSettings, outResult),
	mPlane(inSettings.mPlane),
	mMaterial(inSettings.mMaterial),
	mHalfExtent(inSettings.mHalfExtent)
{
	if (!mPlane.GetNormal().IsNormalized())
	{
		outResult.SetError("Plane normal needs to be normalized!");
		return;
	}

	CalculateLocalBounds();

	outResult.Set(this);
}

MassProperties PlaneShape::GetMassProperties() const
{
	// Object should always be static, return default mass properties
	return MassProperties();
}

void PlaneShape::GetSupportingFace(const SubShapeID &inSubShapeID, Vec4Arg inDirection, Vec4Arg inScale, RMat44Arg inCenterOfMassTransform, SupportingFace &outVertices) const
{
	// Get the vertices of the plane (a 3-cube patch)
	Vec4 vertices[8];
	GetVertices(vertices);

	// TODO(4D): winding / face ordering for the 3D contact polytope is part of the deferred 4D
	// contact-manifold pipeline (conceptual gap 1). For now we emit all 8 corners unordered.

	// Transform them to world space
	outVertices.clear();
	RMat44 com = inCenterOfMassTransform.PreScaled(inScale);
	for (const Vec4 &v : vertices)
		outVertices.push_back(Vec4(com * v));
}

#ifdef JPH_DEBUG_RENDERER
void PlaneShape::Draw(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inUseMaterialColors, bool inDrawWireframe) const
{
	// TODO(4D): Debug rendering is deferred (v1 is headless). A 4D plane patch renders as a 3-cube.
	(void)inRenderer; (void)inCenterOfMassTransform; (void)inScale; (void)inColor; (void)inUseMaterialColors; (void)inDrawWireframe;
}
#endif // JPH_DEBUG_RENDERER

bool PlaneShape::CastRay(const RayCast &inRay, const SubShapeIDCreator &inSubShapeIDCreator, RayCastResult &ioHit) const
{
	JPH_PROFILE_FUNCTION();

	// Test starting inside of negative half space
	float distance = mPlane.SignedDistance(inRay.mOrigin);
	if (distance <= 0.0f)
	{
		ioHit.mFraction = 0.0f;
		ioHit.mSubShapeID2 = inSubShapeIDCreator.GetID();
		return true;
	}

	// Test ray parallel to plane
	float dot = inRay.mDirection.Dot(mPlane.GetNormal());
	if (dot == 0.0f)
		return false;

	// Calculate hit fraction
	float fraction = -distance / dot;
	if (fraction >= 0.0f && fraction < ioHit.mFraction)
	{
		ioHit.mFraction = fraction;
		ioHit.mSubShapeID2 = inSubShapeIDCreator.GetID();
		return true;
	}

	return false;
}

void PlaneShape::CastRay(const RayCast &inRay, const RayCastSettings &inRayCastSettings, const SubShapeIDCreator &inSubShapeIDCreator, CastRayCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	JPH_PROFILE_FUNCTION();

	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	// Inside solid half space?
	float distance = mPlane.SignedDistance(inRay.mOrigin);
	if (inRayCastSettings.mTreatConvexAsSolid
		&& distance <= 0.0f // Inside plane
		&& ioCollector.GetEarlyOutFraction() > 0.0f) // Willing to accept hits at fraction 0
	{
		// Hit at fraction 0
		RayCastResult hit;
		hit.mBodyID = TransformedShape::sGetBodyID(ioCollector.GetContext());
		hit.mFraction = 0.0f;
		hit.mSubShapeID2 = inSubShapeIDCreator.GetID();
		ioCollector.AddHit(hit);
	}

	float dot = inRay.mDirection.Dot(mPlane.GetNormal());
	if (dot != 0.0f // Parallel ray will not hit plane
		&& (inRayCastSettings.mBackFaceModeConvex == EBackFaceMode::CollideWithBackFaces || dot < 0.0f)) // Back face culling
	{
		// Calculate hit with plane
		float fraction = -distance / dot;
		if (fraction >= 0.0f && fraction < ioCollector.GetEarlyOutFraction())
		{
			RayCastResult hit;
			hit.mBodyID = TransformedShape::sGetBodyID(ioCollector.GetContext());
			hit.mFraction = fraction;
			hit.mSubShapeID2 = inSubShapeIDCreator.GetID();
			ioCollector.AddHit(hit);
		}
	}
}

void PlaneShape::CollidePoint(Vec4Arg inPoint, const SubShapeIDCreator &inSubShapeIDCreator, CollidePointCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	JPH_PROFILE_FUNCTION();

	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	// Check if the point is inside the plane
	if (mPlane.SignedDistance(inPoint) < 0.0f)
		ioCollector.AddHit({ TransformedShape::sGetBodyID(ioCollector.GetContext()), inSubShapeIDCreator.GetID() });
}

void PlaneShape::CollideSoftBodyVertices(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const CollideSoftBodyVertexIterator &inVertices, uint inNumVertices, int inCollidingShapeIndex) const
{
	// TODO(4D): CollideSoftBodyVertexIterator still exposes 3D positions. Stubbed until the soft-body
	// subsystem is migrated to 4D (SoftBody is cut for v1).
	(void)inCenterOfMassTransform; (void)inScale; (void)inVertices; (void)inNumVertices; (void)inCollidingShapeIndex;
}

// This is a version of GetSupportingFace that returns a face large enough to cover the shape we're colliding with but
// not as large as the regular GetSupportingFace to avoid numerical precision issues (in 4D the face is a 3-cube).
inline static void sGetSupportingFace(const ConvexShape *inShape, Vec4Arg inShapeCOM, const Plane &inPlane, RMat44Arg inPlaneToWorld, ConvexShape::SupportingFace &outPlaneFace)
{
	// Project COM of shape onto plane
	Plane world_plane = inPlane.GetTransformed(inPlaneToWorld.GetRotation(), Vec4(inPlaneToWorld.GetTranslation()));
	Vec4 center = world_plane.ProjectPointOnPlane(inShapeCOM);

	// Create orthonormal basis for the plane
	Vec4 normal = world_plane.GetNormal();
	Vec4 perp1, perp2, perp3;
	normal.GetNormalizedPerpendicularBasis(perp1, perp2, perp3);

	// Base the size of the face on the bounding box of the shape, ensuring that it is large enough to cover the entire shape
	float size = inShape->GetLocalBounds().GetSize().Length();
	perp1 *= size;
	perp2 *= size;
	perp3 *= size;

	// Emit the 8 corners of the 3-cube face
	outPlaneFace.resize(8);
	for (int i = 0; i < 8; ++i)
		outPlaneFace[i] = center + ((i & 1)? perp1 : -perp1) + ((i & 2)? perp2 : -perp2) + ((i & 4)? perp3 : -perp3);
}

void PlaneShape::sCastConvexVsPlane(const ShapeCast &inShapeCast, const ShapeCastSettings &inShapeCastSettings, const Shape *inShape, Vec4Arg inScale, [[maybe_unused]] const ShapeFilter &inShapeFilter, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, CastShapeCollector &ioCollector)
{
	JPH_PROFILE_FUNCTION();

	// Get the shapes
	JPH_ASSERT(inShapeCast.mShape->GetType() == EShapeType::Convex);
	JPH_ASSERT(inShape->GetType() == EShapeType::Plane);
	const ConvexShape *convex_shape = static_cast<const ConvexShape *>(inShapeCast.mShape);
	const PlaneShape *plane_shape = static_cast<const PlaneShape *>(inShape);

	// Shape cast is provided relative to COM of inShape, so all we need to do is transform our plane with inScale
	Plane plane = plane_shape->mPlane.Scaled(inScale);
	Vec4 normal = plane.GetNormal();

	// Get support function
	ConvexShape::SupportBuffer shape1_support_buffer;
	const ConvexShape::Support *shape1_support = convex_shape->GetSupportFunction(ConvexShape::ESupportMode::Default, shape1_support_buffer, inShapeCast.mScale);

	// Get the support point of the convex shape in the opposite direction of the plane normal in our local space
	Vec4 normal_in_convex_shape_space = inShapeCast.mCenterOfMassStart.Multiply3x3Transposed(normal);
	Vec4 support_point = Vec4(inShapeCast.mCenterOfMassStart * shape1_support->GetSupport(-normal_in_convex_shape_space));
	float signed_distance = plane.SignedDistance(support_point);
	float convex_radius = shape1_support->GetConvexRadius();
	float penetration_depth = -signed_distance + convex_radius;
	float dot = inShapeCast.mDirection.Dot(normal);

	// Collision output
	RMat44 com_hit;
	Vec4 point1, point2;
	float fraction;

	// Do we start in collision?
	if (penetration_depth > 0.0f)
	{
		// Back face culling?
		if (inShapeCastSettings.mBackFaceModeConvex == EBackFaceMode::IgnoreBackFaces && dot > 0.0f)
			return;

		// Shallower hit?
		if (penetration_depth <= -ioCollector.GetEarlyOutFraction())
			return;

		// We're hitting at fraction 0
		fraction = 0.0f;

		// Get contact point
		com_hit = inCenterOfMassTransform2;
		point1 = Vec4(inCenterOfMassTransform2 * (support_point - normal * convex_radius));
		point2 = Vec4(inCenterOfMassTransform2 * (support_point - normal * signed_distance));
	}
	else if (dot < 0.0f) // Moving towards the plane?
	{
		// Calculate hit fraction
		fraction = penetration_depth / dot;
		JPH_ASSERT(fraction >= 0.0f);

		// Further than early out fraction?
		if (fraction >= ioCollector.GetEarlyOutFraction())
			return;

		// Get contact point
		com_hit = inCenterOfMassTransform2.PostTranslated(fraction * inShapeCast.mDirection);
		point1 = point2 = Vec4(com_hit * (support_point - normal * convex_radius));
	}
	else
	{
		// Moving away from the plane
		return;
	}

	// Create cast result
	Vec4 penetration_axis_world = com_hit.Multiply3x3(-normal);
	bool back_facing = dot > 0.0f;
	ShapeCastResult result(fraction, point1, point2, penetration_axis_world, back_facing, inSubShapeIDCreator1.GetID(), inSubShapeIDCreator2.GetID(), TransformedShape::sGetBodyID(ioCollector.GetContext()));

	// Gather faces
	if (inShapeCastSettings.mCollectFacesMode == ECollectFacesMode::CollectFaces)
	{
		// Get supporting face of convex shape
		RMat44 shape_to_world = com_hit * inShapeCast.mCenterOfMassStart;
		convex_shape->GetSupportingFace(SubShapeID(), normal_in_convex_shape_space, inShapeCast.mScale, shape_to_world, result.mShape1Face);

		// Get supporting face of plane
		if (!result.mShape1Face.empty())
			sGetSupportingFace(convex_shape, Vec4(shape_to_world.GetTranslation()), plane, inCenterOfMassTransform2, result.mShape2Face);
	}

	// Notify the collector
	JPH_IF_TRACK_NARROWPHASE_STATS(TrackNarrowPhaseCollector track;)
	ioCollector.AddHit(result);
}

void PlaneShape::GetTrianglesStart(GetTrianglesContext &ioContext, const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale) const
{
	// TODO(4D): The GetTriangles pipeline is still 3D. Stubbed until the 4D (tetrahedra-based)
	// boundary tessellation is defined.
	(void)ioContext; (void)inBox; (void)inPositionCOM; (void)inRotation; (void)inScale;
}

int PlaneShape::GetTrianglesNext(GetTrianglesContext &ioContext, int inMaxTrianglesRequested, Float4 *outTriangleVertices, const PhysicsMaterial **outMaterials) const
{
	// TODO(4D): See GetTrianglesStart. Returns 0 to indicate no tessellation data is currently produced.
	(void)ioContext; (void)inMaxTrianglesRequested; (void)outTriangleVertices; (void)outMaterials;
	return 0;
}

void PlaneShape::sCollideConvexVsPlane(const Shape *inShape1, const Shape *inShape2, Vec4Arg inScale1, Vec4Arg inScale2, RMat44Arg inCenterOfMassTransform1, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, [[maybe_unused]] const ShapeFilter &inShapeFilter)
{
	JPH_PROFILE_FUNCTION();

	// Get the shapes
	JPH_ASSERT(inShape1->GetType() == EShapeType::Convex);
	JPH_ASSERT(inShape2->GetType() == EShapeType::Plane);
	const ConvexShape *shape1 = static_cast<const ConvexShape *>(inShape1);
	const PlaneShape *shape2 = static_cast<const PlaneShape *>(inShape2);

	// Transform the plane to the space of the convex shape
	Plane scaled_plane = shape2->mPlane.Scaled(inScale2);
	RMat44 transform_2_to_1 = inCenterOfMassTransform1.InversedRotationTranslation() * inCenterOfMassTransform2;
	Plane plane = scaled_plane.GetTransformed(transform_2_to_1.GetRotation(), Vec4(transform_2_to_1.GetTranslation()));
	Vec4 normal = plane.GetNormal();

	// Get support function
	ConvexShape::SupportBuffer shape1_support_buffer;
	const ConvexShape::Support *shape1_support = shape1->GetSupportFunction(ConvexShape::ESupportMode::Default, shape1_support_buffer, inScale1);

	// Get the support point of the convex shape in the opposite direction of the plane normal
	Vec4 support_point = shape1_support->GetSupport(-normal);
	float signed_distance = plane.SignedDistance(support_point);
	float convex_radius = shape1_support->GetConvexRadius();
	float penetration_depth = -signed_distance + convex_radius;
	if (penetration_depth > -inCollideShapeSettings.mMaxSeparationDistance)
	{
		// Get contact point
		Vec4 point1 = Vec4(inCenterOfMassTransform1 * (support_point - normal * convex_radius));
		Vec4 point2 = Vec4(inCenterOfMassTransform1 * (support_point - normal * signed_distance));
		Vec4 penetration_axis_world = inCenterOfMassTransform1.Multiply3x3(-normal);

		// Create collision result
		CollideShapeResult result(point1, point2, penetration_axis_world, penetration_depth, inSubShapeIDCreator1.GetID(), inSubShapeIDCreator2.GetID(), TransformedShape::sGetBodyID(ioCollector.GetContext()));

		// Gather faces
		if (inCollideShapeSettings.mCollectFacesMode == ECollectFacesMode::CollectFaces)
		{
			// Get supporting face of shape 1
			shape1->GetSupportingFace(SubShapeID(), normal, inScale1, inCenterOfMassTransform1, result.mShape1Face);

			// Get supporting face of shape 2
			if (!result.mShape1Face.empty())
				sGetSupportingFace(shape1, Vec4(inCenterOfMassTransform1.GetTranslation()), scaled_plane, inCenterOfMassTransform2, result.mShape2Face);
		}

		// Notify the collector
		JPH_IF_TRACK_NARROWPHASE_STATS(TrackNarrowPhaseCollector track;)
		ioCollector.AddHit(result);
	}
}

void PlaneShape::SaveBinaryState(StreamOut &inStream) const
{
	Shape::SaveBinaryState(inStream);

	inStream.Write(mPlane);
	inStream.Write(mHalfExtent);
}

void PlaneShape::RestoreBinaryState(StreamIn &inStream)
{
	Shape::RestoreBinaryState(inStream);

	inStream.Read(mPlane);
	inStream.Read(mHalfExtent);

	CalculateLocalBounds();
}

void PlaneShape::SaveMaterialState(PhysicsMaterialList &outMaterials) const
{
	outMaterials = { mMaterial };
}

void PlaneShape::RestoreMaterialState(const PhysicsMaterialRefC *inMaterials, uint inNumMaterials)
{
	JPH_ASSERT(inNumMaterials == 1);
	mMaterial = inMaterials[0];
}

void PlaneShape::sRegister()
{
	ShapeFunctions &f = ShapeFunctions::sGet(EShapeSubType::Plane);
	f.mConstruct = []() -> Shape * { return new PlaneShape; };
	f.mColor = Color::sDarkRed;

	for (EShapeSubType s : sConvexSubShapeTypes)
	{
		CollisionDispatch::sRegisterCollideShape(s, EShapeSubType::Plane, sCollideConvexVsPlane);
		CollisionDispatch::sRegisterCastShape(s, EShapeSubType::Plane, sCastConvexVsPlane);

		CollisionDispatch::sRegisterCastShape(EShapeSubType::Plane, s, CollisionDispatch::sReversedCastShape);
		CollisionDispatch::sRegisterCollideShape(EShapeSubType::Plane, s, CollisionDispatch::sReversedCollideShape);
	}
}

JPH_NAMESPACE_END
