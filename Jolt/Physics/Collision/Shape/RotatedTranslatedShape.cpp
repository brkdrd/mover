// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Collision/Shape/RotatedTranslatedShape.h>
#include <Jolt/Physics/Collision/CollisionDispatch.h>
#include <Jolt/Physics/Collision/RayCast.h>
#include <Jolt/Physics/Collision/ShapeCast.h>
#include <Jolt/Physics/Collision/TransformedShape.h>
#include <Jolt/Core/StreamIn.h>
#include <Jolt/Core/StreamOut.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>

JPH_NAMESPACE_BEGIN

JPH_IMPLEMENT_SERIALIZABLE_VIRTUAL(RotatedTranslatedShapeSettings)
{
	JPH_ADD_BASE_CLASS(RotatedTranslatedShapeSettings, DecoratedShapeSettings)

	JPH_ADD_ATTRIBUTE(RotatedTranslatedShapeSettings, mPosition)
	JPH_ADD_ATTRIBUTE(RotatedTranslatedShapeSettings, mRotation)
}

ShapeSettings::ShapeResult RotatedTranslatedShapeSettings::Create() const
{
	if (mCachedResult.IsEmpty())
		Ref<Shape> shape = new RotatedTranslatedShape(*this, mCachedResult);
	return mCachedResult;
}

RotatedTranslatedShape::RotatedTranslatedShape(const RotatedTranslatedShapeSettings &inSettings, ShapeResult &outResult) :
	DecoratedShape(EShapeSubType::RotatedTranslated, inSettings, outResult)
{
	if (outResult.HasError())
		return;

	// Calculate center of mass position
	mCenterOfMass = inSettings.mPosition + inSettings.mRotation * mInnerShape->GetCenterOfMass();

	// Store rotation (position is always zero because we center around the center of mass)
	mRotation = inSettings.mRotation;
	mIsRotationIdentity = mRotation.IsClose(Rotor::sIdentity());

	outResult.Set(this);
}

RotatedTranslatedShape::RotatedTranslatedShape(Vec4Arg inPosition, RotorArg inRotation, const Shape *inShape) :
	DecoratedShape(EShapeSubType::RotatedTranslated, inShape)
{
	// Calculate center of mass position
	mCenterOfMass = inPosition + inRotation * mInnerShape->GetCenterOfMass();

	// Store rotation (position is always zero because we center around the center of mass)
	mRotation = inRotation;
	mIsRotationIdentity = mRotation.IsClose(Rotor::sIdentity());
}

MassProperties RotatedTranslatedShape::GetMassProperties() const
{
	// Rotate inertia of child into place
	MassProperties p = mInnerShape->GetMassProperties();
	p.Rotate(Mat44::sRotation(mRotation));
	return p;
}

AABox RotatedTranslatedShape::GetLocalBounds() const
{
	return mInnerShape->GetLocalBounds().Transformed(Mat44::sRotation(mRotation));
}

AABox RotatedTranslatedShape::GetWorldSpaceBounds(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale) const
{
	RMat44 transform = inCenterOfMassTransform * Mat44::sRotation(mRotation);
	return mInnerShape->GetWorldSpaceBounds(transform, TransformScale(inScale));
}

TransformedShape RotatedTranslatedShape::GetSubShapeTransformedShape(const SubShapeID &inSubShapeID, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale, SubShapeID &outRemainder) const
{
	// We don't use any bits in the sub shape ID
	outRemainder = inSubShapeID;

	TransformedShape ts(RVec4(inPositionCOM), inRotation * mRotation, mInnerShape, BodyID());
	ts.SetShapeScale(TransformScale(inScale));
	return ts;
}

Vec4 RotatedTranslatedShape::GetSurfaceNormal(const SubShapeID &inSubShapeID, Vec4Arg inLocalSurfacePosition) const
{
	// Transform surface position to local space of the child (R^T) and pass call on
	Mat44 rot = Mat44::sRotation(mRotation);
	Vec4 normal = mInnerShape->GetSurfaceNormal(inSubShapeID, rot.Transposed() * inLocalSurfacePosition);

	// Transform normal back to this shape's space (R)
	return rot * normal;
}

void RotatedTranslatedShape::GetSupportingFace(const SubShapeID &inSubShapeID, Vec4Arg inDirection, Vec4Arg inScale, RMat44Arg inCenterOfMassTransform, SupportingFace &outVertices) const
{
	Mat44 rot = Mat44::sRotation(mRotation);
	mInnerShape->GetSupportingFace(inSubShapeID, rot.Transposed() * inDirection, TransformScale(inScale), inCenterOfMassTransform * rot, outVertices);
}

void RotatedTranslatedShape::GetSubmergedVolume(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const Plane &inSurface, float &outTotalVolume, float &outSubmergedVolume, Vec4 &outCenterOfBuoyancy, RVec4Arg inBaseOffset) const
{
	// Get center of mass transform of child
	RMat44 transform = inCenterOfMassTransform * Mat44::sRotation(mRotation);

	// Recurse to child
	mInnerShape->GetSubmergedVolume(transform, TransformScale(inScale), inSurface, outTotalVolume, outSubmergedVolume, outCenterOfBuoyancy, inBaseOffset);
}

#ifdef JPH_DEBUG_RENDERER
void RotatedTranslatedShape::Draw(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inUseMaterialColors, bool inDrawWireframe) const
{
	mInnerShape->Draw(inRenderer, inCenterOfMassTransform * Mat44::sRotation(mRotation), TransformScale(inScale), inColor, inUseMaterialColors, inDrawWireframe);
}

void RotatedTranslatedShape::DrawGetSupportFunction(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inDrawSupportDirection) const
{
	mInnerShape->DrawGetSupportFunction(inRenderer, inCenterOfMassTransform * Mat44::sRotation(mRotation), TransformScale(inScale), inColor, inDrawSupportDirection);
}

void RotatedTranslatedShape::DrawGetSupportingFace(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale) const
{
	mInnerShape->DrawGetSupportingFace(inRenderer, inCenterOfMassTransform * Mat44::sRotation(mRotation), TransformScale(inScale));
}
#endif // JPH_DEBUG_RENDERER

bool RotatedTranslatedShape::CastRay(const RayCast &inRay, const SubShapeIDCreator &inSubShapeIDCreator, RayCastResult &ioHit) const
{
	// Transform the ray into the space of the child (inverse rotation)
	RayCast ray = inRay.Transformed(RMat44::sRotation(mRotation.Reversed()));

	return mInnerShape->CastRay(ray, inSubShapeIDCreator, ioHit);
}

void RotatedTranslatedShape::CastRay(const RayCast &inRay, const RayCastSettings &inRayCastSettings, const SubShapeIDCreator &inSubShapeIDCreator, CastRayCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	// Transform the ray into the space of the child (inverse rotation)
	RayCast ray = inRay.Transformed(RMat44::sRotation(mRotation.Reversed()));

	return mInnerShape->CastRay(ray, inRayCastSettings, inSubShapeIDCreator, ioCollector, inShapeFilter);
}

void RotatedTranslatedShape::CollidePoint(Vec4Arg inPoint, const SubShapeIDCreator &inSubShapeIDCreator, CollidePointCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	// Transform the point into the space of the child (inverse rotation)
	Mat44 transform = Mat44::sRotation(mRotation).Transposed();
	mInnerShape->CollidePoint(transform * inPoint, inSubShapeIDCreator, ioCollector, inShapeFilter);
}

void RotatedTranslatedShape::CollideSoftBodyVertices(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const CollideSoftBodyVertexIterator &inVertices, uint inNumVertices, int inCollidingShapeIndex) const
{
	mInnerShape->CollideSoftBodyVertices(inCenterOfMassTransform * Mat44::sRotation(mRotation), inScale, inVertices, inNumVertices, inCollidingShapeIndex);
}

void RotatedTranslatedShape::CollectTransformedShapes(const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale, const SubShapeIDCreator &inSubShapeIDCreator, TransformedShapeCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	mInnerShape->CollectTransformedShapes(inBox, inPositionCOM, inRotation * mRotation, TransformScale(inScale), inSubShapeIDCreator, ioCollector, inShapeFilter);
}

void RotatedTranslatedShape::TransformShape(RMat44Arg inCenterOfMassTransform, TransformedShapeCollector &ioCollector) const
{
	mInnerShape->TransformShape(inCenterOfMassTransform * Mat44::sRotation(mRotation), ioCollector);
}

void RotatedTranslatedShape::sCollideRotatedTranslatedVsShape(const Shape *inShape1, const Shape *inShape2, Vec4Arg inScale1, Vec4Arg inScale2, RMat44Arg inCenterOfMassTransform1, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, const ShapeFilter &inShapeFilter)
{
	JPH_ASSERT(inShape1->GetSubType() == EShapeSubType::RotatedTranslated);
	const RotatedTranslatedShape *shape1 = static_cast<const RotatedTranslatedShape *>(inShape1);

	// Get world transform of 1
	RMat44 transform1 = inCenterOfMassTransform1 * Mat44::sRotation(shape1->mRotation);

	CollisionDispatch::sCollideShapeVsShape(shape1->mInnerShape, inShape2, shape1->TransformScale(inScale1), inScale2, transform1, inCenterOfMassTransform2, inSubShapeIDCreator1, inSubShapeIDCreator2, inCollideShapeSettings, ioCollector, inShapeFilter);
}

void RotatedTranslatedShape::sCollideShapeVsRotatedTranslated(const Shape *inShape1, const Shape *inShape2, Vec4Arg inScale1, Vec4Arg inScale2, RMat44Arg inCenterOfMassTransform1, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, const ShapeFilter &inShapeFilter)
{
	JPH_ASSERT(inShape2->GetSubType() == EShapeSubType::RotatedTranslated);
	const RotatedTranslatedShape *shape2 = static_cast<const RotatedTranslatedShape *>(inShape2);

	// Get world transform of 2
	RMat44 transform2 = inCenterOfMassTransform2 * Mat44::sRotation(shape2->mRotation);

	CollisionDispatch::sCollideShapeVsShape(inShape1, shape2->mInnerShape, inScale1, shape2->TransformScale(inScale2), inCenterOfMassTransform1, transform2, inSubShapeIDCreator1, inSubShapeIDCreator2, inCollideShapeSettings, ioCollector, inShapeFilter);
}

void RotatedTranslatedShape::sCollideRotatedTranslatedVsRotatedTranslated(const Shape *inShape1, const Shape *inShape2, Vec4Arg inScale1, Vec4Arg inScale2, RMat44Arg inCenterOfMassTransform1, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, const ShapeFilter &inShapeFilter)
{
	JPH_ASSERT(inShape1->GetSubType() == EShapeSubType::RotatedTranslated);
	const RotatedTranslatedShape *shape1 = static_cast<const RotatedTranslatedShape *>(inShape1);
	JPH_ASSERT(inShape2->GetSubType() == EShapeSubType::RotatedTranslated);
	const RotatedTranslatedShape *shape2 = static_cast<const RotatedTranslatedShape *>(inShape2);

	// Get world transform of 1 and 2
	RMat44 transform1 = inCenterOfMassTransform1 * Mat44::sRotation(shape1->mRotation);
	RMat44 transform2 = inCenterOfMassTransform2 * Mat44::sRotation(shape2->mRotation);

	CollisionDispatch::sCollideShapeVsShape(shape1->mInnerShape, shape2->mInnerShape, shape1->TransformScale(inScale1), shape2->TransformScale(inScale2), transform1, transform2, inSubShapeIDCreator1, inSubShapeIDCreator2, inCollideShapeSettings, ioCollector, inShapeFilter);
}

void RotatedTranslatedShape::sCastRotatedTranslatedVsShape(const ShapeCast &inShapeCast, const ShapeCastSettings &inShapeCastSettings, const Shape *inShape, Vec4Arg inScale, const ShapeFilter &inShapeFilter, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, CastShapeCollector &ioCollector)
{
	// Fetch rotated translated shape from cast shape
	JPH_ASSERT(inShapeCast.mShape->GetSubType() == EShapeSubType::RotatedTranslated);
	const RotatedTranslatedShape *shape1 = static_cast<const RotatedTranslatedShape *>(inShapeCast.mShape);

	// Transform the shape cast and update the shape
	RMat44 transform = inShapeCast.mCenterOfMassStart * Mat44::sRotation(shape1->mRotation);
	Vec4 scale = shape1->TransformScale(inShapeCast.mScale);
	ShapeCast shape_cast(shape1->mInnerShape, scale, transform, inShapeCast.mDirection);

	CollisionDispatch::sCastShapeVsShapeLocalSpace(shape_cast, inShapeCastSettings, inShape, inScale, inShapeFilter, inCenterOfMassTransform2, inSubShapeIDCreator1, inSubShapeIDCreator2, ioCollector);
}

void RotatedTranslatedShape::sCastShapeVsRotatedTranslated(const ShapeCast &inShapeCast, const ShapeCastSettings &inShapeCastSettings, const Shape *inShape, Vec4Arg inScale, const ShapeFilter &inShapeFilter, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, CastShapeCollector &ioCollector)
{
	JPH_ASSERT(inShape->GetSubType() == EShapeSubType::RotatedTranslated);
	const RotatedTranslatedShape *shape = static_cast<const RotatedTranslatedShape *>(inShape);

	// Determine the local transform (pure rotation, so its inverse is the transpose)
	Mat44 local_transform = Mat44::sRotation(shape->mRotation);

	// Transform the shape cast into the space of the child (inverse rotation)
	ShapeCast shape_cast = inShapeCast.PostTransformed(RMat44::sRotation(shape->mRotation.Reversed()));

	CollisionDispatch::sCastShapeVsShapeLocalSpace(shape_cast, inShapeCastSettings, shape->mInnerShape, shape->TransformScale(inScale), inShapeFilter, inCenterOfMassTransform2 * local_transform, inSubShapeIDCreator1, inSubShapeIDCreator2, ioCollector);
}

void RotatedTranslatedShape::sCastRotatedTranslatedVsRotatedTranslated(const ShapeCast &inShapeCast, const ShapeCastSettings &inShapeCastSettings, const Shape *inShape, Vec4Arg inScale, const ShapeFilter &inShapeFilter, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, CastShapeCollector &ioCollector)
{
	JPH_ASSERT(inShapeCast.mShape->GetSubType() == EShapeSubType::RotatedTranslated);
	const RotatedTranslatedShape *shape1 = static_cast<const RotatedTranslatedShape *>(inShapeCast.mShape);
	JPH_ASSERT(inShape->GetSubType() == EShapeSubType::RotatedTranslated);
	const RotatedTranslatedShape *shape2 = static_cast<const RotatedTranslatedShape *>(inShape);

	// Determine the local transform of shape 2 and its inverse (pure rotation)
	Mat44 local_transform2 = Mat44::sRotation(shape2->mRotation);
	Mat44 local_transform2_transposed = local_transform2.Transposed();
	RMat44 inv_local_transform2 = RMat44::sRotation(shape2->mRotation.Reversed());

	// Transform the shape cast and update the shape (express relative to shape 2's child frame)
	RMat44 transform = (inv_local_transform2 * inShapeCast.mCenterOfMassStart) * Mat44::sRotation(shape1->mRotation);
	Vec4 scale = shape1->TransformScale(inShapeCast.mScale);
	ShapeCast shape_cast(shape1->mInnerShape, scale, transform, local_transform2_transposed * inShapeCast.mDirection);

	CollisionDispatch::sCastShapeVsShapeLocalSpace(shape_cast, inShapeCastSettings, shape2->mInnerShape, shape2->TransformScale(inScale), inShapeFilter, inCenterOfMassTransform2 * local_transform2, inSubShapeIDCreator1, inSubShapeIDCreator2, ioCollector);
}

void RotatedTranslatedShape::SaveBinaryState(StreamOut &inStream) const
{
	DecoratedShape::SaveBinaryState(inStream);

	inStream.Write(mCenterOfMass);
	inStream.Write(mRotation);
}

void RotatedTranslatedShape::RestoreBinaryState(StreamIn &inStream)
{
	DecoratedShape::RestoreBinaryState(inStream);

	inStream.Read(mCenterOfMass);
	inStream.Read(mRotation);
	mIsRotationIdentity = mRotation.IsClose(Rotor::sIdentity());
}

bool RotatedTranslatedShape::IsValidScale(Vec4Arg inScale) const
{
	if (!Shape::IsValidScale(inScale))
		return false;

	if (mIsRotationIdentity || ScaleHelpers::IsUniformScale(inScale))
		return mInnerShape->IsValidScale(inScale);

	if (!ScaleHelpers::CanScaleBeRotated(mRotation, inScale))
		return false;

	return mInnerShape->IsValidScale(ScaleHelpers::RotateScale(mRotation, inScale));
}

Vec4 RotatedTranslatedShape::MakeScaleValid(Vec4Arg inScale) const
{
	Vec4 scale = ScaleHelpers::MakeNonZeroScale(inScale);

	if (mIsRotationIdentity || ScaleHelpers::IsUniformScale(scale))
		return mInnerShape->MakeScaleValid(scale);

	if (ScaleHelpers::CanScaleBeRotated(mRotation, scale))
		return ScaleHelpers::RotateScale(mRotation.Reversed(), mInnerShape->MakeScaleValid(ScaleHelpers::RotateScale(mRotation, scale)));

	Vec4 abs_uniform_scale = ScaleHelpers::MakeUniformScale(scale.Abs());
	Vec4 uniform_scale = scale.GetSign() * abs_uniform_scale;
	if (ScaleHelpers::CanScaleBeRotated(mRotation, uniform_scale))
		return uniform_scale;

	return Sign(scale.GetX()) * abs_uniform_scale;
}

void RotatedTranslatedShape::sRegister()
{
	ShapeFunctions &f = ShapeFunctions::sGet(EShapeSubType::RotatedTranslated);
	f.mConstruct = []() -> Shape * { return new RotatedTranslatedShape; };
	f.mColor = Color::sBlue;

	for (EShapeSubType s : sAllSubShapeTypes)
	{
		CollisionDispatch::sRegisterCollideShape(EShapeSubType::RotatedTranslated, s, sCollideRotatedTranslatedVsShape);
		CollisionDispatch::sRegisterCollideShape(s, EShapeSubType::RotatedTranslated, sCollideShapeVsRotatedTranslated);
		CollisionDispatch::sRegisterCastShape(EShapeSubType::RotatedTranslated, s, sCastRotatedTranslatedVsShape);
		CollisionDispatch::sRegisterCastShape(s, EShapeSubType::RotatedTranslated, sCastShapeVsRotatedTranslated);
	}

	CollisionDispatch::sRegisterCollideShape(EShapeSubType::RotatedTranslated, EShapeSubType::RotatedTranslated, sCollideRotatedTranslatedVsRotatedTranslated);
	CollisionDispatch::sRegisterCastShape(EShapeSubType::RotatedTranslated, EShapeSubType::RotatedTranslated, sCastRotatedTranslatedVsRotatedTranslated);
}

JPH_NAMESPACE_END
