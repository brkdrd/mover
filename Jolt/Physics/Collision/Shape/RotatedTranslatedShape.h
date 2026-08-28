// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Collision/Shape/DecoratedShape.h>
#include <Jolt/Physics/Collision/Shape/ScaleHelpers.h>

JPH_NAMESPACE_BEGIN

class CollideShapeSettings;

/// Class that constructs a RotatedTranslatedShape
class JPH_EXPORT RotatedTranslatedShapeSettings final : public DecoratedShapeSettings
{
	JPH_DECLARE_SERIALIZABLE_VIRTUAL(JPH_EXPORT, RotatedTranslatedShapeSettings)

public:
	/// Constructor
									RotatedTranslatedShapeSettings() = default;

	/// Construct with shape settings, can be serialized.
									RotatedTranslatedShapeSettings(Vec4Arg inPosition, RotorArg inRotation, const ShapeSettings *inShape) : DecoratedShapeSettings(inShape), mPosition(inPosition), mRotation(inRotation) { }

	/// Variant that uses a concrete shape, which means this object cannot be serialized.
									RotatedTranslatedShapeSettings(Vec4Arg inPosition, RotorArg inRotation, const Shape *inShape): DecoratedShapeSettings(inShape), mPosition(inPosition), mRotation(inRotation) { }

	// See: ShapeSettings
	virtual ShapeResult				Create() const override;

	Vec4							mPosition;												///< Position of the sub shape (4D)
	Rotor							mRotation;												///< Rotation of the sub shape (Cl(4,0) rotor)
};

/// A rotated translated shape will rotate and translate a child shape.
/// Shifts the child object so that it is centered around the center of mass.
class JPH_EXPORT RotatedTranslatedShape final : public DecoratedShape
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructor
									RotatedTranslatedShape() : DecoratedShape(EShapeSubType::RotatedTranslated) { }
									RotatedTranslatedShape(const RotatedTranslatedShapeSettings &inSettings, ShapeResult &outResult);
									RotatedTranslatedShape(Vec4Arg inPosition, RotorArg inRotation, const Shape *inShape);

	/// Access the rotation that is applied to the inner shape
	Rotor							GetRotation() const										{ return mRotation; }

	/// Access the translation that has been applied to the inner shape
	Vec4							GetPosition() const										{ return mCenterOfMass - mRotation * mInnerShape->GetCenterOfMass(); }

	// See Shape::GetCenterOfMass
	virtual Vec4					GetCenterOfMass() const override						{ return mCenterOfMass; }

	// See Shape::GetLocalBounds
	virtual AABox					GetLocalBounds() const override;

	// See Shape::GetWorldSpaceBounds
	virtual AABox					GetWorldSpaceBounds(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale) const override;
	using Shape::GetWorldSpaceBounds;

	// See Shape::GetInnerRadius
	virtual float					GetInnerRadius() const override							{ return mInnerShape->GetInnerRadius(); }

	// See Shape::GetMassProperties
	virtual MassProperties			GetMassProperties() const override;

	// See Shape::GetSubShapeTransformedShape
	virtual TransformedShape		GetSubShapeTransformedShape(const SubShapeID &inSubShapeID, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale, SubShapeID &outRemainder) const override;

	// See Shape::GetSurfaceNormal
	virtual Vec4					GetSurfaceNormal(const SubShapeID &inSubShapeID, Vec4Arg inLocalSurfacePosition) const override;

	// See Shape::GetSupportingFace
	virtual void					GetSupportingFace(const SubShapeID &inSubShapeID, Vec4Arg inDirection, Vec4Arg inScale, RMat44Arg inCenterOfMassTransform, SupportingFace &outVertices) const override;

	// See Shape::GetSubmergedVolume
	virtual void					GetSubmergedVolume(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const Plane &inSurface, float &outTotalVolume, float &outSubmergedVolume, Vec4 &outCenterOfBuoyancy, RVec4Arg inBaseOffset) const override;

#ifdef JPH_DEBUG_RENDERER
	// See Shape::Draw
	virtual void					Draw(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inUseMaterialColors, bool inDrawWireframe) const override;

	// See Shape::DrawGetSupportFunction
	virtual void					DrawGetSupportFunction(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inDrawSupportDirection) const override;

	// See Shape::DrawGetSupportingFace
	virtual void					DrawGetSupportingFace(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale) const override;
#endif // JPH_DEBUG_RENDERER

	// See Shape::CastRay
	virtual bool					CastRay(const RayCast &inRay, const SubShapeIDCreator &inSubShapeIDCreator, RayCastResult &ioHit) const override;
	virtual void					CastRay(const RayCast &inRay, const RayCastSettings &inRayCastSettings, const SubShapeIDCreator &inSubShapeIDCreator, CastRayCollector &ioCollector, const ShapeFilter &inShapeFilter = { }) const override;

	// See: Shape::CollidePoint
	virtual void					CollidePoint(Vec4Arg inPoint, const SubShapeIDCreator &inSubShapeIDCreator, CollidePointCollector &ioCollector, const ShapeFilter &inShapeFilter = { }) const override;

	// See: Shape::CollideSoftBodyVertices
	virtual void					CollideSoftBodyVertices(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const CollideSoftBodyVertexIterator &inVertices, uint inNumVertices, int inCollidingShapeIndex) const override;

	// See Shape::CollectTransformedShapes
	virtual void					CollectTransformedShapes(const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale, const SubShapeIDCreator &inSubShapeIDCreator, TransformedShapeCollector &ioCollector, const ShapeFilter &inShapeFilter) const override;

	// See Shape::TransformShape
	virtual void					TransformShape(RMat44Arg inCenterOfMassTransform, TransformedShapeCollector &ioCollector) const override;

	// See Shape::GetTrianglesStart
	virtual void					GetTrianglesStart(GetTrianglesContext &ioContext, const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale) const override { JPH_ASSERT(false, "Cannot call on non-leaf shapes, use CollectTransformedShapes to collect the leaves first!"); }

	// See Shape::GetTrianglesNext
	virtual int						GetTrianglesNext(GetTrianglesContext &ioContext, int inMaxTrianglesRequested, Float4 *outTriangleVertices, const PhysicsMaterial **outMaterials = nullptr) const override { JPH_ASSERT(false, "Cannot call on non-leaf shapes, use CollectTransformedShapes to collect the leaves first!"); return 0; }

	// See Shape
	virtual void					SaveBinaryState(StreamOut &inStream) const override;

	// See Shape::GetStats
	virtual Stats					GetStats() const override								{ return Stats(sizeof(*this), 0); }

	// See Shape::GetVolume
	virtual float					GetVolume() const override								{ return mInnerShape->GetVolume(); }

	// See Shape::IsValidScale
	virtual bool					IsValidScale(Vec4Arg inScale) const override;

	// See Shape::MakeScaleValid
	virtual Vec4					MakeScaleValid(Vec4Arg inScale) const override;

	/// Transform the scale to the local space of the child shape
	inline Vec4						TransformScale(Vec4Arg inScale) const
	{
		// We don't need to transform uniform scale or if the rotation is identity
		if (mIsRotationIdentity || ScaleHelpers::IsUniformScale(inScale))
			return inScale;

		return ScaleHelpers::RotateScale(mRotation, inScale);
	}

	// Register shape functions with the registry
	static void						sRegister();

protected:
	// See: Shape::RestoreBinaryState
	virtual void					RestoreBinaryState(StreamIn &inStream) override;

private:
	// Helper functions called by CollisionDispatch
	static void						sCollideRotatedTranslatedVsShape(const Shape *inShape1, const Shape *inShape2, Vec4Arg inScale1, Vec4Arg inScale2, RMat44Arg inCenterOfMassTransform1, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, const ShapeFilter &inShapeFilter);
	static void						sCollideShapeVsRotatedTranslated(const Shape *inShape1, const Shape *inShape2, Vec4Arg inScale1, Vec4Arg inScale2, RMat44Arg inCenterOfMassTransform1, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, const ShapeFilter &inShapeFilter);
	static void						sCollideRotatedTranslatedVsRotatedTranslated(const Shape *inShape1, const Shape *inShape2, Vec4Arg inScale1, Vec4Arg inScale2, RMat44Arg inCenterOfMassTransform1, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, const CollideShapeSettings &inCollideShapeSettings, CollideShapeCollector &ioCollector, const ShapeFilter &inShapeFilter);
	static void						sCastRotatedTranslatedVsShape(const ShapeCast &inShapeCast, const ShapeCastSettings &inShapeCastSettings, const Shape *inShape, Vec4Arg inScale, const ShapeFilter &inShapeFilter, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, CastShapeCollector &ioCollector);
	static void						sCastShapeVsRotatedTranslated(const ShapeCast &inShapeCast, const ShapeCastSettings &inShapeCastSettings, const Shape *inShape, Vec4Arg inScale, const ShapeFilter &inShapeFilter, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, CastShapeCollector &ioCollector);
	static void						sCastRotatedTranslatedVsRotatedTranslated(const ShapeCast &inShapeCast, const ShapeCastSettings &inShapeCastSettings, const Shape *inShape, Vec4Arg inScale, const ShapeFilter &inShapeFilter, RMat44Arg inCenterOfMassTransform2, const SubShapeIDCreator &inSubShapeIDCreator1, const SubShapeIDCreator &inSubShapeIDCreator2, CastShapeCollector &ioCollector);

	bool							mIsRotationIdentity;									///< If mRotation is close to identity (put here because it falls in padding bytes)
	Vec4							mCenterOfMass;											///< Position of the center of mass
	Rotor							mRotation;												///< Rotation of the child shape
};

JPH_NAMESPACE_END
