// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Collision/Shape/ConvexShape.h>

JPH_NAMESPACE_BEGIN

/// Class that constructs a CapsuleShape
class JPH_EXPORT CapsuleShapeSettings final : public ConvexShapeSettings
{
	JPH_DECLARE_SERIALIZABLE_VIRTUAL(JPH_EXPORT, CapsuleShapeSettings)

public:
	/// Default constructor for deserialization
							CapsuleShapeSettings() = default;

	/// Create a capsule centered around the origin with one sphere cap at (0, -inHalfHeightOfCylinder, 0) and the other at (0, inHalfHeightOfCylinder, 0)
							CapsuleShapeSettings(float inHalfHeightOfCylinder, float inRadius, const PhysicsMaterial *inMaterial = nullptr) : ConvexShapeSettings(inMaterial), mRadius(inRadius), mHalfHeightOfCylinder(inHalfHeightOfCylinder) { }

	/// Check if this is a valid capsule shape
	bool					IsValid() const															{ return mRadius > 0.0f && mHalfHeightOfCylinder >= 0.0f; }

	/// Checks if the settings of this capsule make this shape a sphere
	bool					IsSphere() const														{ return mHalfHeightOfCylinder == 0.0f; }

	// See: ShapeSettings
	virtual ShapeResult		Create() const override;

	float					mRadius = 0.0f;
	float					mHalfHeightOfCylinder = 0.0f;
};

/// A capsule, implemented as a line segment with convex radius
class JPH_EXPORT CapsuleShape final : public ConvexShape
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructor
							CapsuleShape() : ConvexShape(EShapeSubType::Capsule) { }
							CapsuleShape(const CapsuleShapeSettings &inSettings, ShapeResult &outResult);

	/// Create a capsule centered around the origin with one sphere cap at (0, -inHalfHeightOfCylinder, 0) and the other at (0, inHalfHeightOfCylinder, 0)
							CapsuleShape(float inHalfHeightOfCylinder, float inRadius, const PhysicsMaterial *inMaterial = nullptr) : ConvexShape(EShapeSubType::Capsule, inMaterial), mRadius(inRadius), mHalfHeightOfCylinder(inHalfHeightOfCylinder) { JPH_ASSERT(inHalfHeightOfCylinder > 0.0f); JPH_ASSERT(inRadius > 0.0f); }

	/// Radius of the cylinder
	float					GetRadius() const														{ return mRadius; }

	/// Get half of the height of the cylinder
	float					GetHalfHeightOfCylinder() const											{ return mHalfHeightOfCylinder; }

	// See Shape::GetLocalBounds
	virtual AABox			GetLocalBounds() const override;

	// See Shape::GetWorldSpaceBounds
	virtual AABox			GetWorldSpaceBounds(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale) const override;
	using Shape::GetWorldSpaceBounds;

	// See Shape::GetInnerRadius
	virtual float			GetInnerRadius() const override											{ return mRadius; }

	// See Shape::GetMassProperties
	virtual MassProperties	GetMassProperties() const override;

	// See Shape::GetSurfaceNormal
	virtual Vec4			GetSurfaceNormal(const SubShapeID &inSubShapeID, Vec4Arg inLocalSurfacePosition) const override;

	// See Shape::GetSupportingFace
	virtual void			GetSupportingFace(const SubShapeID &inSubShapeID, Vec4Arg inDirection, Vec4Arg inScale, RMat44Arg inCenterOfMassTransform, SupportingFace &outVertices) const override;

	// See ConvexShape::GetSupportFunction
	virtual const Support *	GetSupportFunction(ESupportMode inMode, SupportBuffer &inBuffer, Vec4Arg inScale) const override;

#ifdef JPH_DEBUG_RENDERER
	// See Shape::Draw
	virtual void			Draw(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inUseMaterialColors, bool inDrawWireframe) const override;
#endif // JPH_DEBUG_RENDERER

	// See Shape::CastRay
	using ConvexShape::CastRay;
	virtual bool			CastRay(const RayCast &inRay, const SubShapeIDCreator &inSubShapeIDCreator, RayCastResult &ioHit) const override;

	// See: Shape::CollidePoint
	virtual void			CollidePoint(Vec4Arg inPoint, const SubShapeIDCreator &inSubShapeIDCreator, CollidePointCollector &ioCollector, const ShapeFilter &inShapeFilter = { }) const override;

	// See: Shape::CollideSoftBodyVertices
	virtual void			CollideSoftBodyVertices(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const CollideSoftBodyVertexIterator &inVertices, uint inNumVertices, int inCollidingShapeIndex) const override;

	// See Shape::GetTrianglesStart
	virtual void			GetTrianglesStart(GetTrianglesContext &ioContext, const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale) const override;

	// See Shape::GetTrianglesNext
	virtual int				GetTrianglesNext(GetTrianglesContext &ioContext, int inMaxTrianglesRequested, Float4 *outTriangleVertices, const PhysicsMaterial **outMaterials = nullptr) const override;

	// See Shape
	virtual void			SaveBinaryState(StreamOut &inStream) const override;

	// See Shape::GetStats
	virtual Stats			GetStats() const override												{ return Stats(sizeof(*this), 0); }

	// See Shape::GetVolume
	// 4D capsule = 4D cylinder (3-ball cross section swept along Y) + one 4-ball (two half-4-ball caps).
	virtual float			GetVolume() const override												{ return (8.0f / 3.0f) * JPH_PI * Cubed(mRadius) * mHalfHeightOfCylinder + 0.5f * JPH_PI * JPH_PI * Square(Square(mRadius)); }

	// See Shape::IsValidScale
	virtual bool			IsValidScale(Vec4Arg inScale) const override;

	// See Shape::MakeScaleValid
	virtual Vec4			MakeScaleValid(Vec4Arg inScale) const override;

	// Register shape functions with the registry
	static void				sRegister();

protected:
	// See: Shape::RestoreBinaryState
	virtual void			RestoreBinaryState(StreamIn &inStream) override;

private:
	// Classes for GetSupportFunction
	class					CapsuleNoConvex;
	class					CapsuleWithConvex;

	float					mRadius = 0.0f;
	float					mHalfHeightOfCylinder = 0.0f;
};

JPH_NAMESPACE_END
