// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2024 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Collision/Shape/Shape.h>
#include <Jolt/Physics/Collision/PhysicsMaterial.h>

JPH_NAMESPACE_BEGIN

/// Class that constructs an EmptyShape
class JPH_EXPORT EmptyShapeSettings final : public ShapeSettings
{
	JPH_DECLARE_SERIALIZABLE_VIRTUAL(JPH_EXPORT, EmptyShapeSettings)

public:
							EmptyShapeSettings() = default;
	explicit				EmptyShapeSettings(Vec4Arg inCenterOfMass) : mCenterOfMass(inCenterOfMass) { }

	ShapeResult				Create() const override;

	Vec4					mCenterOfMass = Vec4::sZero();									///< Determines the center of mass for this shape
};

/// An empty shape that has no volume and collides with nothing.
///
/// Possible use cases:
/// - As a placeholder for a shape that will be created later. E.g. if you first need to create a body and only then know what shape it will have.
/// - If you need a kinematic body to attach a constraint to, but you don't want the body to collide with anything.
///
/// Note that, if possible, you should also put your body in an ObjectLayer that doesn't collide with anything.
/// This ensures that collisions will be filtered out at broad phase level instead of at narrow phase level, this is more efficient.
class JPH_EXPORT EmptyShape final : public Shape
{
public:
	// Constructor
							EmptyShape() : Shape(EShapeType::Empty, EShapeSubType::Empty) { }
	explicit				EmptyShape(Vec4Arg inCenterOfMass) : Shape(EShapeType::Empty, EShapeSubType::Empty), mCenterOfMass(inCenterOfMass) { }
							EmptyShape(const EmptyShapeSettings &inSettings, ShapeResult &outResult) : Shape(EShapeType::Empty, EShapeSubType::Empty, inSettings, outResult), mCenterOfMass(inSettings.mCenterOfMass) { outResult.Set(this); }

	// See: Shape
	Vec4					GetCenterOfMass() const override								{ return mCenterOfMass; }
	AABox					GetLocalBounds() const override									{ return { Vec4::sZero(), Vec4::sZero() }; }
	uint					GetSubShapeIDBitsRecursive() const override						{ return 0; }
	float					GetInnerRadius() const override									{ return 0.0f; }
	MassProperties			GetMassProperties() const override;
	const PhysicsMaterial *	GetMaterial([[maybe_unused]] const SubShapeID &inSubShapeID) const override { return PhysicsMaterial::sDefault; }
	virtual Vec4			GetSurfaceNormal(const SubShapeID &inSubShapeID, Vec4Arg inLocalSurfacePosition) const override { return Vec4::sZero(); }
	virtual void			GetSubmergedVolume(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const Plane &inSurface, float &outTotalVolume, float &outSubmergedVolume, Vec4 &outCenterOfBuoyancy, RVec4Arg inBaseOffset) const override { outTotalVolume = 0.0f; outSubmergedVolume = 0.0f; outCenterOfBuoyancy = Vec4::sZero(); }
#ifdef JPH_DEBUG_RENDERER
	virtual void			Draw([[maybe_unused]] DebugRenderer *inRenderer, [[maybe_unused]] RMat44Arg inCenterOfMassTransform, [[maybe_unused]] Vec4Arg inScale, [[maybe_unused]] ColorArg inColor, [[maybe_unused]] bool inUseMaterialColors, [[maybe_unused]] bool inDrawWireframe) const override;
#endif // JPH_DEBUG_RENDERER
	virtual bool			CastRay([[maybe_unused]] const RayCast &inRay, [[maybe_unused]] const SubShapeIDCreator &inSubShapeIDCreator, [[maybe_unused]] RayCastResult &ioHit) const override { return false; }
	virtual void			CastRay([[maybe_unused]] const RayCast &inRay, [[maybe_unused]] const RayCastSettings &inRayCastSettings, [[maybe_unused]] const SubShapeIDCreator &inSubShapeIDCreator, [[maybe_unused]] CastRayCollector &ioCollector, [[maybe_unused]] const ShapeFilter &inShapeFilter = { }) const override { /* Do nothing */ }
	virtual void			CollidePoint([[maybe_unused]] Vec4Arg inPoint, [[maybe_unused]] const SubShapeIDCreator &inSubShapeIDCreator, [[maybe_unused]] CollidePointCollector &ioCollector, [[maybe_unused]] const ShapeFilter &inShapeFilter = { }) const override { /* Do nothing */ }
	virtual void			CollideSoftBodyVertices([[maybe_unused]] RMat44Arg inCenterOfMassTransform, [[maybe_unused]] Vec4Arg inScale, [[maybe_unused]] const CollideSoftBodyVertexIterator &inVertices, [[maybe_unused]] uint inNumVertices, [[maybe_unused]] int inCollidingShapeIndex) const override { /* Do nothing */ }
	virtual void			GetTrianglesStart([[maybe_unused]] GetTrianglesContext &ioContext, [[maybe_unused]] const AABox &inBox, [[maybe_unused]] Vec4Arg inPositionCOM, [[maybe_unused]] RotorArg inRotation, [[maybe_unused]] Vec4Arg inScale) const override { /* Do nothing */ }
	virtual int				GetTrianglesNext([[maybe_unused]] GetTrianglesContext &ioContext, [[maybe_unused]] int inMaxTrianglesRequested, [[maybe_unused]] Float4 *outTriangleVertices, [[maybe_unused]] const PhysicsMaterial **outMaterials = nullptr) const override { return 0; }
	Stats					GetStats() const override										{ return { sizeof(*this), 0 }; }
	float					GetVolume() const override										{ return 0.0f; }
	bool					IsValidScale([[maybe_unused]] Vec4Arg inScale) const override	{ return true; }

	// Register shape functions with the registry
	static void				sRegister();

private:
	Vec4					mCenterOfMass = Vec4::sZero();
};

JPH_NAMESPACE_END
