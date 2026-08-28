// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/RegisterTypes.h>
#include <Jolt/Core/Factory.h>
#include <Jolt/Core/RTTI.h>
#include <Jolt/Core/TickCounter.h>
#include <Jolt/Physics/Collision/CollisionDispatch.h>
// TODO(4D): TriangleShape/TaperedCapsuleShape/CylinderShape/TaperedCylinderShape/MeshShape/
// HeightFieldShape are cut from the v1 build (see CLAUDE.md "Scope cuts"); their headers are still 3D.
#include <Jolt/Physics/Collision/Shape/PlaneShape.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/CapsuleShape.h>
#include <Jolt/Physics/Collision/Shape/ScaledShape.h>
#include <Jolt/Physics/Collision/Shape/ConvexHullShape.h>
#include <Jolt/Physics/Collision/Shape/RotatedTranslatedShape.h>
#include <Jolt/Physics/Collision/Shape/OffsetCenterOfMassShape.h>
#include <Jolt/Physics/Collision/Shape/MutableCompoundShape.h>
#include <Jolt/Physics/Collision/Shape/StaticCompoundShape.h>
#include <Jolt/Physics/Collision/Shape/EmptyShape.h>
#include <Jolt/Physics/Collision/PhysicsMaterialSimple.h>
// TODO(4D): SoftBodyShape cut for v1 (SoftBody is deeply 3D-coupled).

//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, Skeleton)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, SkeletalAnimation)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, RagdollSettings)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, PointConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, SixDOFConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, SliderConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, SwingTwistConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, DistanceConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, HingeConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, FixedConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, ConeConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, PathConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, PathConstraintPath)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, PathConstraintPathHermite)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, VehicleConstraintSettings)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, WheeledVehicleControllerSettings)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, WheelSettingsWV)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, TrackedVehicleControllerSettings)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, WheelSettingsTV)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, MotorcycleControllerSettings)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, RackAndPinionConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, GearConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, PulleyConstraintSettings)	// TODO(4D): joint constraints deferred to v1.1
JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, MotorSettings)
JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, PhysicsScene)
JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, PhysicsMaterial)
JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, GroupFilter)
JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, GroupFilterTable)
JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, BodyCreationSettings)
#ifdef JPH_ENABLE_SOFT_BODY // TODO(4D): SoftBody cut for v1, its sources are out of Jolt.cmake
JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, SoftBodyCreationSettings)
#endif // JPH_ENABLE_SOFT_BODY
//JPH_DECLARE_RTTI_WITH_NAMESPACE_FOR_FACTORY(JPH_EXPORT, JPH, HairSettings)	// TODO(4D): v1 scope cut, sources out of Jolt.cmake

JPH_NAMESPACE_BEGIN

bool VerifyJoltVersionIDInternal(uint64 inVersionID)
{
	return inVersionID == JPH_VERSION_ID;
}

void RegisterTypesInternal(uint64 inVersionID)
{
	// Version check
	if (!VerifyJoltVersionIDInternal(inVersionID))
	{
		Trace("Version mismatch, make sure you compile the client code with the same Jolt version and compiler definitions!");
		uint64 mismatch = JPH_VERSION_ID ^ inVersionID;
		if (mismatch & 0xffffff)
			Trace("Client reported version %d.%d.%d, library version is %d.%d.%d.",
				(inVersionID >> 16) & 0xff, (inVersionID >> 8) & 0xff, inVersionID & 0xff,
				JPH_VERSION_MAJOR, JPH_VERSION_MINOR, JPH_VERSION_PATCH);
		auto check_bit = [mismatch](int inBit, const char *inLabel) { if (mismatch & (uint64(1) << (inBit + 23))) Trace("Mismatching define %s.", inLabel); };
		check_bit(1, "JPH_DOUBLE_PRECISION");
		check_bit(2, "JPH_CROSS_PLATFORM_DETERMINISTIC");
		check_bit(3, "JPH_FLOATING_POINT_EXCEPTIONS_ENABLED");
		check_bit(4, "JPH_PROFILE_ENABLED");
		check_bit(5, "JPH_EXTERNAL_PROFILE");
		check_bit(6, "JPH_DEBUG_RENDERER");
		check_bit(7, "JPH_DISABLE_TEMP_ALLOCATOR");
		check_bit(8, "JPH_DISABLE_CUSTOM_ALLOCATOR");
		check_bit(9, "JPH_OBJECT_LAYER_BITS");
		check_bit(10, "JPH_ENABLE_ASSERTS");
		check_bit(11, "JPH_OBJECT_STREAM");
		std::abort();
	}

#ifndef JPH_DISABLE_CUSTOM_ALLOCATOR
	JPH_ASSERT(Allocate != nullptr && Reallocate != nullptr && Free != nullptr && AlignedAllocate != nullptr && AlignedFree != nullptr, "Need to supply an allocator first or call RegisterDefaultAllocator()");
#endif // !JPH_DISABLE_CUSTOM_ALLOCATOR

	JPH_ASSERT(Factory::sInstance != nullptr, "Need to create a factory first!");

	// Initialize dispatcher
	CollisionDispatch::sInit();

	// Register base classes first so that we can specialize them later
	CompoundShape::sRegister();
	ConvexShape::sRegister();

	// Register compounds before others so that we can specialize them later (register them in reverse order of collision complexity)
	MutableCompoundShape::sRegister();
	StaticCompoundShape::sRegister();

	// Leaf classes
	// TODO(4D): Triangle/TaperedCapsule/Cylinder/TaperedCylinder/Mesh/HeightField/SoftBody cut for v1.
	PlaneShape::sRegister();
	SphereShape::sRegister();
	BoxShape::sRegister();
	CapsuleShape::sRegister();
	ConvexHullShape::sRegister();

	// Register these last because their collision functions are simple so we want to execute them first (register them in reverse order of collision complexity)
	RotatedTranslatedShape::sRegister();
	OffsetCenterOfMassShape::sRegister();
	ScaledShape::sRegister();
	EmptyShape::sRegister();

	// Create list of all types
	const RTTI *types[] = {
//		JPH_RTTI(SkeletalAnimation),	// TODO(4D): v1 scope cut
//		JPH_RTTI(Skeleton),	// TODO(4D): v1 scope cut
		JPH_RTTI(CompoundShapeSettings),
		JPH_RTTI(StaticCompoundShapeSettings),
		JPH_RTTI(MutableCompoundShapeSettings),
		JPH_RTTI(PlaneShapeSettings),
		JPH_RTTI(SphereShapeSettings),
		JPH_RTTI(BoxShapeSettings),
		JPH_RTTI(CapsuleShapeSettings),
		JPH_RTTI(ScaledShapeSettings),
		JPH_RTTI(ConvexHullShapeSettings),
		JPH_RTTI(RotatedTranslatedShapeSettings),
		JPH_RTTI(OffsetCenterOfMassShapeSettings),
		JPH_RTTI(EmptyShapeSettings),
//		JPH_RTTI(RagdollSettings),	// TODO(4D): v1 scope cut
//		JPH_RTTI(PointConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(SixDOFConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(SliderConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(SwingTwistConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(DistanceConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(HingeConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(FixedConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(ConeConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(PathConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(VehicleConstraintSettings),	// TODO(4D): v1 scope cut
//		JPH_RTTI(WheeledVehicleControllerSettings),	// TODO(4D): v1 scope cut
//		JPH_RTTI(WheelSettingsWV),	// TODO(4D): v1 scope cut
//		JPH_RTTI(TrackedVehicleControllerSettings),	// TODO(4D): v1 scope cut
//		JPH_RTTI(WheelSettingsTV),	// TODO(4D): v1 scope cut
//		JPH_RTTI(MotorcycleControllerSettings),	// TODO(4D): v1 scope cut
//		JPH_RTTI(PathConstraintPath),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(PathConstraintPathHermite),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(RackAndPinionConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(GearConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
//		JPH_RTTI(PulleyConstraintSettings),	// TODO(4D): joint constraints deferred to v1.1
		JPH_RTTI(MotorSettings),
		JPH_RTTI(PhysicsScene),
		JPH_RTTI(PhysicsMaterial),
		JPH_RTTI(PhysicsMaterialSimple),
		JPH_RTTI(GroupFilter),
		JPH_RTTI(GroupFilterTable),
		JPH_RTTI(BodyCreationSettings),
#ifdef JPH_ENABLE_SOFT_BODY // TODO(4D): SoftBody cut for v1
		JPH_RTTI(SoftBodyCreationSettings),
#endif // JPH_ENABLE_SOFT_BODY
	};

	// Register them all
	Factory::sInstance->Register(types, (uint)std::size(types));

	// Initialize default physics material
	if (PhysicsMaterial::sDefault == nullptr)
		PhysicsMaterial::sDefault = new PhysicsMaterialSimple("Default", Color::sGrey);
}

// TODO(4D): Hair is a v1 scope cut (CLAUDE.md "Scope cuts for v1"); its sources are out of
// Jolt.cmake, so RegisterHair is a no-op until Hair is ported.
void RegisterHair()
{
}

void UnregisterTypes()
{
	// Unregister all types
	if (Factory::sInstance != nullptr)
		Factory::sInstance->Clear();

	// Delete default physics material
	PhysicsMaterial::sDefault = nullptr;
}

JPH_NAMESPACE_END
