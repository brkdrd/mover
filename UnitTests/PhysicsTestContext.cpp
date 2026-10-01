// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include <Jolt/Physics/Constraints/ContactConstraintManager.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>
#include <Jolt/Core/JobSystemThreadPool.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Core/StreamWrapper.h>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRendererRecorder.h>
#endif

PhysicsTestContext::PhysicsTestContext(float inDeltaTime, int inCollisionSteps, int inWorkerThreads, uint inMaxBodies, uint inMaxBodyPairs, uint inMaxContactConstraints) :
#ifdef JPH_DISABLE_TEMP_ALLOCATOR
	mTempAllocator(new TempAllocatorMalloc()),
#else
	// TODO(4D): was 4 MB. A ContactConstraint is 4128 bytes in 4D (measured) against a few hundred
	// in 3D: it holds 6 manifold points instead of 4, each with 4 AxisConstraintParts instead of 3
	// (friction has a 3D tangent space here), each holding bivector Jacobians of 8 floats instead of
	// Vec3s. The default 1024 constraints alone are 4.2 MB, so the old pool could not serve a single
	// allocation.
	mTempAllocator(new TempAllocatorImpl(32 * 1024 * 1024)),
#endif
	mJobSystem(new JobSystemThreadPool(cMaxPhysicsJobs, cMaxPhysicsBarriers, inWorkerThreads)),
	mDeltaTime(inDeltaTime),
	mCollisionSteps(inCollisionSteps)
{
	// Create physics system
	mSystem = new PhysicsSystem();
	mSystem->Init(inMaxBodies, 0, inMaxBodyPairs, inMaxContactConstraints, mBroadPhaseLayerInterface, mObjectVsBroadPhaseLayerFilter, mObjectVsObjectLayerFilter);
}

PhysicsTestContext::~PhysicsTestContext()
{
#ifdef JPH_DEBUG_RENDERER
	delete mDebugRenderer;
	delete mStreamWrapper;
	delete mStream;
#endif // JPH_DEBUG_RENDERER
	delete mSystem;
	delete mJobSystem;
	delete mTempAllocator;
}

void PhysicsTestContext::ZeroGravity()
{
	mSystem->SetGravity(Vec4::sZero());
}

Body &PhysicsTestContext::CreateFloor()
{
	BodyCreationSettings settings;
	// TODO(4D): the floor's boundary is a 3D volume, so the slab is wide in X, Z *and* W
	settings.SetShape(new BoxShape(Vec4(100.0f, 1.0f, 100.0f, 100.0f), 0.0f));
	settings.mPosition = RVec4(0.0f, -1.0f, 0.0f, 0.0f);
	settings.mMotionType = EMotionType::Static;
	settings.mObjectLayer = Layers::NON_MOVING;

	Body &floor = *mSystem->GetBodyInterface().CreateBody(settings);
	mSystem->GetBodyInterface().AddBody(floor.GetID(), EActivation::DontActivate);
	return floor;
}

Body &PhysicsTestContext::CreateBody(const BodyCreationSettings &inSettings, EActivation inActivation)
{
	Body &body = *mSystem->GetBodyInterface().CreateBody(inSettings);
	mSystem->GetBodyInterface().AddBody(body.GetID(), inActivation);
	return body;
}

Body &PhysicsTestContext::CreateBody(const ShapeSettings *inShapeSettings, RVec4Arg inPosition, RotorArg inRotation, EMotionType inMotionType, EMotionQuality inMotionQuality, ObjectLayer inLayer, EActivation inActivation)
{
	BodyCreationSettings settings;
	settings.SetShapeSettings(inShapeSettings);
	settings.mPosition = inPosition;
	settings.mRotation = inRotation;
	settings.mMotionType = inMotionType;
	settings.mMotionQuality = inMotionQuality;
	settings.mObjectLayer = inLayer;
	settings.mLinearDamping = 0.0f;
	settings.mAngularDamping = 0.0f;
	settings.mCollisionGroup.SetGroupID(0);
	return CreateBody(settings, inActivation);
}

Body &PhysicsTestContext::CreateBox(RVec4Arg inPosition, RotorArg inRotation, EMotionType inMotionType, EMotionQuality inMotionQuality, ObjectLayer inLayer, Vec4Arg inHalfExtent, EActivation inActivation)
{
	return CreateBody(new BoxShapeSettings(inHalfExtent), inPosition, inRotation, inMotionType, inMotionQuality, inLayer, inActivation);
}

Body &PhysicsTestContext::CreateSphere(RVec4Arg inPosition, float inRadius, EMotionType inMotionType, EMotionQuality inMotionQuality, ObjectLayer inLayer, EActivation inActivation)
{
	return CreateBody(new SphereShapeSettings(inRadius), inPosition, Rotor::sIdentity(), inMotionType, inMotionQuality, inLayer, inActivation);
}

EPhysicsUpdateError PhysicsTestContext::SimulateNoDeltaTime()
{
	EPhysicsUpdateError errors = mSystem->Update(0.0f, mCollisionSteps, mTempAllocator, mJobSystem);
#ifndef JPH_DISABLE_TEMP_ALLOCATOR
	JPH_ASSERT(static_cast<TempAllocatorImpl *>(mTempAllocator)->IsEmpty());
#endif // JPH_DISABLE_TEMP_ALLOCATOR
	return errors;
}

EPhysicsUpdateError PhysicsTestContext::SimulateSingleStep()
{
	EPhysicsUpdateError errors = mSystem->Update(mDeltaTime, mCollisionSteps, mTempAllocator, mJobSystem);
#ifndef JPH_DISABLE_TEMP_ALLOCATOR
	JPH_ASSERT(static_cast<TempAllocatorImpl *>(mTempAllocator)->IsEmpty());
#endif // JPH_DISABLE_TEMP_ALLOCATOR
#ifdef JPH_DEBUG_RENDERER
	if (mDebugRenderer != nullptr)
	{
		mSystem->DrawBodies(BodyManager::DrawSettings(), mDebugRenderer);
		mSystem->DrawConstraints(mDebugRenderer);
		mDebugRenderer->EndFrame();
	}
#endif // JPH_DEBUG_RENDERER
	return errors;
}

EPhysicsUpdateError PhysicsTestContext::Simulate(float inTotalTime, function<void()> inPreStepCallback)
{
	EPhysicsUpdateError errors = EPhysicsUpdateError::None;

	const int cNumSteps = int(round(inTotalTime / mDeltaTime));
	for (int s = 0; s < cNumSteps; ++s)
	{
		inPreStepCallback();
		errors |= SimulateSingleStep();
	}

	return errors;
}

RVec4 PhysicsTestContext::PredictPosition(RVec4Arg inPosition, Vec4Arg inVelocity, Vec4Arg inAcceleration, float inTotalTime) const
{
	// Integrate position using a Symplectic Euler step (just like the PhysicsSystem)
	RVec4 pos = inPosition;
	Vec4 vel = inVelocity;

	const float delta_time = GetStepDeltaTime();
	const int cNumSteps = int(round(inTotalTime / delta_time));
	for (int s = 0; s < cNumSteps; ++s)
	{
		vel += inAcceleration * delta_time;
		pos += vel * delta_time;
	}
	return pos;
}

// Predict rotation assuming ballistic motion using initial orientation, angular velocity angular acceleration and time
// TODO(4D): a 4D angular velocity is a bivector that need not describe a single rotation plane, so
// there is no axis/angle to feed to Quat::sRotation. This mirrors Body::AddRotationStep exactly -
// R <- exp(Omega * dt / 2) * R via the general rotor exponential - rather than reimplementing it,
// so a test comparing against this prediction cannot pass by making the same mistake twice in a
// different way. It is still an independent check of the *solver*: it asserts that an unconstrained
// body does nothing but integrate.
Rotor PhysicsTestContext::PredictOrientation(RotorArg inRotation, BivecArg inAngularVelocity, BivecArg inAngularAcceleration, float inTotalTime) const
{
	// Integrate orientation using a Symplectic Euler step (just like the PhysicsSystem)
	Rotor rot = inRotation;
	Bivec vel = inAngularVelocity;

	const float delta_time = GetStepDeltaTime();
	const int cNumSteps = int(round(inTotalTime / delta_time));
	for (int s = 0; s < cNumSteps; ++s)
	{
		vel += inAngularAcceleration * delta_time;
		rot = (Rotor::sExp(0.5f * (vel * delta_time)) * rot).Normalized();
	}
	return rot;
}

#ifdef JPH_DEBUG_RENDERER

void PhysicsTestContext::RecordDebugOutput(const char *inFileName)
{
	mStream = new ofstream;
	mStream->open(inFileName, ofstream::out | ofstream::binary | ofstream::trunc);
	if (mStream->is_open())
	{
		mStreamWrapper = new StreamOutWrapper(*mStream);
		mDebugRenderer = new DebugRendererRecorder(*mStreamWrapper);
	}
	else
	{
		delete mStream;
		mStream = nullptr;
	}
}

#endif // JPH_DEBUG_RENDERER
