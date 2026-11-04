// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Physics/Collision/ContactListener.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>

TEST_SUITE("CCDTempAllocatorTests")
{
	static void sLinearCastWithDefaultAllocator(bool inGyroscopicForce)
	{
		struct Listener final : ContactListener
		{
			int mContacts = 0;
			void OnContactAdded(const Body &, const Body &, const ContactManifold &, ContactSettings &) override
			{
				++mContacts;
			}
		} listener;

		PhysicsTestContext context(1.0f / 120.0f);
		context.ZeroGravity();
		context.GetSystem()->SetContactListener(&listener);
		BodyCreationSettings wall(new BoxShape(Vec4(0.05f, 20, 20, 20), 0), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Static, Layers::NON_MOVING);
		wall.mFriction = 0;
		wall.mRestitution = 0.1f;
		context.CreateBody(wall, EActivation::DontActivate);

		BodyCreationSettings settings(new BoxShape(Vec4(0.35f, 0.65f, 1.05f, 1.6f), 0), RVec4(-4, 0, 0, 0),
			Rotor::sExp(Bivec(0.2f, -0.1f, 0.15f, 0.12f, -0.08f, 0.18f)), EMotionType::Dynamic, Layers::MOVING);
		settings.mMotionQuality = EMotionQuality::LinearCast;
		settings.mApplyGyroscopicForce = inGyroscopicForce;
		settings.mLinearVelocity = Vec4(800, 0, 0, 0);
		settings.mAngularVelocity = Bivec(1.1f, -0.7f, 0.9f, 0.6f, -1.2f, 0.8f);
		settings.mLinearDamping = settings.mAngularDamping = 0;
		settings.mAllowSleeping = false;
		settings.mMaxLinearVelocity = settings.mMaxAngularVelocity = 1.0e6f;
		settings.mFriction = 0;
		settings.mRestitution = 0.1f;
		settings.mOverrideMassProperties = EOverrideMassProperties::CalculateInertia;
		settings.mMassPropertiesOverride.mMass = 1;
		Body &body = context.CreateBody(settings, EActivation::Activate);

		// With one active LinearCast body, the old 16-byte stack stride put the
		// CCDBody after an odd 16-byte chunk. Its ContactSettings contains Bivec,
		// requiring 32-byte alignment on AVX: placement-new crashed on a YMM store.
		// Use the real default allocator directly, never an aligned test substitute.
		CHECK(context.SimulateSingleStep() == EPhysicsUpdateError::None);
		CHECK(listener.mContacts == 1);
		CHECK(body.GetPosition().GetX() < 0);
		CHECK(body.GetLinearVelocity().GetX() < 700);
		CHECK_FALSE(body.GetAngularVelocity().IsNaN());
		CHECK_FALSE(body.GetRotation().IsNaN());
		CHECK(body.GetRotation().IsNormalized());
#ifndef JPH_DISABLE_TEMP_ALLOCATOR
		CHECK(static_cast<TempAllocatorImpl *>(context.GetTempAllocator())->IsEmpty());
#endif
		context.GetSystem()->SetContactListener(nullptr);
	}

	TEST_CASE("DefaultAllocatorLinearCastWithGyroscopicForce")
	{
		sLinearCastWithDefaultAllocator(true);
	}

	TEST_CASE("DefaultAllocatorLinearCastWithoutGyroscopicForce")
	{
		sLinearCastWithDefaultAllocator(false);
	}
}
