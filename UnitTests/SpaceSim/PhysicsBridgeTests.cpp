// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Core/FPException.h>
#include <Jolt/SpaceSim/PhysicsBridge.h>
#include <type_traits>

TEST_SUITE("SpaceSimPhysicsBridge")
{
    TEST_CASE("BridgeOwnershipCannotBeCopiedOrMoved")
    {
        CHECK_FALSE(std::is_copy_constructible<JPH::SpaceSim::PhysicsBridge>::value);
        CHECK_FALSE(std::is_copy_assignable<JPH::SpaceSim::PhysicsBridge>::value);
        CHECK_FALSE(std::is_move_constructible<JPH::SpaceSim::PhysicsBridge>::value);
        CHECK_FALSE(std::is_move_assignable<JPH::SpaceSim::PhysicsBridge>::value);
    }

    TEST_CASE("PrequeuedForceOrTorqueOnAnyBodyRejectsTheEntireStep")
    {
        for (bool queue_on_recipient : { false, true })
            for (bool torque : { false, true })
                for (bool deliver : { false, true })
                {
                    PhysicsTestContext context;
                    context.ZeroGravity();
                    JPH::Body &recipient = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
                        JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
                    JPH::Body &other = context.CreateBox(JPH::RVec4(10, 0, 0, 0), JPH::Rotor::sIdentity(),
                        JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
                    JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
                        *context.GetJobSystem(), context.GetDeltaTime());
                    JPH::SpaceSim::WorldImpulseIntegral delivery;
                    delivery.mBody = recipient.GetID();
                    delivery.mLinear = JPH::Vec4::sAxisW();
                    JPH::Body &queued = queue_on_recipient ? recipient : other;
                    const JPH::Bivec queued_torque(0, 0, 0, 0, 0, 1);
                    if (torque) context.GetBodyInterface().AddTorque(queued.GetID(), queued_torque);
                    else context.GetBodyInterface().AddForce(queued.GetID(), JPH::Vec4::sAxisX());
                    CHECK(bridge.Step(deliver ? &delivery : nullptr, deliver ? 1 : 0) == JPH::SpaceSim::BridgeResult::InvalidDelivery);
                    CHECK(recipient.GetAccumulatedForce() == (!torque && queue_on_recipient ? JPH::Vec4::sAxisX() : JPH::Vec4::sZero()));
                    CHECK(recipient.GetLinearVelocity() == JPH::Vec4::sZero());
                    CHECK(queued.GetAccumulatedForce() == (torque ? JPH::Vec4::sZero() : JPH::Vec4::sAxisX()));
                    CHECK(queued.GetAccumulatedTorque() == (torque ? queued_torque : JPH::Bivec::sZero()));
                    CHECK(bridge.GetTick() == 0);
                    CHECK_FALSE(bridge.IsFailedStopped());
                }
    }

    TEST_CASE("FreeSpacePresetAndIntegratedDeliveryUseTheRealPhysicsStep")
    {
        PhysicsTestContext context(1.0f / 120.0f);
        context.ZeroGravity();
        JPH::BodyCreationSettings settings(new JPH::BoxShape(JPH::Vec4::sReplicate(0.5f)),
            JPH::RVec4::sZero(), JPH::Rotor::sIdentity(), JPH::EMotionType::Dynamic, Layers::MOVING);
        CHECK(JPH::SpaceSim::ConfigureFreeSpaceBody(settings, 1000.0f, 100.0f));
        CHECK(settings.mLinearDamping == 0.0f);
        CHECK(settings.mAngularDamping == 0.0f);
        CHECK_FALSE(settings.mAllowSleeping);
        CHECK(settings.mApplyGyroscopicForce);
        JPH::Body &body = context.CreateBody(settings, JPH::EActivation::Activate);
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime());
        JPH::SpaceSim::WorldImpulseIntegral delivery;
        delivery.mBody = body.GetID();
        delivery.mLinear = JPH::Vec4(0, 0, 0, 2);
        float inverse_mass = body.GetMotionProperties()->GetInverseMass();
        CHECK(bridge.Step(&delivery, 1) == JPH::SpaceSim::BridgeResult::Success);
        CHECK(body.GetLinearVelocity().GetW() == doctest::Approx(2.0f * inverse_mass));
        CHECK(bridge.GetTick() == 1);
        CHECK(body.GetAccumulatedForce() == JPH::Vec4::sZero());
        CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::Success);
        CHECK(body.GetLinearVelocity().GetW() == doctest::Approx(2.0f * inverse_mass));
        CHECK(bridge.GetTick() == 2);
    }

    TEST_CASE("ReciprocalRoundedAngularOverflowIsRejectedBeforeMutation")
    {
        PhysicsTestContext context(0.10618418455123901f);
        context.ZeroGravity();
        JPH::BodyCreationSettings settings(new JPH::BoxShape(JPH::Vec4::sReplicate(0.5f)),
            JPH::RVec4::sZero(), JPH::Rotor::sIdentity(), JPH::EMotionType::Dynamic, Layers::MOVING);
        CHECK(JPH::SpaceSim::ConfigureFreeSpaceBody(settings, 1000.0f, 1.0e9f));
        JPH::Body &body = context.CreateBody(settings, JPH::EActivation::Activate);
        body.GetMotionProperties()->SetInverseInertia(JPH::Bivec(1e-30f, 1e-30f, 1e-30f, 1e-30f, 1e-30f, 1e-30f), JPH::Rotor::sIdentity());
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime());
        JPH::SpaceSim::WorldImpulseIntegral delivery;
        delivery.mBody = body.GetID();
        delivery.mAngular = JPH::Bivec(3.613260311336416e37f, 0, 0, 0, 0, 0);
        {
            // This fixture deliberately evaluates the rounded overflow. Restore
            // native traps before calling Step: validation must protect itself.
            [[maybe_unused]] JPH::FPExceptionDisableOverflow disable_overflow;
            CHECK(std::isfinite(delivery.mAngular[0] / context.GetDeltaTime()));
            CHECK_FALSE(std::isfinite((delivery.mAngular * (1.0f / context.GetDeltaTime()))[0]));
        }
        CHECK(bridge.Step(&delivery, 1) == JPH::SpaceSim::BridgeResult::InvalidDelivery);
        CHECK(body.GetAngularVelocity() == JPH::Bivec::sZero());
        CHECK(body.GetAccumulatedTorque() == JPH::Bivec::sZero());
        CHECK(body.GetRotation() == JPH::Rotor::sIdentity());
        CHECK(bridge.GetTick() == 0);
        CHECK_FALSE(bridge.IsFailedStopped());
    }

    TEST_CASE("FinitePreparedDeliveryWithOverflowingPredictionRejectsBeforeMutation")
    {
        for (bool angular : { false, true })
        {
            CAPTURE(angular);
            PhysicsTestContext context;
            context.ZeroGravity();
            JPH::Body &body = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
                JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
            body.GetMotionProperties()->SetInverseMass(1.0e30f);
            body.GetMotionProperties()->SetInverseInertia(JPH::Bivec(1.0e30f, 1.0e30f, 1.0e30f, 1.0e30f, 1.0e30f, 1.0e30f), JPH::Rotor::sIdentity());
            JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
                *context.GetJobSystem(), context.GetDeltaTime());
            JPH::SpaceSim::WorldImpulseIntegral delivery;
            delivery.mBody = body.GetID();
            if (angular) delivery.mAngular = JPH::Bivec(1.0e10f, 0, 0, 0, 0, 0);
            else delivery.mLinear = JPH::Vec4(1.0e10f, 0, 0, 0);
            CHECK(bridge.Step(&delivery, 1) == JPH::SpaceSim::BridgeResult::InvalidDelivery);
            CHECK(body.GetAccumulatedForce() == JPH::Vec4::sZero());
            CHECK(body.GetAccumulatedTorque() == JPH::Bivec::sZero());
            CHECK(body.GetLinearVelocity() == JPH::Vec4::sZero());
            CHECK(body.GetAngularVelocity() == JPH::Bivec::sZero());
            CHECK(bridge.GetTick() == 0);
            CHECK_FALSE(bridge.IsFailedStopped());
        }
    }

    TEST_CASE("GyroscopicClampFailStopsEvenWithoutAnyDelivery")
    {
        for (bool deliver : { false, true })
        {
            PhysicsTestContext context(0.1f);
            context.ZeroGravity();
            JPH::BodyCreationSettings settings(new JPH::BoxShape(JPH::Vec4(3.0f, 0.5f, 1.0f, 2.0f)),
                JPH::RVec4::sZero(), JPH::Rotor::sIdentity(), JPH::EMotionType::Dynamic, Layers::MOVING);
            CHECK(JPH::SpaceSim::ConfigureFreeSpaceBody(settings, 1000.0f, 1.4143f));
            settings.mAngularVelocity = JPH::Bivec(1, 1, 0, 0, 0, 0);
            JPH::Body &body = context.CreateBody(settings, JPH::EActivation::Activate);
            JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
                *context.GetJobSystem(), context.GetDeltaTime());
            JPH::SpaceSim::WorldImpulseIntegral delivery;
            delivery.mBody = body.GetID();
            CHECK(bridge.Step(deliver ? &delivery : nullptr, deliver ? 1 : 0) == JPH::SpaceSim::BridgeResult::NumericalFailure);
            CHECK(bridge.GetLastNumericalError() == JPH::SpaceSim::BridgeNumericalError::AngularVelocityClamped);
            CHECK(bridge.GetLastNumericalBody() == body.GetID());
            CHECK(bridge.IsFailedStopped());
            CHECK(bridge.GetTick() == 0);
            CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::FailedStopped);
        }
    }

    TEST_CASE("DeliveryBeyondTheDeclaredSpeedLimitIsRejectedBeforeMutation")
    {
        PhysicsTestContext context(1.0f / 120.0f);
        context.ZeroGravity();
        JPH::Body &body = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
            JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING,
            JPH::Vec4::sReplicate(0.5f));
        body.GetMotionProperties()->SetMaxLinearVelocity(1.0f);
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime());
        JPH::SpaceSim::WorldImpulseIntegral delivery;
        delivery.mBody = body.GetID();
        delivery.mLinear = JPH::Vec4::sAxisW() * (2.0f / body.GetMotionProperties()->GetInverseMass());
        CHECK(bridge.Step(&delivery, 1) == JPH::SpaceSim::BridgeResult::InvalidDelivery);
        CHECK(body.GetLinearVelocity() == JPH::Vec4::sZero());
        CHECK(bridge.GetTick() == 0);
    }

    TEST_CASE("UnrepresentableStepIsRejectedWithoutAdvancing")
    {
        PhysicsTestContext context;
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), std::numeric_limits<float>::min());
        CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::InvalidConfiguration);
        CHECK(bridge.GetTick() == 0);
    }

    TEST_CASE("AllSixAngularPlanesAreDeliveredExactlyOnce")
    {
        PhysicsTestContext context(1.0f / 120.0f);
        context.ZeroGravity();
        JPH::BodyCreationSettings settings(new JPH::BoxShape(JPH::Vec4::sReplicate(0.5f)),
            JPH::RVec4::sZero(), JPH::Rotor::sIdentity(), JPH::EMotionType::Dynamic, Layers::MOVING);
        CHECK(JPH::SpaceSim::ConfigureFreeSpaceBody(settings, 1000.0f, 100.0f));
        JPH::Body &body = context.CreateBody(settings, JPH::EActivation::Activate);
        body.GetMotionProperties()->SetInverseInertia(JPH::Bivec(1, 1, 1, 1, 1, 1), JPH::Rotor::sIdentity());
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime(), 1);
        JPH::SpaceSim::WorldImpulseIntegral delivery;
        delivery.mBody = body.GetID();
        delivery.mAngular = JPH::Bivec(0.01f, -0.02f, 0.03f, -0.04f, 0.05f, -0.06f);
        CHECK(bridge.Step(&delivery, 1) == JPH::SpaceSim::BridgeResult::Success);
        for (int axis = 0; axis < 6; ++axis)
            CHECK(body.GetAngularVelocity()[axis] == doctest::Approx(delivery.mAngular[axis]).epsilon(1e-5));
        CHECK(body.GetAccumulatedTorque() == JPH::Bivec::sZero());
        CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::Success);
        for (int axis = 0; axis < 6; ++axis)
            CHECK(body.GetAngularVelocity()[axis] == doctest::Approx(delivery.mAngular[axis]).epsilon(1e-5));
        CHECK(bridge.GetTick() == 2);
    }

    TEST_CASE("EveryNonFiniteOrOverflowingDeliveryLaneRejectsTheWholeBatch")
    {
        PhysicsTestContext context(0.10618418455123901f);
        context.ZeroGravity();
        JPH::Body &first = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
            JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
        JPH::Body &second = context.CreateBox(JPH::RVec4(10, 0, 0, 0), JPH::Rotor::sIdentity(),
            JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime(), 2);
        for (float invalid : { std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity(),
            -std::numeric_limits<float>::infinity(), std::numeric_limits<float>::max(), 3.613260311336416e37f })
            for (int lane = 0; lane < 10; ++lane)
            {
                JPH::SpaceSim::WorldImpulseIntegral batch[2];
                batch[0].mBody = first.GetID();
                batch[0].mLinear = JPH::Vec4::sAxisW();
                batch[1].mBody = second.GetID();
                if (lane < 4)
                {
                    float value[4] = {};
                    value[lane] = invalid;
                    batch[1].mLinear = JPH::Vec4(value[0], value[1], value[2], value[3]);
                }
                else
                {
                    float value[6] = {};
                    value[lane - 4] = invalid;
                    batch[1].mAngular = JPH::Bivec(value[0], value[1], value[2], value[3], value[4], value[5]);
                }
                CHECK(bridge.Step(batch, 2) == JPH::SpaceSim::BridgeResult::InvalidDelivery);
                CHECK(first.GetAccumulatedForce() == JPH::Vec4::sZero());
                CHECK(second.GetAccumulatedTorque() == JPH::Bivec::sZero());
                CHECK(first.GetLinearVelocity() == JPH::Vec4::sZero());
                CHECK(second.GetAngularVelocity() == JPH::Bivec::sZero());
                CHECK(bridge.GetTick() == 0);
            }
    }

    TEST_CASE("CapacityNullAndStaleIDsRejectWithoutPartialDelivery")
    {
        PhysicsTestContext context;
        context.ZeroGravity();
        JPH::Body &body = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
            JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
        JPH::Body &doomed = context.CreateBox(JPH::RVec4(10, 0, 0, 0), JPH::Rotor::sIdentity(),
            JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
        JPH::BodyID stale = doomed.GetID();
        context.GetBodyInterface().RemoveBody(stale);
        context.GetBodyInterface().DestroyBody(stale);
        JPH::SpaceSim::WorldImpulseIntegral batch[2];
        batch[0].mBody = body.GetID();
        batch[0].mLinear = JPH::Vec4::sAxisW();
        batch[1].mBody = stale;
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime(), 2);
        CHECK(bridge.Step(nullptr, 1) == JPH::SpaceSim::BridgeResult::InvalidDelivery);
        CHECK(bridge.Step(batch, 3) == JPH::SpaceSim::BridgeResult::InvalidDelivery);
        CHECK(bridge.Step(batch, 2) == JPH::SpaceSim::BridgeResult::InvalidBody);
        CHECK(body.GetAccumulatedForce() == JPH::Vec4::sZero());
        CHECK(body.GetLinearVelocity() == JPH::Vec4::sZero());
        CHECK(bridge.GetTick() == 0);
        CHECK(bridge.Step(batch, 1) == JPH::SpaceSim::BridgeResult::Success);
    }

    TEST_CASE("ClampDiagnosticsAreStickyButResetBeforeTheNextOwnedStep")
    {
        PhysicsTestContext context;
        context.ZeroGravity();
        JPH::Body &body = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
            JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
        JPH::MotionProperties &motion = *body.GetMotionProperties();
        motion.SetMaxLinearVelocity(1);
        motion.SetMaxAngularVelocity(1);
        motion.SetLinearVelocityClamped(JPH::Vec4(2, 0, 0, 0));
        motion.SetAngularVelocityClamped(JPH::Bivec(0, 0, 0, 0, 0, 2));
        motion.SetLinearVelocityClamped(JPH::Vec4::sZero());
        motion.SetAngularVelocityClamped(JPH::Bivec::sZero());
        CHECK(motion.WasLinearVelocityClamped());
        CHECK(motion.WasAngularVelocityClamped());
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime());
        CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::Success);
        CHECK_FALSE(motion.WasLinearVelocityClamped());
        CHECK_FALSE(motion.WasAngularVelocityClamped());
        CHECK_FALSE(motion.HadNonFiniteVelocity());
        CHECK(bridge.GetLastNumericalError() == JPH::SpaceSim::BridgeNumericalError::None);
    }

    TEST_CASE("GravityClampStillFailStopsAfterContactReducesTheSpeed")
    {
        PhysicsTestContext context(0.1f);
        context.CreateFloor();
        JPH::BodyCreationSettings settings(new JPH::BoxShape(JPH::Vec4::sReplicate(0.5f)),
            JPH::RVec4(0, 0.49f, 0, 0), JPH::Rotor::sIdentity(), JPH::EMotionType::Dynamic, Layers::MOVING);
        settings.mMaxLinearVelocity = 0.1f;
        settings.mAllowSleeping = false;
        settings.mLinearDamping = 0;
        JPH::Body &body = context.CreateBody(settings, JPH::EActivation::Activate);
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime());
        CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::NumericalFailure);
        CHECK(body.GetLinearVelocity().Length() < 0.099f);
        CHECK(body.GetMotionProperties()->WasLinearVelocityClamped());
        CHECK(bridge.GetLastNumericalError() == JPH::SpaceSim::BridgeNumericalError::LinearVelocityClamped);
        CHECK(bridge.GetLastNumericalBody() == body.GetID());
        CHECK(bridge.GetTick() == 0);
        CHECK(bridge.IsFailedStopped());
    }

    TEST_CASE("TransientNonFiniteVelocityCannotBeHiddenByClampingToZero")
    {
        PhysicsTestContext context(0.1f);
        context.GetSystem()->SetGravity(JPH::Vec4(1.0e30f, 0, 0, 0));
        JPH::Body &body = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
            JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime());
        JPH::SpaceSim::BridgeResult result;
        {
            // Unlike delivery validation, this test deliberately overflows the
            // core's float speed-squared calculation to exercise its diagnostic.
            [[maybe_unused]] JPH::FPExceptionDisableOverflow disable_overflow;
            result = bridge.Step(nullptr, 0);
        }
        CHECK(result == JPH::SpaceSim::BridgeResult::NumericalFailure);
        CHECK(body.GetLinearVelocity() == JPH::Vec4::sZero());
        CHECK(body.GetMotionProperties()->HadNonFiniteVelocity());
        CHECK(bridge.GetLastNumericalError() == JPH::SpaceSim::BridgeNumericalError::NonFiniteState);
        CHECK(bridge.GetTick() == 0);
        CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::FailedStopped);
    }

    TEST_CASE("InactiveUndeliveredNonFiniteAngularStateFailStopsInEveryPlane")
    {
        for (int axis = 0; axis < 6; ++axis)
        {
            PhysicsTestContext context;
            context.ZeroGravity();
            JPH::Body &body = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
                JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING,
                JPH::Vec4::sReplicate(0.5f), JPH::EActivation::DontActivate);
            float value[6] = {};
            value[axis] = std::numeric_limits<float>::quiet_NaN();
            {
                [[maybe_unused]] JPH::FPExceptionDisableInvalid disable_invalid;
                body.GetMotionProperties()->SetAngularVelocityClamped(JPH::Bivec(value[0], value[1], value[2], value[3], value[4], value[5]));
            }
            JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
                *context.GetJobSystem(), context.GetDeltaTime());
            CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::NumericalFailure);
            CHECK(bridge.GetLastNumericalError() == JPH::SpaceSim::BridgeNumericalError::NonFiniteState);
            CHECK(bridge.GetLastNumericalBody() == body.GetID());
            CHECK(bridge.IsFailedStopped());
            CHECK(bridge.GetTick() == 0);
            CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::FailedStopped);
        }
    }

    TEST_CASE("ActiveGyroscopicNonFiniteAngularStateFailStopsBeforeResetInEveryPlane")
    {
        for (int axis = 0; axis < 6; ++axis)
            for (bool deliver : { false, true })
            {
                CAPTURE(axis);
                CAPTURE(deliver);
                PhysicsTestContext context;
                context.ZeroGravity();
                JPH::BodyCreationSettings settings(new JPH::BoxShape(JPH::Vec4(3, 0.5f, 1, 2)),
                    JPH::RVec4::sZero(), JPH::Rotor::sIdentity(), JPH::EMotionType::Dynamic, Layers::MOVING);
                CHECK(JPH::SpaceSim::ConfigureFreeSpaceBody(settings, 1000, 100));
                JPH::Body &body = context.CreateBody(settings, JPH::EActivation::Activate);
                JPH::Body &recipient = context.CreateBox(JPH::RVec4(10, 0, 0, 0), JPH::Rotor::sIdentity(),
                    JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING, JPH::Vec4::sReplicate(0.5f));
                float value[6] = {};
                value[axis] = std::numeric_limits<float>::quiet_NaN();
                {
                    // Only the intentional invalid-state injection is masked.
                    [[maybe_unused]] JPH::FPExceptionDisableInvalid disable_invalid;
                    body.GetMotionProperties()->SetAngularVelocityClamped(JPH::Bivec(value[0], value[1], value[2], value[3], value[4], value[5]));
                }
                CHECK(body.IsActive());
                CHECK(body.GetMotionProperties()->HadNonFiniteVelocity());
                JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
                    *context.GetJobSystem(), context.GetDeltaTime());
                JPH::SpaceSim::WorldImpulseIntegral delivery;
                delivery.mBody = recipient.GetID();
                delivery.mLinear = JPH::Vec4::sAxisW();
                CHECK(bridge.Step(deliver ? &delivery : nullptr, deliver ? 1 : 0) == JPH::SpaceSim::BridgeResult::NumericalFailure);
                CHECK(bridge.GetLastNumericalError() == JPH::SpaceSim::BridgeNumericalError::NonFiniteState);
                CHECK(bridge.GetLastNumericalBody() == body.GetID());
                CHECK(body.GetMotionProperties()->HadNonFiniteVelocity());
                CHECK_FALSE(std::isfinite(body.GetAngularVelocity()[axis]));
                CHECK(body.GetRotation() == JPH::Rotor::sIdentity());
                CHECK(recipient.GetAccumulatedForce() == JPH::Vec4::sZero());
                CHECK(recipient.GetLinearVelocity() == JPH::Vec4::sZero());
                CHECK(bridge.IsFailedStopped());
                CHECK(bridge.GetTick() == 0);
                CHECK(bridge.Step(nullptr, 0) == JPH::SpaceSim::BridgeResult::FailedStopped);
            }
    }

    TEST_CASE("InvalidAndDuplicateBatchesDoNotPartiallyApply")
    {
        PhysicsTestContext context(1.0f / 120.0f);
        context.ZeroGravity();
        JPH::Body &body = context.CreateBox(JPH::RVec4::sZero(), JPH::Rotor::sIdentity(),
            JPH::EMotionType::Dynamic, JPH::EMotionQuality::Discrete, Layers::MOVING,
            JPH::Vec4::sReplicate(0.5f));
        JPH::SpaceSim::PhysicsBridge bridge(*context.GetSystem(), *context.GetTempAllocator(),
            *context.GetJobSystem(), context.GetDeltaTime());
        JPH::SpaceSim::WorldImpulseIntegral batch[2];
        batch[0].mBody = body.GetID();
        batch[0].mLinear = JPH::Vec4::sAxisX();
        batch[1] = batch[0];
        CHECK(bridge.Step(batch, 2) == JPH::SpaceSim::BridgeResult::InvalidDelivery);
        CHECK(body.GetLinearVelocity() == JPH::Vec4::sZero());
        CHECK(body.GetAccumulatedForce() == JPH::Vec4::sZero());
        CHECK(bridge.GetTick() == 0);
        batch[1].mBody = JPH::BodyID();
        CHECK(bridge.Step(batch, 2) == JPH::SpaceSim::BridgeResult::InvalidBody);
        CHECK(body.GetAccumulatedForce() == JPH::Vec4::sZero());
        CHECK(bridge.GetTick() == 0);
    }
}
