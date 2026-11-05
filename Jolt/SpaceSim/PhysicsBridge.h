// SPDX-License-Identifier: MIT
#pragma once

#include <Jolt/Jolt.h>
#include <Jolt/Core/FPException.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Body/BodyLock.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

JPH_NAMESPACE_BEGIN
namespace SpaceSim
{
/// Explicit vacuum configuration. Collision response/friction remain material choices.
/// Velocity limits are declared safety limits, not conservation-test error correction.
inline bool ConfigureFreeSpaceBody(BodyCreationSettings &ioSettings, float inMaxLinearSpeed, float inMaxAngularSpeed)
{
    if (!std::isfinite(inMaxLinearSpeed) || !std::isfinite(inMaxAngularSpeed)
        || inMaxLinearSpeed <= 0 || inMaxAngularSpeed <= 0)
        return false;
    ioSettings.mLinearDamping = 0;
    ioSettings.mAngularDamping = 0;
    ioSettings.mGravityFactor = 0;
    ioSettings.mAllowSleeping = false;
    ioSettings.mApplyGyroscopicForce = true;
    ioSettings.mMaxLinearVelocity = inMaxLinearSpeed;
    ioSettings.mMaxAngularVelocity = inMaxAngularSpeed;
    return true;
}

/// An actual delivered impulse integral over ONE fixed interval, in the solve frame.
/// Angular impulse is about the body's current COM, using F wedge r. Aggregate all
/// components of a body before submitting it; a body may occur only once per batch.
struct WorldImpulseIntegral
{
    BodyID mBody;
    Vec4 mLinear = Vec4::sZero();
    Bivec mAngular = Bivec::sZero();
};

enum class BridgeResult { Success, InvalidConfiguration, InvalidDelivery, InvalidBody, PhysicsFailure, FailedStopped, NumericalFailure };
enum class BridgeNumericalError { None, LinearVelocityClamped, AngularVelocityClamped, NonFiniteState, LinearSpeedLimitExceeded, AngularSpeedLimitExceeded, ImpulseDeliveryMismatch };

/// Serial boundary adapter, not a second integrator. The coordinator must own the
/// world exclusively during Step, including validation, delivery and Update. No
/// other producer may have queued forces. Thread-safe Jolt jobs run INSIDE Update.
///
/// Integrated actuation is delivered as constant average force/torque. This preserves
/// its integral but approximates its time distribution/orientation over the tick;
/// it is NOT a within-tick depletion/TOI event integrator. Use refinement tests and
/// reject workloads needing unimplemented event splitting. Gravity must be configured
/// explicitly. An Update error, any core velocity clamp, or nonfinite body state
/// leaves the live core potentially partially advanced: fail-stop, do not publish
/// a snapshot, and restore a complete checkpoint before resuming. The tick only
/// advances on Success; numerical diagnostics identify the first affected body.
class PhysicsBridge
{
public:
    PhysicsBridge(PhysicsSystem &inSystem, TempAllocator &inAllocator, JobSystem &inJobs,
        float inFixedStep, uint32 inMaxDeliveries = 4096)
        : mSystem(inSystem), mAllocator(inAllocator), mJobs(inJobs), mFixedStep(inFixedStep),
          mMaxDeliveries(inMaxDeliveries)
    {
        // Declared adapter envelope: 1 microsecond through 1 second. Extremely
        // small intervals overflow reciprocal/squared-step terms in the core.
        mValid = std::isfinite(inFixedStep) && inFixedStep >= 1.0e-6f && inFixedStep <= 1.0f && inMaxDeliveries > 0
            && inMaxDeliveries <= 1048576;
        if (mValid)
        {
            mIDs.reserve(inMaxDeliveries);
            mPrepared.resize(inMaxDeliveries);
            mBodyIDs.reserve(mSystem.GetMaxBodies());
        }
    }

    PhysicsBridge(const PhysicsBridge &) = delete;
    PhysicsBridge &operator=(const PhysicsBridge &) = delete;
    PhysicsBridge(PhysicsBridge &&) = delete;
    PhysicsBridge &operator=(PhysicsBridge &&) = delete;

    BridgeResult Step(const WorldImpulseIntegral *inDeliveries, size_t inCount)
    {
        if (mFailed) return BridgeResult::FailedStopped;
        if (!mValid || mTick == std::numeric_limits<uint64>::max()) return BridgeResult::InvalidConfiguration;
        if (inCount > mMaxDeliveries || (inCount != 0 && inDeliveries == nullptr)) return BridgeResult::InvalidDelivery;
        mIDs.clear();
        for (size_t i = 0; i < inCount; ++i)
        {
            const auto &d = inDeliveries[i];
            // Compute once, validate the exact SIMD values, then submit those same
            // values. Scalar division and multiplication by a rounded reciprocal
            // can disagree at the finite/overflow boundary.
            for (int axis = 0; axis < 4; ++axis)
                if (!std::isfinite(d.mLinear[axis])) return BridgeResult::InvalidDelivery;
            for (int axis = 0; axis < 6; ++axis)
                if (!std::isfinite(d.mAngular[axis])) return BridgeResult::InvalidDelivery;
            auto &prepared = mPrepared[i];
            {
                // Finite integrals can overflow the exact SIMD conversion.
                // Mask only this validation arithmetic, never delivery/Update.
                [[maybe_unused]] FPExceptionDisableOverflow disable_overflow;
                [[maybe_unused]] FPExceptionDisableInvalid disable_invalid;
                prepared.mLinear = d.mLinear / mFixedStep;
                prepared.mAngular = d.mAngular * (1.0f / mFixedStep);
            }
            for (int axis = 0; axis < 4; ++axis)
                if (!std::isfinite(prepared.mLinear[axis])) return BridgeResult::InvalidDelivery;
            for (int axis = 0; axis < 6; ++axis)
                if (!std::isfinite(prepared.mAngular[axis])) return BridgeResult::InvalidDelivery;
            mIDs.push_back(d.mBody.GetIndexAndSequenceNumber());
        }
        std::sort(mIDs.begin(), mIDs.end());
        if (std::adjacent_find(mIDs.begin(), mIDs.end()) != mIDs.end()) return BridgeResult::InvalidDelivery;
        // Inspect every dynamic body before predictions, diagnostic reset or
        // Update. Active gyro integration can replace a NaN with zero, erasing
        // both the invalid input and the old sticky diagnostic.
        if (ValidateBoundary() != BridgeResult::Success) return BridgeResult::NumericalFailure;
        // Validate the ENTIRE batch before the first world mutation.
        for (size_t i = 0; i < inCount; ++i)
        {
            BodyLockRead lock(mSystem.GetBodyLockInterface(), inDeliveries[i].mBody);
            if (!lock.Succeeded() || !lock.GetBody().IsDynamic() || !lock.GetBody().IsInBroadPhase())
                return BridgeResult::InvalidBody;
            if (lock.GetBody().GetAccumulatedForce() != Vec4::sZero()
                || lock.GetBody().GetAccumulatedTorque() != Bivec::sZero())
                return BridgeResult::InvalidDelivery;
            const Body &body = lock.GetBody();
            const MotionProperties &motion = *body.GetMotionProperties();
            Vec4 linear;
            Bivec angular;
            {
                // Prediction is validation too: representable force/torque can
                // overflow when multiplied by inverse mass/inertia. Intermediate
                // infinities can then cause invalid SIMD arithmetic (0 * inf).
                [[maybe_unused]] FPExceptionDisableOverflow disable_overflow;
                [[maybe_unused]] FPExceptionDisableInvalid disable_invalid;
                linear = body.GetLinearVelocity() + mFixedStep * (mPrepared[i].mLinear * motion.GetInverseMass());
                angular = body.GetAngularVelocity() + mFixedStep * motion.MultiplyWorldSpaceInverseInertiaByVector(
                    body.GetRotation(), mPrepared[i].mAngular);
            }
            double linear_sq = 0, angular_sq = 0;
            for (int axis = 0; axis < 4; ++axis) linear_sq += double(linear[axis]) * double(linear[axis]);
            for (int axis = 0; axis < 6; ++axis) angular_sq += double(angular[axis]) * double(angular[axis]);
            if (!std::isfinite(linear_sq) || !std::isfinite(angular_sq)
                || linear_sq > double(motion.GetMaxLinearVelocity()) * motion.GetMaxLinearVelocity()
                || angular_sq > double(motion.GetMaxAngularVelocity()) * motion.GetMaxAngularVelocity())
                return BridgeResult::InvalidDelivery;
        }
        // Reject prequeued input anywhere in the owned world, even when this
        // batch is empty, before mutating recipients or clearing diagnostics.
        mSystem.GetBodies(mBodyIDs);
        for (BodyID id : mBodyIDs)
        {
            BodyLockRead lock(mSystem.GetBodyLockInterface(), id);
            if (lock.Succeeded() && lock.GetBody().IsDynamic()
                && (lock.GetBody().GetAccumulatedForce() != Vec4::sZero()
                    || lock.GetBody().GetAccumulatedTorque() != Bivec::sZero()))
                return BridgeResult::InvalidDelivery;
        }
        // Reset only once before this owned step. Diagnostics remain sticky
        // across every core clamp site, including pre-contact and CCD.
        for (BodyID id : mBodyIDs)
        {
            BodyLockWrite lock(mSystem.GetBodyLockInterface(), id);
            if (lock.Succeeded() && lock.GetBody().IsDynamic())
                lock.GetBody().GetMotionProperties()->ResetVelocityDiagnostics();
        }
        for (size_t i = 0; i < inCount; ++i)
            mSystem.GetBodyInterface().AddForceAndTorque(inDeliveries[i].mBody,
                mPrepared[i].mLinear, mPrepared[i].mAngular, EActivation::Activate);
        mLastError = mSystem.Update(mFixedStep, 1, &mAllocator, &mJobs);
        // Inspect ALL dynamic bodies, not only recipients: contacts/gyro/gravity
        // can clamp an undelivered body or leave an endpoint beyond its limits.
        if (ValidateBoundary() != BridgeResult::Success) return BridgeResult::NumericalFailure;
        if (mLastError != EPhysicsUpdateError::None)
        {
            mFailed = true;
            return BridgeResult::PhysicsFailure;
        }
        ++mTick;
        return BridgeResult::Success;
    }

    /// Validate immediately before publication as well as before Step. Sticky
    /// diagnostics are evidence of an unacknowledged numerical failure, NOT reset
    /// permission. Recovery requires restoring a complete world and a new bridge.
    BridgeResult ValidateBoundary()
    {
        if (mFailed) return BridgeResult::FailedStopped;
        mSystem.GetBodies(mBodyIDs);
        for (BodyID id : mBodyIDs)
        {
            BodyLockRead lock(mSystem.GetBodyLockInterface(), id);
            if (!lock.Succeeded() || !lock.GetBody().IsDynamic()) continue;
            const Body &body = lock.GetBody();
            const auto &motion = *body.GetMotionProperties();
            BridgeNumericalError error = BridgeNumericalError::None;
            if (motion.WasLinearVelocityClamped()) error = BridgeNumericalError::LinearVelocityClamped;
            if (motion.WasAngularVelocityClamped()) error = BridgeNumericalError::AngularVelocityClamped;
            // Check finite state before speed arithmetic and before any reset.
            if (motion.HadNonFiniteVelocity() || !IsFiniteState(body)) error = BridgeNumericalError::NonFiniteState;
            if (error == BridgeNumericalError::None)
            {
                if (!SpeedFits(body.GetLinearVelocity(), 4, motion.GetMaxLinearVelocity())) error = BridgeNumericalError::LinearSpeedLimitExceeded;
                if (!SpeedFits(body.GetAngularVelocity(), 6, motion.GetMaxAngularVelocity())) error = BridgeNumericalError::AngularSpeedLimitExceeded;
            }
            if (error != BridgeNumericalError::None)
            {
                mLastNumericalError = error;
                mLastNumericalBody = id;
                mFailed = true;
                return BridgeResult::NumericalFailure;
            }
        }
        return BridgeResult::Success;
    }

    /// Tick-end instantaneous query response; NOT a second physics integration.
    /// ImpactBatch owns event identity/exactly-once dispatch. Predict the exact core
    /// SIMD operations with the installed inverse mass/inertia, before either kick.
    /// InvalidDelivery leaves the body unchanged; any post-write discrepancy is a
    /// fail-stop numerical error and must never be published by the coordinator.
    BridgeResult ApplyBoundaryImpulse(const WorldImpulseIntegral &delivery)
    {
        if (mFailed) return BridgeResult::FailedStopped;
        if (!mValid) return BridgeResult::InvalidConfiguration;
        if (ValidateBoundary() != BridgeResult::Success) return BridgeResult::NumericalFailure;
        Vec4 linear;
        Bivec angular;
        {
            BodyLockRead lock(mSystem.GetBodyLockInterface(), delivery.mBody);
            if (!lock.Succeeded() || !lock.GetBody().IsDynamic() || !lock.GetBody().IsInBroadPhase()) return BridgeResult::InvalidBody;
            const Body &body = lock.GetBody();
            const auto &motion = *body.GetMotionProperties();
            for (uint i = 0; i < 4; ++i) if (!std::isfinite(delivery.mLinear[i])) return BridgeResult::InvalidDelivery;
            for (uint i = 0; i < 6; ++i) if (!std::isfinite(delivery.mAngular[i])) return BridgeResult::InvalidDelivery;
            {
                // Only pure preview arithmetic is masked. Never mask delivery or Update.
                [[maybe_unused]] FPExceptionDisableOverflow overflow;
                [[maybe_unused]] FPExceptionDisableInvalid invalid;
                linear = body.GetLinearVelocity() + delivery.mLinear * motion.GetInverseMass();
                angular = body.GetAngularVelocity() + motion.MultiplyWorldSpaceInverseInertiaByVector(body.GetRotation(), delivery.mAngular);
            }
            if (!SpeedFits(linear, 4, motion.GetMaxLinearVelocity()) || !SpeedFits(angular, 6, motion.GetMaxAngularVelocity())
                || !DeliveredImpulseFits(body, delivery, linear, angular)) return BridgeResult::InvalidDelivery;
        }
        auto &bi = mSystem.GetBodyInterface();
        bi.AddImpulse(delivery.mBody, delivery.mLinear);
        bi.AddAngularImpulse(delivery.mBody, delivery.mAngular);
        if (ValidateBoundary() != BridgeResult::Success) return BridgeResult::NumericalFailure;
        if (bi.GetLinearVelocity(delivery.mBody) != linear || bi.GetAngularVelocity(delivery.mBody) != angular)
        {
            mLastNumericalError = BridgeNumericalError::ImpulseDeliveryMismatch;
            mLastNumericalBody = delivery.mBody;
            mFailed = true;
            return BridgeResult::NumericalFailure;
        }
        return BridgeResult::Success;
    }

    uint64 GetTick() const { return mTick; }
    float GetFixedStep() const { return mFixedStep; }
    bool IsFailedStopped() const { return mFailed; }
    EPhysicsUpdateError GetLastPhysicsError() const { return mLastError; }
    BridgeNumericalError GetLastNumericalError() const { return mLastNumericalError; }
    BodyID GetLastNumericalBody() const { return mLastNumericalBody; }

private:
    // Declared impulse-ledger budget: 1e-5 absolute per lane, or 1e-5 relative
    // for larger kicks. Even below that error budget, a nonzero integral may
    // NEVER be reported delivered when the entire kick rounded away.
    static bool DeliveredImpulseFits(const Body &body, const WorldImpulseIntegral &delivery, Vec4Arg linear, BivecArg angular)
    {
        const auto &motion = *body.GetMotionProperties();
        double inverse_mass = motion.GetInverseMass();
        if (!std::isfinite(inverse_mass) || inverse_mass <= 0) return false;
        if ((delivery.mLinear != Vec4::sZero() && linear == body.GetLinearVelocity())
            || (delivery.mAngular != Bivec::sZero() && angular == body.GetAngularVelocity())) return false;
        for (uint i = 0; i < 4; ++i)
        {
            double actual = (double(linear[i]) - body.GetLinearVelocity()[i]) / inverse_mass;
            if (std::abs(actual - delivery.mLinear[i]) > 1e-5 * std::max(1.0, std::abs(double(delivery.mLinear[i])))) return false;
        }
        Rotor principal = body.GetRotation() * motion.GetInertiaRotation();
        Bivec local_delta = principal.InverseRotateBivec(angular - body.GetAngularVelocity());
        double actual[6] = {};
        for (uint i = 0; i < 6; ++i)
        {
            double inverse = motion.GetInverseInertiaDiagonal()[i];
            if (!std::isfinite(inverse) || inverse <= 0) return false;
            double momentum = double(local_delta[i]) / inverse;
            float lane[6] = {}; lane[i] = 1;
            Bivec column = principal.RotateBivec(Bivec(lane[0], lane[1], lane[2], lane[3], lane[4], lane[5]));
            for (uint j = 0; j < 6; ++j) actual[j] += double(column[j]) * momentum;
        }
        for (uint i = 0; i < 6; ++i)
            if (!std::isfinite(actual[i]) || std::abs(actual[i] - delivery.mAngular[i]) > 1e-5 * std::max(1.0, std::abs(double(delivery.mAngular[i])))) return false;
        return true;
    }
    template <class V> static bool SpeedFits(const V &velocity, uint lanes, float limit)
    {
        if (!std::isfinite(limit) || limit <= 0 || double(limit) * limit >= double(std::numeric_limits<float>::max())) return false;
        double squared = 0;
        for (uint i = 0; i < lanes; ++i)
        {
            if (!std::isfinite(velocity[i])) return false;
            squared += double(velocity[i]) * velocity[i];
        }
        // Core clamp arithmetic squares/reduces in float. A double norm can
        // admit a lane combination which the exact SIMD predicate will clamp.
        if (squared > double(limit) * limit || squared >= double(std::numeric_limits<float>::max())) return false;
        float native_squared;
        {
            [[maybe_unused]] FPExceptionDisableOverflow overflow;
            [[maybe_unused]] FPExceptionDisableInvalid invalid;
            native_squared = velocity.LengthSq();
        }
        return std::isfinite(native_squared) && native_squared <= limit * limit;
    }
    static bool IsFiniteState(const Body &inBody)
    {
        const auto position = inBody.GetCenterOfMassPosition();
        const Vec4 linear = inBody.GetLinearVelocity();
        const Bivec angular = inBody.GetAngularVelocity();
        const Vec8 rotation = inBody.GetRotation().GetVec8();
        for (int axis = 0; axis < 4; ++axis)
            if (!std::isfinite(position[axis]) || !std::isfinite(linear[axis])) return false;
        for (int axis = 0; axis < 6; ++axis)
            if (!std::isfinite(angular[axis])) return false;
        for (int axis = 0; axis < 8; ++axis)
            if (!std::isfinite(rotation[axis])) return false;
        return true;
    }
    PhysicsSystem &mSystem;
    TempAllocator &mAllocator;
    JobSystem &mJobs;
    float mFixedStep;
    uint32 mMaxDeliveries;
    std::vector<uint32> mIDs;
    struct PreparedDelivery { Vec4 mLinear; Bivec mAngular; };
    std::vector<PreparedDelivery> mPrepared;
    BodyIDVector mBodyIDs;
    uint64 mTick = 0;
    EPhysicsUpdateError mLastError = EPhysicsUpdateError::None;
    BridgeNumericalError mLastNumericalError = BridgeNumericalError::None;
    BodyID mLastNumericalBody;
    bool mValid = false;
    bool mFailed = false;
};
}
JPH_NAMESPACE_END
