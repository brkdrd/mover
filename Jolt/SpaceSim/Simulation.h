// SPDX-License-Identifier: MIT
#pragma once
#include <Jolt/Core/FPException.h>
#include <Jolt/SpaceSim/MassModel.h>
#include <Jolt/SpaceSim/Combat.h>
#include <Jolt/SpaceSim/Propulsion.h>
#include <Jolt/SpaceSim/Authority.h>
#include <Jolt/SpaceSim/World.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Body/BodyLock.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/StaticCompoundShape.h>
#include <Jolt/Physics/Collision/Shape/OffsetCenterOfMassShape.h>

namespace JPH::SpaceSim
{
inline Vec4 PhysicsVector(const MassVector &v) { return Vec4(float(v[0]), float(v[1]), float(v[2]), float(v[3])); }
inline RVec4 PhysicsPosition(const MassVector &v) { return RVec4(Real(v[0]), Real(v[1]), Real(v[2]), Real(v[3])); }
inline Bivec PhysicsBivector(const MassBivector &v) { return Bivec(float(v[0]), float(v[1]), float(v[2]), float(v[3]), float(v[4]), float(v[5])); }

inline constexpr double cPhysicsNarrowingBudget = 1e-5; ///< Absolute position / linear speed error.
inline constexpr double cPhysicsAngularNarrowingBudget = 1e-6;
template <class Scalar> inline bool PhysicsScalarRepresentable(double value, double budget)
{
    // Validate BEFORE narrowing: an overflowing conversion traps on native builds.
    if (!std::isfinite(value) || std::abs(value) > double(std::numeric_limits<Scalar>::max())) return false;
    return std::abs(double(Scalar(value)) - value) <= budget;
}
inline bool PhysicsRotorLanesBounded(RotorArg rotor)
{
    // Necessary for either unit-rotor tolerance used below, without evaluating
    // a float norm first. Leave the existing normalization predicates unchanged.
    for (uint i = 0; i < 8; ++i)
        if (!std::isfinite(rotor.GetVec8()[i]) || std::abs(rotor.GetVec8()[i]) > 2) return false;
    return true;
}
inline bool PhysicsStateRepresentable(const MassBody &state)
{
    if (!PhysicsRotorLanesBounded(state.mRotation)) return false;
    for (const auto &component : state.mComponents)
        if (!PhysicsRotorLanesBounded(component.mLocalRotation)) return false;
    for (uint i = 0; i < 4; ++i)
        // Double positions still belong to a bounded local solve frame. Merely
        // fitting Real is not sufficient for geometry/broad-phase arithmetic.
        if (!PhysicsScalarRepresentable<Real>(state.mOrigin[i], cPhysicsNarrowingBudget)
            || std::abs(state.mOrigin[i]) > cMaxLocalRange
            || !PhysicsScalarRepresentable<float>(state.mLinearVelocity[i], cPhysicsNarrowingBudget)) return false;
    for (uint i = 0; i < 6; ++i)
        if (!PhysicsScalarRepresentable<float>(state.mAngularVelocity[i], cPhysicsAngularNarrowingBudget)) return false;
    return true;
}

/// Prepared immutable boundary payload. No live world writes occur during preparation.
/// Hull assets are authoritative sharp boxes in the fixed geometry-origin frame;
/// mass-only components (e.g. point fuel reservoirs) need not have a collider.
class PreparedPhysicsBody
{
public:
    PreparedPhysicsBody() = default; ///< Empty output slot; cannot be installed.
    const MassBody &GetState() const { return mState; }
    const MassProperties &GetMass() const { return mMass; }
    const Shape *GetShape() const { return mShape.GetPtr(); }
private:
    friend bool PreparePhysicsBody(const MassBody &, const std::vector<CombatBox> &, PreparedPhysicsBody &);
    MassBody mState;
    MassProperties mMass;
    RefConst<Shape> mShape;
};

inline bool PreparePhysicsBody(const MassBody &state, const std::vector<CombatBox> &boxes, PreparedPhysicsBody &out)
{
    BodyMassData data;
    PreparedPhysicsBody result;
    if (!PhysicsStateRepresentable(state) || !MassModel::Derive(state, data) || !MassModel::ToMassProperties(data.mLocal, result.mMass)
        || boxes.empty() || boxes.size() > 256) return false;
    StaticCompoundShapeSettings compound;
    std::vector<uint64> ids;
    for (const auto &box : boxes)
    {
        if (!CombatDetail::Finite(box.mLocalCenter) || !CombatDetail::Finite(box.mHalfExtent)
            || !PhysicsRotorLanesBounded(box.mLocalRotation)) return false;
        // The existing radius budget implies these lane bounds. Check them before
        // float Length / rotor norms so finite hostile inputs cannot overflow.
        for (uint i = 0; i < 4; ++i)
            if (std::abs(box.mLocalCenter[i]) > 256 || std::abs(box.mHalfExtent[i]) > 256) return false;
        if (!CombatDetail::ValidRotor(box.mLocalRotation)
            || box.mLocalCenter.Length() + box.mHalfExtent.Length() > 256) return false;
        for (uint i = 0; i < 4; ++i) if (box.mHalfExtent[i] < 0.005f) return false;
        if (std::find(ids.begin(), ids.end(), box.mComponentID) != ids.end()
            || std::none_of(state.mComponents.begin(), state.mComponents.end(), [&](const ComponentMass &c) { return c.mID == box.mComponentID; })) return false;
        ids.push_back(box.mComponentID);
        Ref<Shape> leaf = new BoxShape(box.mHalfExtent, 0);
        leaf->SetUserData(box.mComponentID); // Full 64-bit semantic ID, never a truncating compound index.
        compound.AddShape(box.mLocalCenter, box.mLocalRotation, leaf);
    }
    auto built = compound.Create();
    if (built.HasError()) return false;
    RefConst<Shape> inner = built.Get();
    result.mShape = new OffsetCenterOfMassShape(inner, PhysicsVector(data.mLocal.mLocalCOM) - inner->GetCenterOfMass());
    for (uint i = 0; i < 4; ++i)
        if (std::abs(double(result.mShape->GetCenterOfMass()[i]) - data.mLocal.mLocalCOM[i]) > cPhysicsNarrowingBudget) return false;
    result.mState = state;
    out = std::move(result);
    return true;
}

inline bool PhysicsVelocityFits(const MassBody &state, float max_linear, float max_angular)
{
    if (!std::isfinite(max_linear) || !std::isfinite(max_angular) || max_linear <= 0 || max_angular <= 0
        || double(max_linear) * max_linear >= double(std::numeric_limits<float>::max())
        || double(max_angular) * max_angular >= double(std::numeric_limits<float>::max())) return false;
    // Divide only bounded finite lanes, never square untrusted doubles (DBL_MAX).
    auto fits = [](const auto &lanes, double limit)
    {
        double squared = 0;
        for (double v : lanes)
        {
            if (!std::isfinite(v) || std::abs(v) > limit) return false;
            double fraction = v / limit;
            squared += fraction * fraction;
        }
        return squared <= 1;
    };
    if (!fits(state.mLinearVelocity, max_linear) || !fits(state.mAngularVelocity, max_angular)) return false;
    float linear_squared, angular_squared;
    {
        // The exact installed float vectors must pass the core's own predicate,
        // not just a double norm of the source state. This is preview only.
        [[maybe_unused]] FPExceptionDisableOverflow overflow;
        [[maybe_unused]] FPExceptionDisableInvalid invalid;
        linear_squared = PhysicsVector(state.mLinearVelocity).LengthSq();
        angular_squared = PhysicsBivector(state.mAngularVelocity).LengthSq();
    }
    return std::isfinite(linear_squared) && std::isfinite(angular_squared)
        && linear_squared <= max_linear * max_linear && angular_squared <= max_angular * max_angular;
}

/// Caller exclusively owns a quiescent PhysicsSystem (no concurrent Update/readers).
/// Native BodyID stays unchanged on install; semantic ID is checked, never reassigned.
/// Mass/COM, geometry, pose, velocities, broad phase and contact invalidation are one
/// logical boundary operation. These are not rollback/transaction primitives: reserve
/// all fragment bodies and all payloads BEFORE installing any participant.
inline bool InstallPhysicsBody(PhysicsSystem &system, BodyID id, const PreparedPhysicsBody &prepared)
{
    if (prepared.GetShape() == nullptr) return false;
    {
        BodyLockRead lock(system.GetBodyLockInterface(), id);
        if (!lock.Succeeded()) return false;
        const Body &body = lock.GetBody();
        if (!body.IsDynamic() || !body.IsInBroadPhase() || body.GetUserData() != prepared.GetState().mID
            || body.GetMotionProperties()->GetAllowedDOFs() != EAllowedDOFs::All
            || !PhysicsVelocityFits(prepared.GetState(), body.GetMotionProperties()->GetMaxLinearVelocity(), body.GetMotionProperties()->GetMaxAngularVelocity())) return false;
    }
    auto &bi = system.GetBodyInterface();
    bi.SetShape(id, prepared.GetShape(), false, EActivation::Activate);
    {
        BodyLockWrite lock(system.GetBodyLockInterface(), id);
        if (!lock.Succeeded()) return false; // Exclusive ownership is required throughout.
        lock.GetBody().GetMotionProperties()->SetMassProperties(EAllowedDOFs::All, prepared.GetMass());
    }
    bi.SetPositionAndRotation(id, PhysicsPosition(prepared.GetState().mOrigin), prepared.GetState().mRotation, EActivation::Activate);
    bi.SetLinearAndAngularVelocity(id, PhysicsVector(prepared.GetState().mLinearVelocity), PhysicsBivector(prepared.GetState().mAngularVelocity));
    bi.InvalidateContactCache(id);
    return true;
}

/// A null return leaves the world unchanged (body-capacity failure is explicit).
/// Use BodyInterface::CreateBody instead when reserving a multi-body transaction;
/// this convenience operation publishes ONE newly owned body immediately.
inline BodyID CreatePhysicsBody(PhysicsSystem &system, const PreparedPhysicsBody &prepared, BodyCreationSettings settings)
{
    if (prepared.GetShape() == nullptr || settings.mMotionType != EMotionType::Dynamic
        || settings.mAllowedDOFs != EAllowedDOFs::All
        || !PhysicsVelocityFits(prepared.GetState(), settings.mMaxLinearVelocity, settings.mMaxAngularVelocity)) return BodyID();
    settings.SetShape(prepared.GetShape());
    settings.mPosition = PhysicsPosition(prepared.GetState().mOrigin);
    settings.mRotation = prepared.GetState().mRotation;
    settings.mLinearVelocity = PhysicsVector(prepared.GetState().mLinearVelocity);
    settings.mAngularVelocity = PhysicsBivector(prepared.GetState().mAngularVelocity);
    settings.mUserData = prepared.GetState().mID;
    settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
    settings.mMassPropertiesOverride = prepared.GetMass();
    return system.GetBodyInterface().CreateAndAddBody(settings, EActivation::Activate);
}

/// Components are immutable model state supplied by caller. Physics owns pose and
/// velocities: recapture before every mass, flux or topology preparation.
inline bool CapturePhysicsBody(PhysicsSystem &system, BodyID id, MassBody &ioState)
{
    BodyLockRead lock(system.GetBodyLockInterface(), id);
    if (!lock.Succeeded() || lock.GetBody().GetUserData() != ioState.mID || !lock.GetBody().IsDynamic()) return false;
    const Body &body = lock.GetBody();
    MassBody state = ioState;
    auto position = body.GetPosition();
    for (uint i = 0; i < 4; ++i) { state.mOrigin[i] = double(position[i]); state.mLinearVelocity[i] = body.GetLinearVelocity()[i]; }
    for (uint i = 0; i < 6; ++i) state.mAngularVelocity[i] = body.GetAngularVelocity()[i];
    state.mRotation = body.GetRotation();
    BodyMassData check;
    if (!MassModel::Derive(state, check)) return false;
    ioState = std::move(state);
    return true;
}

struct PhysicsBinding { uint64 mID = 0; BodyID mBody; };

/// Snapshot conventions: EntityState.position is the GEOMETRY origin, velocity is
/// COM velocity. Principal plane inertia and its actual local rotor come from the
/// installed core, not default unit tensors or a fictitious tensor-as-rotor encoding.
inline bool CaptureEntityPhysicsBody(PhysicsSystem &system, BodyID id, EntityState &ioEntity)
{
    BodyLockRead lock(system.GetBodyLockInterface(), id);
    if (!lock.Succeeded() || !lock.GetBody().IsDynamic() || lock.GetBody().GetUserData() != ioEntity.id) return false;
    const Body &body = lock.GetBody();
    const MotionProperties &motion = *body.GetMotionProperties();
    EntityState entity = ioEntity;
    auto position = body.GetPosition();
    entity.mass = 1.0f / motion.GetInverseMass();
    for (uint i = 0; i < 4; ++i)
    {
        entity.position[i] = double(position[i]);
        entity.velocity[i] = body.GetLinearVelocity()[i];
        entity.center_of_mass[i] = body.GetShape()->GetCenterOfMass()[i];
    }
    for (uint i = 0; i < 6; ++i)
    {
        float inverse = motion.GetInverseInertiaDiagonal()[i];
        if (!(inverse > 0)) return false;
        entity.inertia[i] = 1.0f / inverse;
        entity.angular_velocity[i] = body.GetAngularVelocity()[i];
    }
    for (uint i = 0; i < 8; ++i)
    {
        entity.orientation[i] = body.GetRotation().GetVec8()[i];
        entity.inertia_orientation[i] = motion.GetInertiaRotation().GetVec8()[i];
    }
    if (!Finite(entity.position) || !Finite(entity.velocity) || !Finite(entity.inertia)
        || !Finite(entity.angular_velocity) || !ValidRotor(entity.orientation) || !ValidRotor(entity.inertia_orientation)) return false;
    ioEntity = entity;
    return true;
}

/// World positions and velocities both refer to COM. A fresh physical attitude is
/// committed at the SAME boundary epoch, never analytically advanced by World.
inline bool CaptureWorldPhysicsBody(PhysicsSystem &system, BodyID id, Tick epoch, ObjectState &ioObject)
{
    BodyLockRead lock(system.GetBodyLockInterface(), id);
    if (!lock.Succeeded() || !lock.GetBody().IsDynamic() || lock.GetBody().GetUserData() != ioObject.mID) return false;
    const Body &body = lock.GetBody();
    ObjectState object = ioObject;
    object.mEpoch = epoch;
    auto position = body.GetCenterOfMassPosition();
    for (uint i = 0; i < 4; ++i) { object.mPosition[i] = double(position[i]); object.mVelocity[i] = body.GetLinearVelocity()[i]; }
    object.mOrientation = body.GetRotation();
    object.mAngularVelocity = body.GetAngularVelocity();
    if (!object.IsValid()) return false;
    ioObject = object;
    return true;
}

/// Whole owned solve-world migration, exclusively at a completed physics boundary.
/// ALL live bodies must be registered dynamic bodies and present exactly once.
/// Static/kinematic support is intentionally rejected, never left in the old frame.
/// Stage World + all numeric conversions before any live mutation; rebuild broad
/// phase and invalidate every participant's contacts, then publish World revision.
/// External weapon/history caches must key that revision; this does not migrate
/// pending event payloads or a Coordinator checkpoint on the caller's behalf.
inline bool TransferPhysicsFrame(PhysicsSystem &system, World &ioWorld, const std::vector<PhysicsBinding> &bindings,
    FrameID destination, Tick epoch)
{
    BodyIDVector all;
    system.GetBodies(all);
    if (bindings.empty() || bindings.size() != all.size()) return false;
    World staged = ioWorld;
    std::vector<ObjectID> ids;
    std::vector<Vec4> shape_com;
    for (size_t i = 0; i < bindings.size(); ++i)
    {
        const auto &binding = bindings[i];
        if (std::find(all.begin(), all.end(), binding.mBody) == all.end()) return false;
        for (size_t j = 0; j < i; ++j)
            if (binding.mID == bindings[j].mID || binding.mBody == bindings[j].mBody) return false;
        auto old = staged.FindObject(binding.mID);
        if (old == nullptr || old->mTier != Tier::Detailed || old->mEpoch > epoch) return false;
        ObjectState state = *old;
        if (!CaptureWorldPhysicsBody(system, binding.mBody, epoch, state) || !staged.CommitDetailedState(state)) return false;
        ids.push_back(binding.mID);
        shape_com.push_back(system.GetBodyInterface().GetShape(binding.mBody)->GetCenterOfMass());
    }
    if (!staged.TransformInteractingSet(ids, destination, epoch)) return false;
    std::vector<RVec4> origins;
    std::vector<Vec4> velocities;
    for (size_t i = 0; i < bindings.size(); ++i)
    {
        const auto &next = *staged.FindObject(bindings[i].mID);
        Vec4 offset = next.mOrientation * shape_com[i];
        MassVector origin;
        for (uint k = 0; k < 4; ++k)
        {
            origin[k] = next.mPosition[k] - double(offset[k]);
            if (!PhysicsScalarRepresentable<Real>(origin[k], cPhysicsNarrowingBudget) || std::abs(origin[k]) > cMaxLocalRange
                || !PhysicsScalarRepresentable<float>(next.mVelocity[k], cPhysicsNarrowingBudget)) return false;
        }
        auto position = PhysicsPosition(origin);
        auto velocity = PhysicsVector(next.mVelocity);
        BodyLockRead lock(system.GetBodyLockInterface(), bindings[i].mBody);
        if (!lock.Succeeded() || !lock.GetBody().IsInBroadPhase()) return false;
        MassBody check;
        check.mLinearVelocity = next.mVelocity;
        for (uint k = 0; k < 6; ++k) check.mAngularVelocity[k] = next.mAngularVelocity[k];
        if (!PhysicsVelocityFits(check, lock.GetBody().GetMotionProperties()->GetMaxLinearVelocity(), lock.GetBody().GetMotionProperties()->GetMaxAngularVelocity())) return false;
        origins.push_back(position); velocities.push_back(velocity);
    }
    for (size_t i = 0; i < bindings.size(); ++i)
    {
        const auto &next = *staged.FindObject(bindings[i].mID);
        auto &bi = system.GetBodyInterface();
        bi.SetPositionAndRotation(bindings[i].mBody, origins[i], next.mOrientation, EActivation::Activate);
        bi.SetLinearAndAngularVelocity(bindings[i].mBody, velocities[i], next.mAngularVelocity);
        bi.InvalidateContactCache(bindings[i].mBody);
    }
    ioWorld = std::move(staged);
    return true;
}

struct ReservoirBinding { uint64 mThrusterID = 0, mFuelComponentID = 0; };
struct PreparedReservoirBurn
{
    MassTransaction mMass;
    double mWasteHeat = 0; ///< Includes the propulsion-declared heat, not in addition to it.
};

/// Explicit electric, lumped-point-reservoir rocket model. Each used thruster has
/// one unique, zero-central-inertia fuel component at its mount. All consumption
/// occurs in ONE instantaneous pulse at the sampled boundary, followed by ordinary
/// PhysicsSystem drift/contacts. This is a declared pulse approximation, not a
/// finite-duration variable-mass integrator. No extra thrust impulse may be applied.
/// Exhaust absolute velocity = incoming rigid point velocity - direction * J/dm.
/// Expel supplies changing COM/inertia, spin, exported momentum and work ledgers.
/// Electrical resource debit must pay for kinetic work AND declared heat (no hidden
/// chemical energy). Reject float component-mass quantization above 1e-12 relative
/// flux: silently changing the fuel debit would break ownership/resource closure.
/// Zero-use ticks need no mass transaction; mixed massless actuators are rejected.
inline bool PreparePointReservoirBurn(const MassBody &state, const PropulsionResult &delivery,
    const std::vector<ReservoirBinding> &bindings, PreparedReservoirBurn &out)
{
    if (!std::isfinite(delivery.mPropellantConsumed) || delivery.mPropellantConsumed <= 0
        || !PropulsionDetail::Nonnegative(delivery.mEnergyConsumed) || !PropulsionDetail::Nonnegative(delivery.mHeatGenerated)
        || bindings.size() > 4096) return false;
    BodyMassData data;
    MassMatrix rotation;
    if (!MassModel::Derive(state, data) || !MassDetail::Rotation(state.mRotation, rotation)) return false;
    for (size_t i = 0; i < bindings.size(); ++i)
    {
        if (bindings[i].mThrusterID == 0 || bindings[i].mFuelComponentID == 0) return false;
        for (size_t j = 0; j < i; ++j)
            if (bindings[i].mThrusterID == bindings[j].mThrusterID || bindings[i].mFuelComponentID == bindings[j].mFuelComponentID) return false;
    }
    auto retained = state.mComponents;
    std::vector<LumpedExhaust> exhaust;
    double total = 0;
    for (const auto &actuator : delivery.mActuators)
    {
        if (!PropulsionDetail::Valid(actuator.mThruster) || !PropulsionDetail::Nonnegative(actuator.mDeliveredImpulse)) return false;
        if (actuator.mDeliveredImpulse == 0) continue;
        const auto &engine = actuator.mThruster;
        if (!(engine.mMaxForce > 0) || !(engine.mPropellantRate > 0)) return false;
        auto binding = std::find_if(bindings.begin(), bindings.end(), [&](const ReservoirBinding &b) { return b.mThrusterID == engine.mID; });
        if (binding == bindings.end()) return false;
        auto fuel = std::find_if(retained.begin(), retained.end(), [&](const ComponentMass &c) { return c.mID == binding->mFuelComponentID; });
        if (fuel == retained.end()) return false;
        for (uint i = 0; i < 4; ++i)
        {
            if (fuel->mLocalCOM[i] != double(engine.mMount[i])) return false;
            for (uint j = 0; j < 4; ++j) if (fuel->mProperties.mInertia(i, j) != 0) return false;
        }
        double consumed;
        {
            // Only speculative scalar charge arithmetic is guarded, never a live
            // physics update. Volatile stores keep MSVC /fp:precise from sinking
            // an operation past the guard's destructor. Validate each stage before
            // using it, including underflow (which must not create free impulse).
            [[maybe_unused]] FPExceptionDisableOverflow overflow;
            volatile double duration = actuator.mDeliveredImpulse / engine.mMaxForce;
            if (!std::isnormal(duration) || !(duration > 0)) return false;
            volatile double charge = duration * engine.mPropellantRate;
            if (!std::isnormal(charge) || !(charge > 0)) return false;
            consumed = charge;
        }
        double old_mass = fuel->mProperties.mMass;
        if (consumed > old_mass) return false;
        float remaining = float(old_mass - consumed);
        double removed = old_mass - double(remaining);
        if (!(removed > 0) || std::abs(removed - consumed) > 1e-12 * consumed) return false;
        LumpedExhaust gas;
        gas.mID = engine.mID;
        gas.mMass = removed;
        gas.mWorldPosition = MassDetail::Transform(rotation, fuel->mLocalCOM);
        for (uint i = 0; i < 4; ++i) gas.mWorldPosition[i] += state.mOrigin[i];
        if (!MassModel::PointVelocity(state, gas.mWorldPosition, gas.mWorldVelocity)) return false;
        MassVector direction {engine.mDirection[0], engine.mDirection[1], engine.mDirection[2], engine.mDirection[3]};
        direction = MassDetail::Transform(rotation, direction);
        {
            [[maybe_unused]] FPExceptionDisableOverflow overflow;
            volatile double speed = actuator.mDeliveredImpulse / removed;
            if (!std::isnormal(speed) || !(speed > 0)) return false;
            for (uint i = 0; i < 4; ++i)
            {
                volatile double kick = direction[i] * speed;
                if (!std::isfinite(kick)) return false;
                volatile double velocity = gas.mWorldVelocity[i] - kick;
                if (!std::isfinite(velocity)) return false;
                gas.mWorldVelocity[i] = velocity;
            }
        }
        // Expel enforces this same envelope before its momentum / energy products.
        if (!MassDetail::Bounded(gas.mWorldVelocity)) return false;
        exhaust.push_back(gas);
        total += removed;
        if (remaining == 0) retained.erase(fuel);
        else fuel->mProperties.mMass = remaining;
    }
    if (std::abs(total - delivery.mPropellantConsumed) > 1e-12 * delivery.mPropellantConsumed) return false;
    PreparedReservoirBurn result;
    if (!MassModel::Expel(state, retained, exhaust, result.mMass)) return false;
    {
        [[maybe_unused]] FPExceptionDisableOverflow overflow;
        volatile double available = delivery.mEnergyConsumed + result.mMass.mLedger.mDissipatedEnergy;
        if (!std::isfinite(available)) return false;
        volatile double waste = available - result.mMass.mLedger.mAddedEnergy;
        if (!std::isfinite(waste) || waste < delivery.mHeatGenerated) return false;
        result.mWasteHeat = waste;
    }
    out = std::move(result);
    return true;
}
}
