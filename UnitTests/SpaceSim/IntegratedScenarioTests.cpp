// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include <Jolt/SpaceSim/Simulation.h>
#include <Jolt/SpaceSim/Propulsion.h>
#include <Jolt/SpaceSim/Authority.h>
#include <Jolt/SpaceSim/World.h>
#include <Jolt/SpaceSim/PhysicsBridge.h>
#include <Jolt/Physics/Collision/RayCast.h>
#include <Jolt/Physics/Collision/CastResult.h>
#include <type_traits>
using namespace JPH::SpaceSim;

namespace
{
constexpr float cStep = 1.0f / 128.0f;
ComponentMass BoxMass(uint64 id, Vec4Arg half, float mass, MassVector center = {})
{
    ComponentMass c;
    c.mID = id;
    c.mLocalCOM = center;
    c.mProperties.SetMassAndInertiaOfSolidBox(2 * half, 1);
    c.mProperties.ScaleToMass(mass);
    return c;
}
MassBody Ship()
{
    MassBody ship;
    ship.mID = 10;
    ship.mComponents = { BoxMass(101, Vec4::sReplicate(0.5f), 8),
        BoxMass(102, Vec4::sReplicate(0.25f), 2, {0, 2, 0, 0}) };
    ComponentMass fuel;
    fuel.mID = 103;
    fuel.mLocalCOM = {0, 2, 0, 0};
    fuel.mProperties.mMass = 2;
    fuel.mProperties.mInertia = Mat44::sZero(); // Explicit point reservoir at nozzle, not hull geometry.
    ship.mComponents.push_back(fuel);
    return ship;
}
std::vector<CombatBox> Hull()
{
    CombatBox a, b;
    a.mComponentID = 101;
    b.mComponentID = 102;
    b.mLocalCenter = Vec4(0, 2, 0, 0);
    b.mHalfExtent = Vec4::sReplicate(0.25f);
    return {a, b};
}
BodyCreationSettings Vacuum()
{
    BodyCreationSettings settings;
    settings.mObjectLayer = Layers::MOVING;
    settings.mLinearDamping = settings.mAngularDamping = settings.mGravityFactor = 0;
    settings.mAllowSleeping = false;
    settings.mApplyGyroscopicForce = true;
    settings.mMaxLinearVelocity = 1000;
    settings.mMaxAngularVelocity = 1000;
    settings.mFriction = settings.mRestitution = 0;
    return settings;
}
bool Prepare(const MassBody &body, const std::vector<CombatBox> &boxes, PreparedPhysicsBody &out)
{
    return PreparePhysicsBody(body, boxes, out);
}

template <class T, class = void> struct PublicPreparedState : std::false_type {};
template <class T> struct PublicPreparedState<T, std::void_t<decltype(std::declval<T &>().mState)>> : std::true_type {};
template <class T, class = void> struct PublicPreparedMass : std::false_type {};
template <class T> struct PublicPreparedMass<T, std::void_t<decltype(std::declval<T &>().mMass)>> : std::true_type {};
template <class T, class = void> struct PublicPreparedShape : std::false_type {};
template <class T> struct PublicPreparedShape<T, std::void_t<decltype(std::declval<T &>().mShape)>> : std::true_type {};
template <class T> void ProbeMutablePreparedPayload(T prepared, PhysicsSystem &system)
{
    if constexpr (PublicPreparedState<T>::value)
    {
        prepared.mState.mOrigin[0] = 100000000.25;
        CHECK(CreatePhysicsBody(system, prepared, Vacuum()).IsInvalid());
    }
}

// This fixture is an actual single-owner application adapter. Every callback either
// prepares/installs live core state, runs a module, or serializes a real boundary.
// No callback replaces PhysicsSystem integration with hand-authored motion.
struct Engagement
{
    PhysicsTestContext physics {cStep};
    PhysicsBridge bridge {*physics.GetSystem(), *physics.GetTempAllocator(), *physics.GetJobSystem(), cStep};
    World world {cStep};
    Limits limits;
    std::vector<PhysicsBinding> bindings;
    std::map<uint64, MassBody> models;
    std::map<uint64, std::vector<CombatBox>> geometry;
    Thruster engine;
    PropulsionResources resources {2, 4096, 1024};
    PropulsionResult last_delivery;
    double heat = 0, exported_fuel = 0, propulsion_work = 0, propulsion_debit = 0;
    double impact_heat = 0, last_impact_work = 0, weapon_debit = 0, max_impulse_residual = 0;
    std::map<uint64, double> damage_heat;
    float throttle = 0;
    bool engine_failed = false, detached = false, migrated = false;
    uint32 updates = 0, responses = 0, observations = 0, publications = 0, pair_candidates = 0;
    uint32 prepared_hits = 0;
    MassTransaction split;
    PreparedPhysicsBody retained, debris;
    BodyID reserved_debris;
    std::vector<std::pair<uint64, Vec4>> local_witnesses;
    std::vector<ImpactEvent> observed;

    PhysicsSystem &System() { return *physics.GetSystem(); }
    BodyID Native(uint64 id) const
    {
        for (auto binding : bindings) if (binding.mID == id) return binding.mBody;
        return BodyID();
    }
    bool Spawn(const MassBody &state, const std::vector<CombatBox> &boxes, Checkpoint &checkpoint)
    {
        PreparedPhysicsBody prepared;
        if (!Prepare(state, boxes, prepared)) return false;
        BodyID id = CreatePhysicsBody(System(), prepared, Vacuum());
        if (id.IsInvalid()) return false;
        bindings.push_back({state.mID, id});
        models[state.mID] = state;
        geometry[state.mID] = boxes;
        ObjectState object;
        object.mID = state.mID; object.mFrame = checkpoint.frame.id;
        object.mTier = Tier::Detailed; object.mRotationRadius = 4;
        if (!CaptureWorldPhysicsBody(System(), id, checkpoint.tick, object) || !world.AddObject(object)) return false;
        EntityState entity;
        entity.id = state.mID; entity.owner = state.mID == 10? 77 : 0;
        entity.frame = checkpoint.frame.id;
        if (!CaptureEntityPhysicsBody(System(), id, entity)) return false;
        checkpoint.entities.push_back(entity);
        return true;
    }
    bool Initialize(Checkpoint &checkpoint)
    {
        physics.ZeroGravity();
        checkpoint.physics_hz = 128;
        checkpoint.snapshot_period = 4;
        Frame frame; frame.mID = 1; frame.mVelocity = {100, 0, 0, 0};
        if (!world.AddFrame(frame)) return false;
        frame.mID = 2; frame.mVelocity[0] = 102;
        if (!frame.mOrigin.Translate({10, 0, 0, 0}, frame.mOrigin) || !world.AddFrame(frame)) return false;
        checkpoint.frame.velocity[0] = 100;
        MassBody ship = Ship();
        ship.mLinearVelocity = {4, 0, 0, 0};
        ComponentMass ammo;
        ammo.mID = 104; ammo.mProperties.mMass = 1.0f / 128;
        ammo.mProperties.mInertia = Mat44::sZero();
        ship.mComponents.push_back(ammo);
        if (!Spawn(ship, Hull(), checkpoint)) return false;
        MassBody asteroid;
        asteroid.mID = 20;
        asteroid.mOrigin = {0.15, 0, 0, 2};
        asteroid.mLinearVelocity = {-0.5, 0, 0, -0.25};
        asteroid.mComponents = {BoxMass(201, Vec4(0.8f, 0.4f, 0.4f, 0.5f), 60),
            BoxMass(202, Vec4::sReplicate(0.2f), 4, {0.9, 0, 0, 0})};
        CombatBox a, b;
        a.mComponentID = 201; a.mHalfExtent = Vec4(0.8f, 0.4f, 0.4f, 0.5f);
        b.mComponentID = 202; b.mHalfExtent = Vec4::sReplicate(0.2f); b.mLocalCenter = Vec4(0.9f, 0, 0, 0);
        if (!Spawn(asteroid, {a, b}, checkpoint)) return false;
        engine.mID = 7; engine.mMount = Vec4(0, 2, 0, 0);
        engine.mMaxForce = 16; engine.mPropellantRate = 1;
        engine.mPowerDraw = 1024; engine.mHeatRate = 8;
        return Sync(checkpoint, 0);
    }
    bool Sync(Checkpoint &checkpoint, Tick epoch)
    {
        for (auto binding : bindings)
        {
            if (!CapturePhysicsBody(System(), binding.mBody, models[binding.mID])) return false;
            auto e = std::find_if(checkpoint.entities.begin(), checkpoint.entities.end(), [&](const EntityState &v) { return v.id == binding.mID; });
            if (e == checkpoint.entities.end()) return false;
            if (!CaptureEntityPhysicsBody(System(), binding.mBody, *e)) return false;
            e->frame = checkpoint.frame.id;
            auto object = world.FindObject(binding.mID);
            if (!object) return false;
            ObjectState state = *object;
            if (!CaptureWorldPhysicsBody(System(), binding.mBody, epoch, state) || !world.CommitDetailedState(state)) return false;
            e->heat = float(damage_heat[binding.mID]);
            if (binding.mID == 10)
            {
                e->fuel = float(resources.mPropellant); e->energy = float(resources.mEnergy); e->heat = float(heat + damage_heat[binding.mID]);
                e->health = engine_failed? 0.5f : 1.0f;
                e->actuators[0] = float(engine.mForce);
                e->topology_version = detached? 1 : 0;
                e->geometry = detached? 2 : 1;
            }
        }
        return true;
    }
    CombatCompound Target(uint64 id, uint64 frame, double epoch)
    {
        auto &model = models[id];
        BodyMassData data;
        CHECK(MassModel::Derive(model, data));
        CombatCompound target;
        target.mBodyID = id; target.mBoxes = geometry[id].data(); target.mBoxCount = uint32(geometry[id].size());
        target.mMotion.mFrameID = frame; target.mMotion.mEpoch = epoch;
        // Motion is COM; QueueHit translates the fixed-origin hull assets to COM.
        target.mMotion.mPosition = PhysicsVector(data.mWorldCOM);
        target.mMotion.mRotation = model.mRotation;
        target.mMotion.mVelocity = PhysicsVector(model.mLinearVelocity);
        target.mMotion.mAngularVelocity = PhysicsBivector(model.mAngularVelocity);
        return target;
    }
    bool QueueHit(Checkpoint &draft, const TickContext &context, uint64 target_id,
        const CombatProjectile &projectile, double mass, Vec4Arg outgoing, uint64 event_id)
    {
        auto target = Target(target_id, draft.frame.id, double(context.tick - 1) * cStep);
        BodyMassData data;
        if (!MassModel::Derive(models[target_id], data)) return false;
        // Query and live geometry agree: translate fixed-origin hull coordinates to COM.
        auto boxes = geometry[target_id];
        for (auto &box : boxes) box.mLocalCenter -= PhysicsVector(data.mLocal.mLocalCOM);
        target.mBoxes = boxes.data();
        CombatQuery query; query.mFrameID = draft.frame.id; query.mEpoch = target.mMotion.mEpoch; query.mDuration = cStep;
        auto hit = SweepProjectile(projectile, &target, 1, query);
        if (hit.mStatus != CombatQueryStatus::Hit) return false;
        auto rotation_at_hit = Rotor::sExp(target.mMotion.mAngularVelocity * float(0.5 * hit.mTime)) * target.mMotion.mRotation;
        Vec4 offset = hit.mPosition - (target.mMotion.mPosition + target.mMotion.mVelocity * float(hit.mTime));
        local_witnesses.push_back({event_id, rotation_at_hit.Reversed() * offset});
        auto impact = PrepareCombatImpact(CombatResponseOwner::QueryEvent, mass, projectile.mMotion.mVelocity, outgoing, offset);
        if (!impact.mValid || !impact.mApplyImpulse || impact.mProjectileEnergyLoss < 0) return false;
        ImpactEvent event;
        event.tick = context.tick; event.id = event_id; event.source = 90;
        event.target = target_id; event.component = hit.mComponentID; event.frame = draft.frame.id;
        event.fraction = uint64(std::min(1000000.0, hit.mTime / cStep * 1000000));
        event.owner = ResponseOwner::Query; event.energy = float(impact.mProjectileEnergyLoss);
        for (uint i = 0; i < 4; ++i) event.impulse[i] = impact.mTargetImpulse[i];
        draft.events.push_back(event);
        ++prepared_hits;
        return true;
    }
    bool Reserve(const Checkpoint &, const TickContext &context)
    {
        if (context.tick != 12) return System().GetNumBodies() == bindings.size();
        auto &ship = models[10];
        if (!CapturePhysicsBody(System(), Native(10), ship)
            || !MassModel::Fragment(ship, {{10, {101, 103}}, {11, {102}}}, split) || split.mBodies.size() != 2
            || !Prepare(split.mBodies[0], {Hull()[0]}, retained) || !Prepare(split.mBodies[1], {Hull()[1]}, debris)) return false;
        // Allocate detached body BEFORE any boundary writes, but do not add it yet.
        auto settings = Vacuum(); settings.SetShape(debris.GetShape());
        settings.mPosition = PhysicsPosition(debris.GetState().mOrigin); settings.mRotation = debris.GetState().mRotation;
        settings.mLinearVelocity = PhysicsVector(debris.GetState().mLinearVelocity); settings.mAngularVelocity = PhysicsBivector(debris.GetState().mAngularVelocity);
        settings.mUserData = 11; settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
        settings.mMassPropertiesOverride = debris.GetMass();
        Body *body = System().GetBodyInterface().CreateBody(settings);
        if (!body) return false;
        reserved_debris = body->GetID();
        return true;
    }
    bool StageUpdate(Stage stage, Checkpoint &draft, const TickContext &context)
    {
        if (stage == Stage::Boundary)
        {
            if (context.tick == 12)
            {
                if (reserved_debris.IsInvalid() || !InstallPhysicsBody(System(), Native(10), retained)) return false;
                System().GetBodyInterface().AddBody(reserved_debris, EActivation::Activate);
                bindings.push_back({11, reserved_debris}); reserved_debris = BodyID();
                models[10] = split.mBodies[0]; models[11] = split.mBodies[1];
                geometry[10] = {Hull()[0]}; geometry[11] = {Hull()[1]};
                ObjectState object; object.mID = 11; object.mFrame = draft.frame.id; object.mTier = Tier::Detailed; object.mRotationRadius = 3;
                if (!CaptureWorldPhysicsBody(System(), Native(11), context.tick - 1, object) || !world.AddObject(object)) return false;
                EntityState entity; entity.id = 11; entity.frame = draft.frame.id;
                if (!CaptureEntityPhysicsBody(System(), Native(11), entity)) return false;
                draft.entities.push_back(entity);
                detached = true;
            }
            if (context.tick == 16)
            {
                if (!draft.events.empty() || !TransferPhysicsFrame(System(), world, bindings, 2, context.tick - 1)) return false;
                draft.frame.id = 2; draft.frame.epoch = 0; draft.frame.origin = {10, 0, 0, 0}; draft.frame.velocity = {102, 0, 0, 0};
                for (auto &entity : draft.entities) entity.frame = 2;
                migrated = true;
            }
            return Sync(draft, context.tick - 1);
        }
        if (stage == Stage::Region)
        {
            ActivationResult result;
            if (!world.FindInteractions(context.tick - 1, 1, result)) return false;
            pair_candidates += uint32(result.mPairs.size());
        }
        if (stage == Stage::Resources)
        {
            for (const auto &command : context.commands) if (command.entity == 10) throttle = command.translation[0];
            auto &ship = models[10]; BodyMassData data;
            if (!CapturePhysicsBody(System(), Native(10), ship) || !MassModel::Derive(ship, data)) return false;
            PropulsionRequest request; request.mForce = Vec4(16 * throttle, 0, 0, 0);
            for (uint i = 4; i < 10; ++i) request.mWeights[i] = 0;
            if (!StepPropulsion(request, {engine}, PhysicsVector(data.mLocal.mLocalCOM), resources, cStep, last_delivery)) return false;
            if (last_delivery.mPropellantConsumed > 0)
            {
                PreparedReservoirBurn burn; PreparedPhysicsBody prepared;
                if (!PreparePointReservoirBurn(ship, last_delivery, {{7, 103}}, burn)
                    || !Prepare(burn.mMass.mBodies[0], geometry[10], prepared) || !InstallPhysicsBody(System(), Native(10), prepared)) return false;
                ship = burn.mMass.mBodies[0];
                exported_fuel += burn.mMass.mLedger.mExported.mMass;
                propulsion_work += burn.mMass.mLedger.mAddedEnergy;
                heat += burn.mWasteHeat;
            }
            propulsion_debit += last_delivery.mEnergyConsumed;
            resources = last_delivery.mResources;
            engine = last_delivery.mActuators[0].mThruster;
            return Sync(draft, context.tick - 1);
        }
        if (stage == Stage::Forces)
        {
            local_witnesses.clear();
            if (context.tick == 4)
            {
                // Fire a carried point projectile by exporting its actual mass and
                // momentum; Expel supplies recoil once. No extra launch impulse.
                auto &ship = models[10];
                auto retained_components = ship.mComponents;
                auto ammo = std::find_if(retained_components.begin(), retained_components.end(), [](const ComponentMass &c) { return c.mID == 104; });
                if (ammo == retained_components.end()) return false;
                LumpedExhaust shot; shot.mID = 104; shot.mMass = ammo->mProperties.mMass; shot.mWorldPosition = ship.mOrigin;
                if (!MassModel::PointVelocity(ship, shot.mWorldPosition, shot.mWorldVelocity)) return false;
                auto delta = ship.mRotation * Vec4(0, 0, 0, 512);
                for (uint i = 0; i < 4; ++i) shot.mWorldVelocity[i] += delta[i];
                retained_components.erase(ammo);
                MassTransaction launch; PreparedPhysicsBody prepared;
                if (!MassModel::Expel(ship, retained_components, {shot}, launch)
                    || !Prepare(launch.mBodies[0], geometry[10], prepared)) return false;
                // Instantaneous weapon discharge draws from the actual onboard
                // electrical store; authorize BEFORE installing projectile recoil.
                if (launch.mLedger.mAddedEnergy > resources.mEnergy
                    || !InstallPhysicsBody(System(), Native(10), prepared)) return false;
                resources.mEnergy -= launch.mLedger.mAddedEnergy;
                weapon_debit += launch.mLedger.mAddedEnergy;
                ship = launch.mBodies[0];
                CombatProjectile projectile; projectile.mRadius = 0.02;
                projectile.mMotion.mFrameID = draft.frame.id; projectile.mMotion.mEpoch = double(context.tick - 1) * cStep;
                projectile.mMotion.mPosition = PhysicsVector(shot.mWorldPosition); projectile.mMotion.mVelocity = PhysicsVector(shot.mWorldVelocity);
                if (!QueueHit(draft, context, 20, projectile, shot.mMass, projectile.mMotion.mVelocity - 0.5f * delta, 1001)) return false;
            }
            if (context.tick == 8)
            {
                // An externally owned incoming round scatters through the engine.
                // Its incoming/outgoing momentum is explicit, not an HP impulse.
                auto target = Target(10, draft.frame.id, double(context.tick - 1) * cStep);
                BodyMassData data; if (!MassModel::Derive(models[10], data)) return false;
                Vec4 engine_offset = PhysicsVector(MassDetail::Sub(MassVector{0, 2, 0, 0}, data.mLocal.mLocalCOM));
                CombatProjectile round; round.mRadius = 0.02;
                round.mMotion.mFrameID = draft.frame.id; round.mMotion.mEpoch = target.mMotion.mEpoch;
                Vec4 delta = target.mMotion.mRotation * Vec4(0, 0, 0, 512);
                round.mMotion.mPosition = target.mMotion.mPosition + target.mMotion.mRotation * (engine_offset - Vec4(0, 0, 0, 2));
                round.mMotion.mVelocity = target.mMotion.mVelocity + delta;
                if (!QueueHit(draft, context, 10, round, 0.001, round.mMotion.mVelocity - 0.5f * delta, 1002)) return false;
            }
            return Sync(draft, context.tick - 1);
        }
        if (stage == Stage::Damage)
        {
            for (const auto &event : observed)
                if (event.target == 10 && event.component == 102 && impact_heat > 10)
                    engine.mFailed = engine_failed = true;
            observed.clear();
        }
        if (stage == Stage::Validate) return bridge.ValidateBoundary() == BridgeResult::Success && Sync(draft, context.tick);
        return true;
    }
    Callbacks Adapters()
    {
        Callbacks cb;
        cb.reserve = [&](const Checkpoint &s, const TickContext &t) { return Reserve(s, t); };
        cb.release = [&](bool) { if (!reserved_debris.IsInvalid()) { System().GetBodyInterface().DestroyBody(reserved_debris); reserved_debris = BodyID(); } };
        cb.stage = [&](Stage stage, Checkpoint &s, const TickContext &t) { return StageUpdate(stage, s, t); };
        cb.physics_update = [&](float seconds, Checkpoint &draft)
        {
            if (seconds != cStep) return uint32(1);
            ++updates;
            auto error = bridge.Step(nullptr, 0);
            if (error != BridgeResult::Success) return uint32(1);
            // Query response is intentionally deferred to tick-end. Advect the
            // surface witness using the ACTUAL integrated endpoint orientation.
            for (auto &event : draft.events)
                for (const auto &witness : local_witnesses) if (witness.first == event.id)
                {
                    Vec4 impulse(event.impulse[0], event.impulse[1], event.impulse[2], event.impulse[3]);
                    Bivec angular = Bivec::sWedge(impulse, System().GetBodyInterface().GetRotation(Native(event.target)) * witness.second);
                    for (uint i = 0; i < 6; ++i) event.angular_impulse[i] = angular[i];
                }
            return Sync(draft, draft.tick + 1)? uint32(0) : uint32(1);
        };
        cb.apply_query = [&](Checkpoint &draft, const ImpactEvent &event)
        {
            BodyID id = Native(event.target);
            if (id.IsInvalid()) return false;
            auto &model = models[event.target];
            BodyMassData before, after;
            if (!CapturePhysicsBody(System(), id, model) || !MassModel::Derive(model, before)) return false;
            Vec4 impulse(event.impulse[0], event.impulse[1], event.impulse[2], event.impulse[3]);
            Bivec angular(event.angular_impulse[0], event.angular_impulse[1], event.angular_impulse[2], event.angular_impulse[3], event.angular_impulse[4], event.angular_impulse[5]);
            if (bridge.ApplyBoundaryImpulse({id, impulse, angular}) != BridgeResult::Success) return false;
            if (!CapturePhysicsBody(System(), id, model) || !MassModel::Derive(model, after)) return false;
            last_impact_work = after.mMomentum.mKineticEnergy - before.mMomentum.mKineticEnergy;
            MassVector reference = before.mWorldCOM;
            // Impulses change only velocities at this boundary: spin residual is
            // measured about the common COM and not a moving global reference.
            for (uint i = 0; i < 4; ++i)
            {
                max_impulse_residual = std::max(max_impulse_residual,
                    std::abs(after.mMomentum.mLinear[i] - before.mMomentum.mLinear[i] - double(impulse[i])));
            }
            MassBivector orbital_before = MassDetail::Wedge(before.mMomentum.mLinear, reference);
            MassBivector orbital_after = MassDetail::Wedge(after.mMomentum.mLinear, reference);
            for (uint i = 0; i < 6; ++i)
                max_impulse_residual = std::max(max_impulse_residual, std::abs(
                    (after.mMomentum.mAngular[i] - orbital_after[i]) - (before.mMomentum.mAngular[i] - orbital_before[i]) - double(angular[i])));
            double deposited = event.energy - last_impact_work;
            if (deposited < -1e-4) return false;
            damage_heat[event.target] += std::max(0.0, deposited);
            if (event.target == 10) impact_heat = deposited;
            ++responses;
            return Sync(draft, draft.tick + 1);
        };
        cb.observe = [&](const ImpactEvent &event) { observed.push_back(event); ++observations; return true; };
        cb.publish = [&](const Checkpoint &state, bool) { std::vector<uint8> bytes; CHECK(EncodeCheckpoint(state, limits, bytes)); ++publications; };
        return cb;
    }
};

// Application-owned, versioned immutable assets and causal sidecar. EntityState
// is the compact network view, NOT a lossless application rollback image. Exact
// resource/actuator doubles and dynamic component membership belong to this
// sampled checkpoint; model geometry/inertia and engine specifications are assets.
struct ReplayAssets
{
    uint64 model_version = 1;
    std::vector<ComponentMass> components;
    std::map<uint64, std::vector<CombatBox>> hulls;
    Thruster engine;
};
struct ComponentState { uint64 id; float mass; };
struct ReplaySnapshot
{
    uint64 model_version = 0;
    Tick tick = 0;
    FrameState frame;
    EntityState entity;
    std::vector<ComponentState> components;
    PropulsionResources resources;
    double force = 0, actuator_heat = 0, heat = 0, damage_heat = 0;
    float controller = 0;
    bool failed = false;
};
ReplaySnapshot SampleReplay(const Engagement &scene, const Checkpoint &checkpoint)
{
    ReplaySnapshot snapshot;
    snapshot.model_version = checkpoint.model_version;
    snapshot.tick = checkpoint.tick; snapshot.frame = checkpoint.frame;
    auto entity = std::find_if(checkpoint.entities.begin(), checkpoint.entities.end(), [](const EntityState &e) { return e.id == 10; });
    if (entity == checkpoint.entities.end()) return snapshot;
    snapshot.entity = *entity;
    for (const auto &c : scene.models.at(10).mComponents) snapshot.components.push_back({c.mID, c.mProperties.mMass});
    snapshot.resources = scene.resources;
    snapshot.force = scene.engine.mForce; snapshot.actuator_heat = scene.engine.mHeat;
    snapshot.failed = scene.engine.mFailed; snapshot.controller = scene.throttle;
    snapshot.heat = scene.heat;
    auto damage = scene.damage_heat.find(10);
    if (damage != scene.damage_heat.end()) snapshot.damage_heat = damage->second;
    return snapshot;
}

// Isolated local-entity replay. Future external query responses are immutable
// authoritative INPUTS for the interval, not borrowed mutable server state. This
// restores no coupled-world/contact caches; the tested interval is contact-free.
struct LocalReplay
{
    const ReplayAssets assets;
    const ReplaySnapshot checkpoint;
    const std::vector<ImpactEvent> impacts;
    PhysicsTestContext physics {cStep};
    PhysicsBridge bridge {*physics.GetSystem(), *physics.GetTempAllocator(), *physics.GetJobSystem(), cStep};
    MassBody model;
    Thruster engine;
    PropulsionResources resources;
    BodyID id;
    Tick tick = 0;
    double heat = 0, damage_heat = 0, consumed = 0;
    float controller = 0;
    uint32 updates = 0, active_ticks = 0, failed_ticks = 0, responses = 0;
    bool valid = false;

    LocalReplay(const ReplayAssets &inAssets, const ReplaySnapshot &snapshot, const std::vector<ImpactEvent> &inputs = {})
        : assets(inAssets), checkpoint(snapshot), impacts(inputs)
    {
        physics.ZeroGravity();
        if (snapshot.model_version != assets.model_version || snapshot.entity.frame != snapshot.frame.id
            || assets.hulls.find(snapshot.entity.geometry) == assets.hulls.end()) return;
        model.mID = snapshot.entity.id; model.mOrigin = snapshot.entity.position;
        for (uint i = 0; i < 4; ++i) model.mLinearVelocity[i] = snapshot.entity.velocity[i];
        for (uint i = 0; i < 6; ++i) model.mAngularVelocity[i] = snapshot.entity.angular_velocity[i];
        const auto &r = snapshot.entity.orientation;
        model.mRotation = Rotor(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
        for (const auto &state : snapshot.components)
        {
            auto asset = std::find_if(assets.components.begin(), assets.components.end(), [&](const ComponentMass &c) { return c.mID == state.id; });
            if (asset == assets.components.end() || !std::isfinite(state.mass) || state.mass <= 0) return;
            auto component = *asset;
            component.mProperties.ScaleToMass(state.mass);
            model.mComponents.push_back(component);
        }
        engine = assets.engine;
        engine.mForce = snapshot.force; engine.mHeat = snapshot.actuator_heat; engine.mFailed = snapshot.failed;
        resources = snapshot.resources; controller = snapshot.controller;
        heat = snapshot.heat; damage_heat = snapshot.damage_heat; tick = snapshot.tick;
        PreparedPhysicsBody prepared;
        if (!Prepare(model, assets.hulls.at(snapshot.entity.geometry), prepared)) return;
        id = CreatePhysicsBody(*physics.GetSystem(), prepared, Vacuum());
        valid = !id.IsInvalid();
    }
    bool Step(EntityState &entity, Tick next_tick, float seconds, const std::vector<Command> &commands)
    {
        if (!valid || next_tick != tick + 1 || seconds != cStep || entity.frame != checkpoint.frame.id) return false;
        for (const auto &command : commands) controller = command.translation[0];
        if (!CapturePhysicsBody(*physics.GetSystem(), id, model)) return false;
        BodyMassData data;
        PropulsionRequest request; request.mForce = Vec4(16 * controller, 0, 0, 0);
        for (uint i = 4; i < 10; ++i) request.mWeights[i] = 0;
        PropulsionResult delivery;
        if (!MassModel::Derive(model, data) || !StepPropulsion(request, {engine}, PhysicsVector(data.mLocal.mLocalCOM), resources, seconds, delivery)) return false;
        if (delivery.mPropellantConsumed > 0)
        {
            PreparedReservoirBurn burn; PreparedPhysicsBody prepared;
            if (!PreparePointReservoirBurn(model, delivery, {{7, 103}}, burn)
                || !Prepare(burn.mMass.mBodies[0], assets.hulls.at(entity.geometry), prepared)
                || !InstallPhysicsBody(*physics.GetSystem(), id, prepared)) return false;
            model = burn.mMass.mBodies[0]; heat += burn.mWasteHeat;
            consumed += delivery.mPropellantConsumed; ++active_ticks;
        }
        if (engine.mFailed) ++failed_ticks;
        engine = delivery.mActuators[0].mThruster; resources = delivery.mResources;
        if (bridge.Step(nullptr, 0) != BridgeResult::Success) return false;
        ++updates;
        for (const auto &event : impacts) if (event.tick == next_tick && event.target == entity.id)
        {
            if (event.frame != entity.frame || event.owner != ResponseOwner::Query) return false;
            BodyMassData before, after;
            if (!CapturePhysicsBody(*physics.GetSystem(), id, model) || !MassModel::Derive(model, before)) return false;
            const auto &j = event.impulse; const auto &k = event.angular_impulse;
            if (bridge.ApplyBoundaryImpulse({id, Vec4(j[0], j[1], j[2], j[3]), Bivec(k[0], k[1], k[2], k[3], k[4], k[5])}) != BridgeResult::Success
                || !CapturePhysicsBody(*physics.GetSystem(), id, model) || !MassModel::Derive(model, after)) return false;
            double deposited = event.energy - (after.mMomentum.mKineticEnergy - before.mMomentum.mKineticEnergy);
            if (deposited < -1e-4) return false;
            damage_heat += std::max(0.0, deposited);
            if (event.component == 102 && deposited > 10) engine.mFailed = true;
            ++responses;
        }
        if (bridge.ValidateBoundary() != BridgeResult::Success || !CapturePhysicsBody(*physics.GetSystem(), id, model)
            || !CaptureEntityPhysicsBody(*physics.GetSystem(), id, entity)) return false;
        entity.fuel = float(resources.mPropellant); entity.energy = float(resources.mEnergy);
        entity.heat = float(heat + damage_heat); entity.health = engine.mFailed? 0.5f : 1;
        entity.actuators[0] = float(engine.mForce);
        tick = next_tick;
        return true;
    }
};
}

TEST_SUITE("IntegratedSpaceSimTests")
{
    TEST_CASE("AuthoritativeEngagementBurnsFiresLosesEngineDetachesMigratesAndReconciles")
    {
        Engagement scene;
        Checkpoint initial;
        CHECK(scene.Initialize(initial));
        Coordinator server(initial, scene.limits);
        CHECK_FALSE(server.Stopped());
        if (server.Stopped()) return;
        auto callbacks = scene.Adapters();
        PredictionBuffer prediction(77, 10, 1, scene.limits);
        Command first; first.client = 77; first.entity = 10; first.sequence = 1; first.tick = 1; first.translation[0] = 1;
        CHECK(prediction.Queue(first));
        CHECK(server.Submit(first) == CommandResult::Accepted);
        const ReplayAssets assets {initial.model_version, scene.models[10].mComponents,
            {{1, Hull()}, {2, {Hull()[0]}}}, scene.engine};
        auto original_id = scene.Native(10);
        auto advance = [&]()
        {
            auto result = server.Advance(7812500, 1, callbacks);
            CAPTURE(server.State().tick); CAPTURE(int(result.error));
            CHECK(result.error == Error::None);
            CHECK(result.steps == 1);
            return result.error == Error::None && result.steps == 1;
        };
        if (!advance()) return;
        CHECK(server.State().streams[0].applied == 1);
        Command delayed = first; delayed.sequence = 2;
        CHECK(prediction.Queue(delayed));
        CHECK(server.Submit(delayed) == CommandResult::Late);
        Command next = first; next.sequence = 3; next.tick = 2;
        CHECK(prediction.Queue(next));
        CHECK(server.Submit(next) == CommandResult::Accepted);
        Command future = next; future.sequence = 4; future.tick = 20;
        CHECK(prediction.Queue(future));
        CHECK(server.Submit(future) == CommandResult::Accepted);
        Command intruder = future; intruder.client = 88; intruder.sequence = 100;
        CHECK(server.Submit(intruder) == CommandResult::Authority);
        while (server.State().tick < 16)
        {
            if (!advance()) return;
            if (server.State().tick == 8)
            {
                CHECK(scene.engine_failed);
                CHECK(scene.impact_heat > 10);
            }
            if (server.State().tick == 9)
            {
                CHECK(scene.last_delivery.mImpulse == Vec4::sZero());
                CHECK(scene.last_delivery.mAllocationResidualForce.GetX() == doctest::Approx(16));
                CHECK(scene.last_delivery.mPropellantConsumed == 0);
            }
        }
        CHECK(scene.prepared_hits == 2);
        CHECK(scene.responses == 2); // Coordinator's real ImpactBatch delivered each query once.
        CHECK(scene.observations == 2);
        CHECK(scene.detached); CHECK(scene.migrated);
        CHECK(scene.Native(10) == original_id);
        CHECK(scene.System().GetNumBodies() == 3);
        CHECK(scene.resources.mPropellant == doctest::Approx(2 - 8.0 / 128));
        CHECK(scene.exported_fuel == doctest::Approx(8.0 / 128));
        CHECK(scene.resources.mEnergy == doctest::Approx(4096 - scene.propulsion_debit - scene.weapon_debit));
        CHECK(scene.weapon_debit > 1024);
        CHECK(scene.bridge.GetTick() == 16);
        CHECK(scene.max_impulse_residual < 1e-4);
        CHECK(scene.propulsion_work + scene.heat == doctest::Approx(scene.propulsion_debit).epsilon(1e-10));
        CHECK(scene.pair_candidates > 0);
        CHECK(server.State().frame.id == 2);
        CHECK(server.State().streams[0].accepted == 4);
        CHECK(server.State().streams[0].applied == 3);
        auto snapshot = server.Published();
        CHECK(snapshot->tick == 16);
        CHECK(server.Historical(15)->frame.id == 1);
        CHECK(server.Historical(16)->frame.id == 2);
        std::vector<uint8> wire;
        CHECK(EncodeCheckpoint(*snapshot, scene.limits, wire));
        Checkpoint received;
        CHECK(DecodeCheckpoint(wire, scene.limits, received));

        // Reconcile from a delayed authoritative image, replaying four REAL core
        // steps in an isolated local-entity world. No collisions occur in this
        // replay interval; this is not coupled-world rollback or cache restoration.
        // The only server reads are the immutable checkpoint/asset sampling here.
        // Restore and every replay tick below have NO live Engagement dependency.
        const ReplaySnapshot replay_checkpoint = SampleReplay(scene, received);
        LocalReplay client_replay(assets, replay_checkpoint);
        CHECK(client_replay.valid);
        auto replay = [&client_replay](EntityState &entity, Tick tick, float seconds, const std::vector<Command> &commands)
        {
            return client_replay.Step(entity, tick, seconds, commands);
        };
        CHECK(prediction.Reconcile(received, 20, replay) == ReconcileResult::Applied);
        CHECK(client_replay.updates == 4);
        // Applied sequence 3 also covers the rejected lower sequence 2.
        CHECK(prediction.LateDiscarded() == 0);
        CHECK(prediction.Pending().size() == 1);
        while (server.State().tick < 20) if (!advance()) return;
        const auto &authoritative = server.State().entities.front();
        for (uint i = 0; i < 4; ++i)
        {
            CHECK(std::abs(prediction.Predicted().position[i] - authoritative.position[i]) < 2e-5);
            CHECK(std::abs(prediction.Predicted().velocity[i] - authoritative.velocity[i]) < 2e-5);
        }
        CHECK(server.State().streams[0].applied == 4);
        CHECK(authoritative.mass == doctest::Approx(8 + scene.resources.mPropellant).epsilon(1e-6));
        CHECK(authoritative.geometry == 2);
        CHECK(authoritative.topology_version == 1);
        CHECK(prediction.Predicted().frame == 2);
        CHECK(prediction.Reconcile(server.State(), 20, replay) == ReconcileResult::Applied);
        CHECK(prediction.Pending().empty());
        CHECK(client_replay.updates == 4);
        Command late_packet = future; late_packet.sequence = 5; late_packet.tick = 21;
        CHECK(prediction.Queue(late_packet));
        if (!advance()) return;
        CHECK(server.Submit(late_packet) == CommandResult::Late);
        CHECK(prediction.Reconcile(server.State(), 21, replay) == ReconcileResult::Applied);
        CHECK(prediction.LateDiscarded() == 1);
        CHECK(prediction.Pending().empty());
        CHECK(scene.updates == 21);
        CHECK(scene.publications == 21);
        CHECK(scene.world.FindObject(10)->mFrame == 2);
        CHECK(scene.world.FindObject(11)->mFrame == 2);
        CHECK(scene.world.FindObject(20)->mFrame == 2);
    }

    TEST_CASE("HistoricalReplayRestoresActiveEngineAndComponentMassNotTheCurrentScene")
    {
        Engagement scene; Checkpoint initial;
        CHECK(scene.Initialize(initial));
        const ReplayAssets assets {initial.model_version, scene.models[10].mComponents,
            {{1, Hull()}, {2, {Hull()[0]}}}, scene.engine};
        Coordinator server(initial, scene.limits);
        auto cb = scene.Adapters();
        std::vector<ImpactEvent> recorded_inputs;
        auto observe = cb.observe;
        cb.observe = [&, observe](const ImpactEvent &event)
        {
            if (event.target == 10) recorded_inputs.push_back(event);
            return observe(event);
        };
        Command command; command.client = 77; command.entity = 10; command.sequence = 1;
        command.tick = 1; command.translation[0] = 1;
        CHECK(server.Submit(command) == CommandResult::Accepted);
        for (uint i = 0; i < 4; ++i) CHECK(server.Advance(7812500, 1, cb).steps == 1);
        const Checkpoint checkpoint = server.State();
        const ReplaySnapshot snapshot = SampleReplay(scene, checkpoint);
        command.sequence = 2; command.tick = 6; command.translation[0] = 0.5f;
        CHECK(server.Submit(command) == CommandResult::Accepted);
        for (uint i = 4; i < 10; ++i) CHECK(server.Advance(7812500, 1, cb).steps == 1);
        CHECK(scene.engine_failed);
        const EntityState expected = server.State().entities.front();
        const auto expected_resources = scene.resources;
        const std::vector<ImpactEvent> replay_inputs = recorded_inputs;
        CHECK(replay_inputs.size() == 1);
        // Destroy the tempting live backing data before BOTH restore and replay.
        scene.models.clear(); scene.geometry.clear(); scene.engine = Thruster();
        scene.resources = {}; scene.throttle = 0;
        LocalReplay client(assets, snapshot, replay_inputs);
        CHECK(client.valid); if (!client.valid) return;
        CHECK(client.model.mComponents.back().mProperties.mMass == snapshot.components.back().mass);
        CHECK_FALSE(client.engine.mFailed);
        CHECK(client.engine.mForce == snapshot.force);
        PredictionBuffer prediction(77, 10, assets.model_version, Limits());
        CHECK(prediction.Queue(command));
        auto replay = [&client](EntityState &entity, Tick tick, float seconds, const std::vector<Command> &commands)
        {
            return client.Step(entity, tick, seconds, commands);
        };
        CHECK(prediction.Reconcile(checkpoint, 10, replay) == ReconcileResult::Applied);
        CHECK(client.updates == 6);
        CHECK(client.active_ticks == 4);
        CHECK(client.failed_ticks == 2);
        CHECK(client.responses == 1);
        CHECK(client.consumed > 0);
        CHECK(client.resources.mPropellant == expected_resources.mPropellant);
        CHECK(client.resources.mEnergy == expected_resources.mEnergy);
        CHECK(prediction.Predicted().fuel == expected.fuel);
        CHECK(prediction.Predicted().energy == expected.energy);
        CHECK(prediction.Predicted().heat == expected.heat);
        CHECK(prediction.Predicted().health == expected.health);
        CHECK(prediction.Predicted().actuators == expected.actuators);
        CHECK(prediction.Predicted().orientation == expected.orientation);
        CHECK(prediction.Predicted().angular_velocity == expected.angular_velocity);
        CHECK(prediction.Predicted().position == expected.position);
        CHECK(prediction.Predicted().velocity == expected.velocity);
        auto stale_snapshot = snapshot; ++stale_snapshot.model_version;
        LocalReplay stale(assets, stale_snapshot, replay_inputs);
        CHECK_FALSE(stale.valid);
    }

    TEST_CASE("CoordinatorRejectsOversizedQueryImpulseBeforePublicationOrClamp")
    {
        Engagement scene;
        Checkpoint initial;
        CHECK(scene.Initialize(initial));
        Coordinator server(initial, scene.limits);
        auto cb = scene.Adapters();
        auto stage = cb.stage;
        cb.stage = [&, stage](Stage s, Checkpoint &draft, const TickContext &t)
        {
            if (!stage(s, draft, t)) return false;
            if (s == Stage::Forces)
            {
                ImpactEvent event;
                event.tick = t.tick; event.id = 9001; event.source = 90; event.target = 10;
                event.component = 101; event.owner = ResponseOwner::Query;
                event.impulse[0] = 1e7f; event.energy = 1e12f;
                draft.events.push_back(event);
            }
            return true;
        };
        auto result = server.Advance(7812500, 1, cb);
        CHECK(result.error == Error::Impact);
        CHECK(result.steps == 0);
        CHECK(server.Stopped());
        CHECK(server.State().tick == 0);
        CHECK(scene.publications == 0);
        CHECK(scene.responses == 0);
        BodyLockRead lock(scene.System().GetBodyLockInterface(), scene.Native(10));
        CHECK_FALSE(lock.GetBody().GetMotionProperties()->WasLinearVelocityClamped());
        CHECK(lock.GetBody().GetLinearVelocity().GetX() == 4);
    }

    TEST_CASE("QueryBoundaryPreflightsTheExactFloatClampPredicate")
    {
        PhysicsTestContext physics(cStep); physics.ZeroGravity();
        PreparedPhysicsBody prepared;
        CHECK(Prepare(Ship(), Hull(), prepared));
        auto id = CreatePhysicsBody(*physics.GetSystem(), prepared, Vacuum());
        CHECK_FALSE(id.IsInvalid()); if (id.IsInvalid()) return;
        const Vec4 impulse(0.5495255589485168f, -0.5596849918365479f, 0.5335712432861328f, -0.23720909655094147f);
        const float limit = 0.9778507351875305f;
        {
            BodyLockWrite lock(physics.GetSystem()->GetBodyLockInterface(), id);
            auto &motion = *lock.GetBody().GetMotionProperties();
            motion.SetInverseMass(1); motion.SetMaxLinearVelocity(limit);
        }
        double norm = 0;
        for (uint i = 0; i < 4; ++i) norm += double(impulse[i]) * impulse[i];
        CHECK(norm <= double(limit) * limit);
        CHECK(impulse.LengthSq() > limit * limit); // Double-only preflight misses this core clamp.
        auto admitted = Ship();
        for (uint i = 0; i < 4; ++i) admitted.mLinearVelocity[i] = impulse[i];
        CHECK_FALSE(PhysicsVelocityFits(admitted, limit, 1000));
        PhysicsBridge bridge(*physics.GetSystem(), *physics.GetTempAllocator(), *physics.GetJobSystem(), cStep);
        CHECK(bridge.ApplyBoundaryImpulse({id, impulse, Bivec::sZero()}) == BridgeResult::InvalidDelivery);
        CHECK(physics.GetBodyInterface().GetLinearVelocity(id) == Vec4::sZero());
        BodyLockRead lock(physics.GetSystem()->GetBodyLockInterface(), id);
        CHECK_FALSE(lock.GetBody().GetMotionProperties()->WasLinearVelocityClamped());
        CHECK_FALSE(bridge.IsFailedStopped());
    }

    TEST_CASE("QueryBoundaryRejectsFiniteButUnrepresentableCoreSpeedLimits")
    {
        PhysicsTestContext physics(cStep); physics.ZeroGravity();
        PreparedPhysicsBody prepared;
        CHECK(Prepare(Ship(), Hull(), prepared));
        auto id = CreatePhysicsBody(*physics.GetSystem(), prepared, Vacuum());
        CHECK_FALSE(id.IsInvalid()); if (id.IsInvalid()) return;
        {
            BodyLockWrite lock(physics.GetSystem()->GetBodyLockInterface(), id);
            lock.GetBody().GetMotionProperties()->SetMaxLinearVelocity(std::numeric_limits<float>::max());
        }
        PhysicsBridge bridge(*physics.GetSystem(), *physics.GetTempAllocator(), *physics.GetJobSystem(), cStep);
        CHECK(bridge.ApplyBoundaryImpulse({id, Vec4::sAxisX(), Bivec::sZero()}) == BridgeResult::NumericalFailure);
        CHECK(physics.GetBodyInterface().GetLinearVelocity(id) == Vec4::sZero());
        CHECK(bridge.IsFailedStopped());
    }

    TEST_CASE("QueryBoundaryRejectsRoundedAwayImpulseInsteadOfClaimingDelivery")
    {
        PhysicsTestContext physics(cStep); physics.ZeroGravity();
        auto state = Ship(); state.mLinearVelocity[0] = 1e6;
        PreparedPhysicsBody prepared;
        CHECK(Prepare(state, Hull(), prepared));
        auto settings = Vacuum(); settings.mMaxLinearVelocity = 1e7f;
        auto id = CreatePhysicsBody(*physics.GetSystem(), prepared, settings);
        CHECK_FALSE(id.IsInvalid()); if (id.IsInvalid()) return;
        PhysicsBridge bridge(*physics.GetSystem(), *physics.GetTempAllocator(), *physics.GetJobSystem(), cStep);
        for (float impulse : {0.01f, 1e-6f})
            CHECK(bridge.ApplyBoundaryImpulse({id, Vec4(impulse, 0, 0, 0), Bivec::sZero()}) == BridgeResult::InvalidDelivery);
        CHECK(physics.GetBodyInterface().GetLinearVelocity(id).GetX() == 1e6f);
        CHECK_FALSE(bridge.IsFailedStopped());
    }

    TEST_CASE("ImpactBatchCannotDoubleApplyLiveImpulseOrReapplySolverOwnedResponse")
    {
        PhysicsTestContext physics(cStep); physics.ZeroGravity();
        auto &system = *physics.GetSystem();
        PreparedPhysicsBody prepared;
        auto first_body = Ship(); first_body.mOrigin[0] = -0.75; first_body.mLinearVelocity[0] = 2;
        CHECK(Prepare(first_body, Hull(), prepared));
        auto id = CreatePhysicsBody(system, prepared, Vacuum());
        CHECK_FALSE(id.IsInvalid()); if (id.IsInvalid()) return;
        auto other = Ship(); other.mID = 20; other.mOrigin[0] = 0.75; other.mLinearVelocity[0] = -2;
        CHECK(Prepare(other, Hull(), prepared));
        auto other_id = CreatePhysicsBody(system, prepared, Vacuum());
        CHECK_FALSE(other_id.IsInvalid()); if (other_id.IsInvalid()) return;
        PhysicsBridge bridge(system, *physics.GetTempAllocator(), *physics.GetJobSystem(), cStep);
        for (uint i = 0; i < 20; ++i) CHECK(bridge.Step(nullptr, 0) == BridgeResult::Success);
        CHECK(system.WereBodiesInContact(id, other_id));
        Vec4 solver_velocity = system.GetBodyInterface().GetLinearVelocity(id);
        CHECK(solver_velocity.GetX() < 1);
        ImpactEvent query;
        query.tick = 1; query.id = 100; query.source = 90; query.target = 10; query.component = 101;
        query.owner = ResponseOwner::Query; query.impulse[0] = 1;
        ImpactEvent solver = query; solver.id = 101; solver.owner = ResponseOwner::Solver;
        // Net actual solver impulse over the closed two-body contact interval.
        for (uint i = 0; i < 4; ++i) solver.impulse[i] = 12 * (solver_velocity[i] - float(first_body.mLinearVelocity[i]));
        uint32 observations = 0;
        auto apply = [&](const ImpactEvent &event)
        {
            if (event.target != 10) return false;
            return bridge.ApplyBoundaryImpulse({id, Vec4(event.impulse[0], event.impulse[1], event.impulse[2], event.impulse[3]), Bivec::sZero()}) == BridgeResult::Success;
        };
        auto observe = [&](const ImpactEvent &) { ++observations; return true; };
        ImpactBatch batch({solver, query}, Limits());
        CHECK(batch.Process(apply, observe));
        auto velocity = system.GetBodyInterface().GetLinearVelocity(id);
        CHECK(velocity.GetX() == doctest::Approx(solver_velocity.GetX() + 1.0 / 12));
        CHECK(observations == 2);
        CHECK_FALSE(batch.Process(apply, observe));
        CHECK(system.GetBodyInterface().GetLinearVelocity(id) == velocity);
        solver.id = query.id;
        ImpactBatch conflicting({query, solver}, Limits());
        CHECK_FALSE(conflicting.Process(apply, observe));
        CHECK(system.GetBodyInterface().GetLinearVelocity(id) == velocity);
        // No second core response or query delivery is allowed for either event.
    }

    TEST_CASE("PreparedBodyCapacityAndPreflightFailuresLeaveLiveStateUntouched")
    {
        PhysicsTestContext physics(cStep, 1, 0, 1); physics.ZeroGravity();
        auto &system = *physics.GetSystem();
        PreparedPhysicsBody prepared, empty;
        CHECK(CreatePhysicsBody(system, empty, Vacuum()).IsInvalid());
        CHECK(Prepare(Ship(), Hull(), prepared));
        auto settings = Vacuum(); settings.mMaxLinearVelocity = settings.mMaxAngularVelocity = 10;
        auto id = CreatePhysicsBody(system, prepared, settings);
        CHECK_FALSE(id.IsInvalid()); if (id.IsInvalid()) return;
        const auto shape = system.GetBodyInterface().GetShape(id);
        CHECK(CreatePhysicsBody(system, prepared, settings).IsInvalid()); // Capacity is one.
        CHECK(system.GetNumBodies() == 1);
        auto rejected = Ship(); rejected.mLinearVelocity[0] = 11;
        CHECK(Prepare(rejected, Hull(), prepared));
        CHECK_FALSE(InstallPhysicsBody(system, id, prepared));
        rejected = Ship(); rejected.mID = 99;
        CHECK(Prepare(rejected, Hull(), prepared));
        CHECK_FALSE(InstallPhysicsBody(system, id, prepared));
        CHECK_FALSE(InstallPhysicsBody(system, id, empty));
        CHECK(system.GetBodyInterface().GetShape(id) == shape);
        CHECK(system.GetBodyInterface().GetPosition(id) == RVec4::sZero());
        CHECK(system.GetBodyInterface().GetLinearVelocity(id) == Vec4::sZero());
        BodyLockRead lock(system.GetBodyLockInterface(), id);
        CHECK(1 / lock.GetBody().GetMotionProperties()->GetInverseMass() == 12);
        CHECK(lock.GetBody().GetUserData() == 10);
    }

    TEST_CASE("QueryBoundaryUsesActualRotatedInverseInertiaAndPreflightsBothKicks")
    {
        PhysicsTestContext physics(cStep); physics.ZeroGravity();
        auto &system = *physics.GetSystem();
        PreparedPhysicsBody prepared;
        auto state = Ship(); state.mRotation = Rotor::sExp(Bivec(0.2f, -0.1f, 0.3f, 0.1f, 0.2f, -0.1f));
        CHECK(Prepare(state, Hull(), prepared));
        auto id = CreatePhysicsBody(system, prepared, Vacuum());
        CHECK_FALSE(id.IsInvalid()); if (id.IsInvalid()) return;
        Vec4 expected_linear;
        Bivec expected_angular;
        const Vec4 impulse(0.1f, -0.2f, 0.3f, -0.4f);
        const Bivec angular(0.01f, -0.02f, 0.03f, -0.04f, 0.05f, -0.06f);
        {
            BodyLockWrite lock(system.GetBodyLockInterface(), id);
            auto &motion = *lock.GetBody().GetMotionProperties();
            motion.SetInverseInertia(Bivec(0.5f, 1, 2, 3, 4, 5), Rotor::sExp(Bivec(0, 0.2f, 0, -0.1f, 0.3f, 0)));
            motion.SetMaxLinearVelocity(1); motion.SetMaxAngularVelocity(1);
            expected_linear = impulse * motion.GetInverseMass();
            expected_angular = motion.MultiplyWorldSpaceInverseInertiaByVector(state.mRotation, angular);
        }
        PhysicsBridge bridge(system, *physics.GetTempAllocator(), *physics.GetJobSystem(), cStep);
        CHECK(bridge.ApplyBoundaryImpulse({id, impulse, Bivec(100, 0, 0, 0, 0, 0)}) == BridgeResult::InvalidDelivery);
        CHECK(system.GetBodyInterface().GetLinearVelocity(id) == Vec4::sZero());
        CHECK(system.GetBodyInterface().GetAngularVelocity(id) == Bivec::sZero());
        for (float bad : {std::numeric_limits<float>::max(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity()})
            for (uint lane = 0; lane < 10; ++lane)
            {
                float j[4] = {}, k[6] = {};
                if (lane < 4) j[lane] = bad; else k[lane - 4] = bad;
                CHECK(bridge.ApplyBoundaryImpulse({id, Vec4(j[0], j[1], j[2], j[3]), Bivec(k[0], k[1], k[2], k[3], k[4], k[5])}) == BridgeResult::InvalidDelivery);
                CHECK(system.GetBodyInterface().GetLinearVelocity(id) == Vec4::sZero());
                CHECK(system.GetBodyInterface().GetAngularVelocity(id) == Bivec::sZero());
            }
        CHECK(bridge.ApplyBoundaryImpulse({id, impulse, angular}) == BridgeResult::Success);
        CHECK(system.GetBodyInterface().GetLinearVelocity(id) == expected_linear);
        CHECK(system.GetBodyInterface().GetAngularVelocity(id) == expected_angular);
        CHECK(bridge.ValidateBoundary() == BridgeResult::Success);
        CHECK(bridge.GetTick() == 0); // Boundary impulse does NOT integrate or advance time.
    }

    TEST_CASE("CoordinatorValidatesStickyBoundaryDiagnosticsBeforePublication")
    {
        Engagement scene; Checkpoint initial;
        CHECK(scene.Initialize(initial));
        Coordinator server(initial, scene.limits);
        auto cb = scene.Adapters(); auto stage = cb.stage;
        cb.stage = [&, stage](Stage s, Checkpoint &draft, const TickContext &t)
        {
            if (!stage(s, draft, t)) return false;
            if (s == Stage::Damage)
            {
                BodyLockWrite lock(scene.System().GetBodyLockInterface(), scene.Native(10));
                auto &motion = *lock.GetBody().GetMotionProperties();
                motion.SetLinearVelocityClamped(Vec4(2000, 0, 0, 0));
                motion.SetLinearVelocityClamped(Vec4::sZero()); // Finite endpoint cannot erase evidence.
            }
            return true;
        };
        auto result = server.Advance(7812500, 1, cb);
        CHECK(result.error == Error::StageFailed);
        CHECK(result.steps == 0);
        CHECK(server.State().tick == 0);
        CHECK(scene.publications == 0);
        CHECK(scene.bridge.IsFailedStopped());
        CHECK(scene.bridge.GetLastNumericalError() == BridgeNumericalError::LinearVelocityClamped);
        CHECK(scene.bridge.Step(nullptr, 0) == BridgeResult::FailedStopped);
    }

    TEST_CASE("PreparedPayloadCannotExposeMutableStateMassOrShape")
    {
        CHECK_FALSE(PublicPreparedState<PreparedPhysicsBody>::value);
        CHECK_FALSE(PublicPreparedMass<PreparedPhysicsBody>::value);
        CHECK_FALSE(PublicPreparedShape<PreparedPhysicsBody>::value);
        PhysicsTestContext physics(cStep); physics.ZeroGravity();
        PreparedPhysicsBody prepared;
        CHECK(Prepare(Ship(), Hull(), prepared));
        ProbeMutablePreparedPayload(prepared, *physics.GetSystem());
    }

    TEST_CASE("BoundaryValidationRejectsDoubleExtremesWithNativeTrapsEnabled")
    {
        for (double bad : {std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
            std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
            for (uint lane = 0; lane < 14; ++lane)
            {
                auto state = Ship();
                if (lane < 4) state.mOrigin[lane] = bad;
                else if (lane < 8) state.mLinearVelocity[lane - 4] = bad;
                else state.mAngularVelocity[lane - 8] = bad;
                CHECK_FALSE(PhysicsStateRepresentable(state));
                if (lane >= 4) CHECK_FALSE(PhysicsVelocityFits(state, 1000, 1000));
                PreparedPhysicsBody prepared;
                CHECK_FALSE(Prepare(state, Hull(), prepared));
            }
        CHECK(PhysicsStateRepresentable(Ship()));
        CHECK(PhysicsVelocityFits(Ship(), 1000, 1000));
        CHECK_FALSE(PhysicsVelocityFits(Ship(), std::numeric_limits<float>::infinity(), 1000));
        CHECK_FALSE(PhysicsVelocityFits(Ship(), std::numeric_limits<float>::max(), 1000));
        CHECK_FALSE(PhysicsVelocityFits(Ship(), 1000, std::numeric_limits<float>::max()));
    }

    TEST_CASE("GeometryPreviewRejectsExtremeLanesBeforeNormAndKeepsPreparedPayload")
    {
        PreparedPhysicsBody prepared;
        CHECK(Prepare(Ship(), Hull(), prepared));
        const auto saved_shape = prepared.GetShape();
        const auto saved_mass = prepared.GetMass().mMass;
        for (float bad : {std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
            std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()})
            for (uint lane = 0; lane < 24; ++lane)
            {
                CAPTURE(lane); CAPTURE(bad);
#if defined(JPH_USE_SSE)
                const auto control = _mm_getcsr() & ~0x3fu;
#endif
                auto state = Ship();
                auto boxes = Hull();
                if (lane < 4) boxes.back().mLocalCenter[lane] = bad;
                else if (lane < 8) boxes.back().mHalfExtent[lane - 4] = bad;
                else
                {
                    float r[8] = {1, 0, 0, 0, 0, 0, 0, 0};
                    r[(lane - 8) % 8] = bad;
                    Rotor rotor(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
                    if (lane < 16) boxes.back().mLocalRotation = rotor;
                    else state.mRotation = rotor;
                }
                CHECK_FALSE(Prepare(state, boxes, prepared));
#if defined(JPH_USE_SSE)
                CHECK((_mm_getcsr() & ~0x3fu) == control);
#endif
                CHECK(prepared.GetShape() == saved_shape);
                CHECK(prepared.GetMass().mMass == saved_mass);
                CHECK(prepared.GetState().mOrigin == MassVector{});
                CHECK(prepared.GetState().mRotation == Rotor::sIdentity());
            }
        // Preserve the original radius <= 256 budget, including its boundary.
        auto boxes = Hull();
        boxes[0].mLocalCenter = Vec4(255, 0, 0, 0);
        CHECK(Prepare(Ship(), boxes, prepared));
        boxes[0].mLocalCenter = Vec4(255.25f, 0, 0, 0);
        CHECK_FALSE(Prepare(Ship(), boxes, prepared));
        boxes = Hull();
        boxes[0].mLocalRotation = Rotor(1.000004f, 0, 0, 0, 0, 0, 0, 0);
        CHECK(Prepare(Ship(), boxes, prepared));
        for (uint lane = 0; lane < 8; ++lane)
        {
            auto state = Ship();
            float r[8] = {1, 0, 0, 0, 0, 0, 0, 0};
            r[lane] = std::numeric_limits<float>::max();
            state.mComponents.back().mLocalRotation = Rotor(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
            CHECK_FALSE(Prepare(state, Hull(), prepared));
        }
    }

    TEST_CASE("BoundaryRejectsUnrepresentableStateWithoutPublishingPartialPayload")
    {
        PreparedPhysicsBody prepared;
        CHECK(Prepare(Ship(), Hull(), prepared));
        auto saved_shape = prepared.GetShape();
        auto large = Ship();
        large.mLinearVelocity[0] = 100000000.25; // float velocity would silently lose .25 m/s.
        CHECK_FALSE(Prepare(large, Hull(), prepared));
        CHECK(prepared.GetState().mLinearVelocity[0] == 0);
        CHECK(prepared.GetShape() == saved_shape);
#ifndef JPH_DOUBLE_PRECISION
        large = Ship(); large.mOrigin[0] = 100000000.25;
        CHECK_FALSE(Prepare(large, Hull(), prepared));
        CHECK(prepared.GetShape() == saved_shape);
#endif
    }

    TEST_CASE("FueledPulseExportsMassMomentumAndWorkInsteadOfApplyingThrustTwice")
    {
        PhysicsTestContext physics(cStep);
        physics.ZeroGravity();
        auto &system = *physics.GetSystem();
        MassBody ship = Ship();
        ship.mLinearVelocity = {4, 0, 0, 0};
        Thruster engine;
        engine.mID = 7;
        engine.mMount = Vec4(0, 2, 0, 0);
        engine.mMaxForce = 16;
        engine.mPropellantRate = 1;
        engine.mPowerDraw = 1024;
        engine.mHeatRate = 8;
        PropulsionRequest request;
        request.mForce = Vec4(16, 0, 0, 0);
        for (uint i = 4; i < 10; ++i) request.mWeights[i] = 0;
        PropulsionResult delivery;
        BodyMassData initial;
        CHECK(MassModel::Derive(ship, initial));
        CHECK(StepPropulsion(request, {engine}, PhysicsVector(initial.mLocal.mLocalCOM), {2, 100, 1024}, cStep, delivery));
        CHECK(delivery.mPropellantConsumed > 0);
        CHECK(delivery.mPropellantConsumed <= double(cStep) * engine.mPropellantRate);
        double charged_fuel = 0, charged_energy = 0;
        for (const auto &actuator : delivery.mActuators)
        {
            double duration = actuator.mDeliveredImpulse / actuator.mThruster.mMaxForce;
            charged_fuel += duration * actuator.mThruster.mPropellantRate;
            charged_energy += duration * actuator.mThruster.mPowerDraw;
        }
        CHECK(delivery.mPropellantConsumed == charged_fuel);
        CHECK(delivery.mEnergyConsumed == charged_energy);
        CHECK(delivery.mEnergyConsumed <= double(cStep) * 1024);
        CHECK(delivery.mResources.mPropellant == 2 - charged_fuel);
        CHECK(delivery.mResources.mEnergy == 100 - charged_energy);
        PreparedReservoirBurn burn;
        CHECK(PreparePointReservoirBurn(ship, delivery, {{7, 103}}, burn));
        if (burn.mMass.mBodies.empty()) return;
        CHECK(burn.mMass.mExhaust.size() == 1);
        CHECK(burn.mMass.mLedger.mAddedEnergy > 0);
        double exported = burn.mMass.mLedger.mExported.mMass;
        CHECK(std::abs(exported - charged_fuel) <= 1e-12 * charged_fuel);
        CHECK(double(ship.mComponents.back().mProperties.mMass) - double(burn.mMass.mBodies.front().mComponents.back().mProperties.mMass) == exported);
        double work = burn.mMass.mLedger.mAddedEnergy - burn.mMass.mLedger.mDissipatedEnergy;
        CHECK(std::abs(burn.mWasteHeat + work - charged_energy) <= 64 * std::numeric_limits<double>::epsilon() * charged_energy);
        CHECK(burn.mWasteHeat >= delivery.mHeatGenerated);
        for (double residual : burn.mMass.mLedger.mResidual.mLinear) CHECK(std::abs(residual) < 1e-10);
        for (double residual : burn.mMass.mLedger.mResidual.mAngular) CHECK(std::abs(residual) < 1e-10);
        PreparedPhysicsBody prepared;
        CHECK(Prepare(ship, Hull(), prepared));
        BodyID id = CreatePhysicsBody(system, prepared, Vacuum());
        CHECK_FALSE(id.IsInvalid());
        if (id.IsInvalid()) return;
        CHECK(Prepare(burn.mMass.mBodies.front(), Hull(), prepared));
        CHECK(InstallPhysicsBody(system, id, prepared));
        // One pulse already changed momentum through Expel; NEVER AddForce/AddImpulse.
        CHECK(physics.SimulateSingleStep() == EPhysicsUpdateError::None);
        MassBody live = burn.mMass.mBodies.front();
        CHECK(CapturePhysicsBody(system, id, live));
        CHECK(live.mLinearVelocity[0] == doctest::Approx((48 + 0.09375) / (12 - 1.0 / 128)).epsilon(1e-6));
        CHECK(live.mAngularVelocity[0] > 0);
        CHECK(live.mOrigin[0] > 0);
        // Reject flux from an unresolved finite reservoir rather than discarding spin/inertia.
        auto invalid = ship;
        invalid.mComponents.back().mProperties.mInertia = Mat44::sIdentity();
        PreparedReservoirBurn unchanged = burn;
        CHECK_FALSE(PreparePointReservoirBurn(invalid, delivery, {{7, 103}}, unchanged));
        CHECK(unchanged.mMass.mExhaust.front().mMass == burn.mMass.mExhaust.front().mMass);
        delivery.mEnergyConsumed = 0;
        CHECK_FALSE(PreparePointReservoirBurn(ship, delivery, {{7, 103}}, unchanged));
    }

    TEST_CASE("ReservoirPreviewRejectsOverflowUnderflowAndUnboundedExhaustAtomically")
    {
        const auto state = Ship();
        Thruster engine;
        engine.mID = 7; engine.mMount = Vec4(0, 2, 0, 0);
        engine.mMaxForce = 16; engine.mPropellantRate = 1;
        PropulsionResult valid;
        valid.mPropellantConsumed = 1.0 / 128;
        valid.mEnergyConsumed = 8;
        valid.mActuators.push_back({engine, 0, 0.125});
        PreparedReservoirBurn burn;
        CHECK(PreparePointReservoirBurn(state, valid, {{7, 103}}, burn));
        if (burn.mMass.mBodies.empty() || burn.mMass.mExhaust.empty()) return;
        const auto saved = burn;
        auto rejects_unchanged = [&](const PropulsionResult &bad)
        {
#if defined(JPH_USE_SSE)
            const auto control = _mm_getcsr() & ~0x3fu;
#endif
            CHECK_FALSE(PreparePointReservoirBurn(state, bad, {{7, 103}}, burn));
#if defined(JPH_USE_SSE)
            CHECK((_mm_getcsr() & ~0x3fu) == control);
#endif
            CHECK(burn.mWasteHeat == saved.mWasteHeat);
            CHECK(burn.mMass.mBodies.size() == saved.mMass.mBodies.size());
            CHECK(burn.mMass.mExhaust.size() == saved.mMass.mExhaust.size());
            CHECK(burn.mMass.mBodies.front().mLinearVelocity == saved.mMass.mBodies.front().mLinearVelocity);
            CHECK(burn.mMass.mBodies.front().mAngularVelocity == saved.mMass.mBodies.front().mAngularVelocity);
            CHECK(burn.mMass.mBodies.front().mComponents.back().mProperties.mMass == saved.mMass.mBodies.front().mComponents.back().mProperties.mMass);
            CHECK(burn.mMass.mExhaust.front().mMass == saved.mMass.mExhaust.front().mMass);
            CHECK(burn.mMass.mExhaust.front().mWorldVelocity == saved.mMass.mExhaust.front().mWorldVelocity);
            CHECK(burn.mMass.mLedger.mAddedEnergy == saved.mMass.mLedger.mAddedEnergy);
        };
        auto bad = valid;
        bad.mActuators[0].mDeliveredImpulse = std::numeric_limits<double>::max();
        bad.mActuators[0].mThruster.mMaxForce = 1e-100;
        rejects_unchanged(bad); // J / force overflows before the fuel check.
        bad.mActuators[0].mThruster.mMaxForce = 1;
        bad.mActuators[0].mThruster.mPropellantRate = 1e100;
        rejects_unchanged(bad); // Finite duration * rate overflows.
        bad = valid;
        bad.mActuators[0].mDeliveredImpulse = 1e-300;
        bad.mActuators[0].mThruster.mMaxForce = 1e100;
        bad.mActuators[0].mThruster.mPropellantRate = 1e100;
        rejects_unchanged(bad); // Never recover a flushed-away duration as free impulse.
        bad.mActuators[0].mThruster.mMaxForce = 1;
        bad.mActuators[0].mThruster.mPropellantRate = 1e-100;
        rejects_unchanged(bad); // Positive duration, underflowing consumption.
        bad = valid;
        bad.mActuators[0].mDeliveredImpulse = 1e-100;
        rejects_unchanged(bad); // Positive consumption, zero represented fuel decrement.
        bad = valid;
        bad.mPropellantConsumed = 1;
        bad.mActuators[0].mDeliveredImpulse = 1e200;
        bad.mActuators[0].mThruster.mMaxForce = 1e100;
        bad.mActuators[0].mThruster.mPropellantRate = 1e-100;
        rejects_unchanged(bad); // Charge fits; exhaust speed does not fit MassModel.
        bad = valid;
        bad.mActuators.push_back(bad.mActuators.front());
        bad.mActuators.back().mDeliveredImpulse = std::numeric_limits<double>::max();
        bad.mActuators.back().mThruster.mMaxForce = 1e-100;
        rejects_unchanged(bad); // Even after staging a valid first fuel debit.
        bad = valid; bad.mEnergyConsumed = 0;
        rejects_unchanged(bad); // Accounting rejection also leaves the prior result intact.
        CHECK(state.mComponents.back().mProperties.mMass == 2);
        CHECK(state.mLinearVelocity == MassVector{});
        CHECK(PreparePointReservoirBurn(state, valid, {{7, 103}}, burn));
        CHECK(burn.mMass.mLedger.mExported.mMass == valid.mPropellantConsumed);
        CHECK(burn.mWasteHeat + burn.mMass.mLedger.mAddedEnergy == doctest::Approx(8).epsilon(1e-12));
    }

    TEST_CASE("FrameBoundaryRejectsDoubleMaximumVelocityBeforeNarrowingOrMutation")
    {
        PhysicsTestContext physics(cStep); physics.ZeroGravity();
        auto &system = *physics.GetSystem();
        PreparedPhysicsBody prepared;
        CHECK(Prepare(Ship(), Hull(), prepared));
        auto id = CreatePhysicsBody(system, prepared, Vacuum());
        CHECK_FALSE(id.IsInvalid()); if (id.IsInvalid()) return;
        World world(cStep);
        Frame source, destination;
        source.mID = 1; destination.mID = 2;
        destination.mVelocity[0] = -std::numeric_limits<double>::max();
        CHECK(world.AddFrame(source)); CHECK(world.AddFrame(destination));
        ObjectState object; object.mID = 10; object.mFrame = 1; object.mTier = Tier::Detailed;
        CHECK(CaptureWorldPhysicsBody(system, id, 0, object)); CHECK(world.AddObject(object));
        auto before = system.GetBodyInterface().GetPosition(id);
        CHECK_FALSE(TransferPhysicsFrame(system, world, {{10, id}}, 2, 0));
        CHECK(world.FindObject(10)->mFrame == 1);
        CHECK(system.GetBodyInterface().GetPosition(id) == before);
        CHECK(system.GetBodyInterface().GetLinearVelocity(id) == Vec4::sZero());
    }

    TEST_CASE("FrameBoundaryMovesEntireLiveSolveSetAndSerializesRealInertia")
    {
        PhysicsTestContext physics(cStep);
        physics.ZeroGravity();
        auto &system = *physics.GetSystem();
        World world(cStep);
        Frame a, b;
        a.mID = 1;
        a.mOrigin.mSector = {INT64_C(1) << 50, 0, 0, 0};
        a.mVelocity = {100, 0, 0, 0};
        b = a;
        b.mID = 2;
        b.mVelocity = {102, 0, 0, 0};
        CHECK(a.mOrigin.Translate({10, 0, 0, 0}, b.mOrigin));
        CHECK(world.AddFrame(a)); CHECK(world.AddFrame(b));
        auto first = Ship(), second = Ship();
        second.mID = 20;
        second.mOrigin = {0, 0, 0, 5};
        first.mLinearVelocity = second.mLinearVelocity = {4, 0, 0, 0};
        std::vector<PhysicsBinding> bindings;
        for (auto state : {first, second})
        {
            PreparedPhysicsBody prepared;
            CHECK(Prepare(state, Hull(), prepared));
            auto id = CreatePhysicsBody(system, prepared, Vacuum());
            CHECK_FALSE(id.IsInvalid());
            if (id.IsInvalid()) return;
            bindings.push_back({state.mID, id});
            ObjectState object;
            object.mID = state.mID; object.mFrame = 1;
            object.mTier = Tier::Detailed; object.mRotationRadius = 3;
            CHECK(CaptureWorldPhysicsBody(system, id, 0, object));
            CHECK(world.AddObject(object));
        }
        GlobalPosition before;
        CHECK(world.EvaluateGlobal(10, 0, before));
        // Leaving even one live collider in the old frame is rejected before mutations.
        CHECK_FALSE(TransferPhysicsFrame(system, world, {bindings[0]}, 2, 0));
        CHECK(world.FindObject(10)->mFrame == 1);
        CHECK(system.GetBodyInterface().GetPosition(bindings[0].mBody).GetX() == 0);
        CHECK(TransferPhysicsFrame(system, world, bindings, 2, 0));
        GlobalPosition after;
        CHECK(world.EvaluateGlobal(10, 0, after));
        JPH::SpaceSim::Vector delta;
        CHECK(after.Difference(before, 1, delta));
        for (double v : delta) CHECK(std::abs(v) < 1e-6);
        for (auto binding : bindings)
        {
            CHECK(system.GetBodyInterface().GetPosition(binding.mBody).GetX() == doctest::Approx(-10));
            CHECK(system.GetBodyInterface().GetLinearVelocity(binding.mBody).GetX() == doctest::Approx(2));
            CHECK(world.FindObject(binding.mID)->mFrame == 2);
            BodyLockRead lock(system.GetBodyLockInterface(), binding.mBody);
            CHECK(lock.GetBody().IsCollisionCacheInvalid());
        }
        EntityState entity;
        entity.id = 10; entity.owner = 77; entity.frame = 2;
        CHECK(CaptureEntityPhysicsBody(system, bindings[0].mBody, entity));
        CHECK(entity.mass == doctest::Approx(12));
        CHECK(entity.center_of_mass[1] == doctest::Approx(2.0 / 3));
        CHECK(entity.inertia[0] != 1);
        Checkpoint checkpoint;
        checkpoint.frame.id = 2;
        checkpoint.entities = {entity};
        std::vector<uint8> bytes;
        CHECK(EncodeCheckpoint(checkpoint, Limits(), bytes));
        CHECK(physics.SimulateSingleStep() == EPhysicsUpdateError::None);
    }

    TEST_CASE("LiveMassBoundaryPreservesBodyIDGeometryOriginAndConservativeDetach")
    {
        PhysicsTestContext physics(cStep);
        physics.ZeroGravity();
        auto &system = *physics.GetSystem();
        MassBody ship = Ship();
        ship.mLinearVelocity = {1, -2, 0.5, 0.25};
        ship.mAngularVelocity = {0.1, 0.2, 0.3, -0.1, 0.05, 0.07};
        PreparedPhysicsBody prepared;
        CHECK(Prepare(ship, Hull(), prepared));
        BodyID id = CreatePhysicsBody(system, prepared, Vacuum());
        CHECK_FALSE(id.IsInvalid());
        if (id.IsInvalid()) return;
        CHECK(system.GetBodyInterface().GetUserData(id) == ship.mID);
        MassBody live = ship;
        CHECK(CapturePhysicsBody(system, id, live));
        for (uint i = 0; i < 4; ++i) CHECK(live.mOrigin[i] == doctest::Approx(ship.mOrigin[i]).epsilon(1e-6));
        MassTransaction split;
        CHECK(MassModel::Fragment(live, {{10, {101, 103}}, {11, {102}}}, split));
        if (split.mBodies.size() != 2) return;
        PreparedPhysicsBody retained, debris;
        CHECK(Prepare(split.mBodies[0], {Hull()[0]}, retained));
        CHECK(Prepare(split.mBodies[1], {Hull()[1]}, debris));
        BodyID debris_id = CreatePhysicsBody(system, debris, Vacuum());
        CHECK_FALSE(debris_id.IsInvalid());
        CHECK(InstallPhysicsBody(system, id, retained));
        CHECK(system.GetBodyInterface().GetUserData(id) == 10);
        CHECK(system.GetBodyInterface().GetUserData(debris_id) == 11);
        {
            BodyLockRead lock(system.GetBodyLockInterface(), id);
            CHECK(lock.Succeeded());
            if (!lock.Succeeded()) return;
            CHECK(1 / lock.GetBody().GetMotionProperties()->GetInverseMass() == doctest::Approx(10));
            CHECK(lock.GetBody().IsCollisionCacheInvalid());
        }
        // The detached off-centre box must not remain in the retained live collider.
        RayCastResult hit;
        CHECK_FALSE(system.GetBodyInterface().GetShape(id)->CastRay(
            RayCast(Vec4(-2, 2, 0, 0) - retained.GetShape()->GetCenterOfMass(), Vec4(4, 0, 0, 0)), SubShapeIDCreator(), hit));
        RayCastResult debris_hit;
        CHECK(debris.GetShape()->CastRay(RayCast(Vec4(-2, 2, 0, 0) - debris.GetShape()->GetCenterOfMass(), Vec4(4, 0, 0, 0)), SubShapeIDCreator(), debris_hit));
        CHECK(debris.GetShape()->GetSubShapeUserData(debris_hit.mSubShapeID2) == 102);
        CHECK(physics.SimulateSingleStep() == EPhysicsUpdateError::None);
        MassMomentum total;
        for (auto pair : {std::make_pair(id, split.mBodies[0]), std::make_pair(debris_id, split.mBodies[1])})
        {
            CHECK(CapturePhysicsBody(system, pair.first, pair.second));
            BodyMassData data;
            CHECK(MassModel::Derive(pair.second, data, split.mLedger.mReferenceOrigin));
            for (uint i = 0; i < 4; ++i) total.mLinear[i] += data.mMomentum.mLinear[i];
            for (uint i = 0; i < 6; ++i) total.mAngular[i] += data.mMomentum.mAngular[i];
        }
        for (uint i = 0; i < 4; ++i) CHECK(total.mLinear[i] == doctest::Approx(split.mLedger.mBefore.mLinear[i]).epsilon(1e-5));
        for (uint i = 0; i < 6; ++i) CHECK(std::abs(total.mAngular[i] - split.mLedger.mBefore.mAngular[i]) < 2e-4);
    }
}
