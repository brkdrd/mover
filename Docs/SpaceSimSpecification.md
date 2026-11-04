# 4D space-sim engine specification

Status: **proposed design**, 2026-11-04. This document sets requirements for a PvP space combat simulator built on 4D Jolt. It does not claim the systems below already exist. See the [README](../README.md) for implementation status and [system architecture](SpaceSimArchitecture.md) for integration contracts.

## Product and workload

Players pilot ships through four spatial dimensions, fight other ships, and interact with asteroid fields, stations, projectiles and debris. The initial scope does not require physical player bodies, walking interiors or a general character controller. Ship handling and combat outcomes should follow coherent mechanics, including component failures and changing mass distribution.

Support both detailed ships and large populations of independent objects. Most objects may be moving without contact; collision bursts, ramming, docking and fragmentation remain essential workloads. A dense box pile is a useful solver stress test, but cannot establish the capacity or accuracy of a space-combat simulation.

The proposed foundation is rigid-body ship dynamics plus component systems, scalable spatial queries and authoritative networking. An intact ship is normally one rigid body with many functional components. Independently moving mechanisms and detached pieces may be additional bodies. Continuous structural deformation is optional, with an explicit cost and accuracy budget.

This follows the useful separation documented by DCS: a vehicle model supplies forces and moments to rigid-body/contact physics. BeamNG's point-mass and beam approach is an alternative for structures whose deformation is central to gameplay, rather than a requirement for detailed propulsion or damage. [DCS flight-model description](https://www.digitalcombatsimulator.com/en/support/faq/general/), [BeamNG physics overview](https://www.beamng.com/game/about/physics/).

## Seven engineering priorities

### 1. Conservation laws and stable numerical integration

Requirements:

- Preserve linear momentum and world-space angular momentum in isolated systems within declared numerical tolerances. A torque-free body's angular velocity need not be constant.
- Account for external impulses, gravity, propulsion, recoil and momentum carried by removed or added mass. Conservation applies to the appropriate closed system, not to a thrusting ship alone.
- Track energy according to the model: propulsion performs work, friction and inelastic impacts dissipate energy, and numerical integration must not cause uncontrolled energy growth.
- Use a fixed authoritative timebase and a documented integration scheme, including force sampling, rotational integration, contact solving and any substeps.
- Disable unmodeled drag/damping for free-space acceptance tests. Velocity clamps, sleeping and positional correction must not silently conceal conservation errors.
- Compare timestep refinements over the same simulated duration against analytical or independently computed references. More ticks alone do not prove accuracy.

Acceptance evidence: free drift; off-center impulses; torque-free asymmetric spin; two-body impacts; powered maneuvers; fuel ejection; fragmentation and docking. Report absolute and scale-normalized errors over short events and long runs. Normalize against meaningful reference scales, not near-zero initial momentum.

### 2. Correct 4D rotational mechanics and changing mass properties

Requirements:

- Model four translational and six rotational degrees of freedom, with general simultaneous rotations. Use a valid Spin(4) rotor or equivalent SO(4) representation; six independent scalar angle integrators are insufficient.
- Keep angular velocity, torque, inertia and orientation conventions consistent with point velocity and virtual work. The current repository uses `F ∧ r` for torque; see the [mathematical contract](SpaceSimArchitecture.md#coordinate-and-dynamics-contract).
- Apply the full physical rotational inertia operator, including coupling and gyroscopic dynamics. It maps six-component bivectors to bivectors; an implementation may exploit structure instead of storing a dense 6×6 matrix.
- Derive center of mass and inertia from ship components. Update them consistently when fuel moves, cargo changes, docking occurs or sections detach.
- Distinguish a change of reference point from a physical impulse. Moving a center-of-mass representation must not teleport geometry or introduce motion.

Acceptance evidence: independent Clifford-algebra references; finite-difference point velocities; power identity; analytical mass distributions; asymmetric free spin; random rotations; conservation across mass and topology changes. Positive effective mass alone cannot validate torque signs.

### 3. Component-based forces, propulsion and control

Requirements:

- Define thruster position, direction, maximum force, response time, propellant usage, power draw, heat and failure state. Aggregate their forces and torques about the current center of mass.
- Make flight assistance allocate requests to available actuators. It must not directly assign an ideal velocity or orientation while bypassing physical limits.
- Support four translation commands and six rotation commands internally. Player input and display may expose a smaller context-dependent control scheme without restricting the underlying dynamics.
- Account for actuator saturation, damaged engines, power shortages and changes in mass/inertia. A controller must report or tolerate unachievable requests.
- Integrate fuel, power and heat consistently across different subsystem rates. Charge consumption for actual delivered actuation; do not grant a full tick of thrust when resources run out partway through it.
- Define gravitational and environmental forces explicitly. Four-dimensional geometry alone does not specify the game's gravity law, atmosphere model or unit conventions for fields.

Acceptance evidence: symmetric thrust gives the expected translation; offset thrust gives the expected torque; failed thrusters reduce control authority; fuel/energy budgets match delivered work and modeled losses; control remains bounded under latency and actuator limits.

### 4. Reliable collisions at combat speeds

Requirements:

- Detect relevant intersections over a time interval, including translational and rotational motion. A start/end overlap test is not sufficient for high-speed combat.
- Use conservative broad-phase bounds covering swept motion, followed by narrow-phase tests appropriate to the geometry. Rotation, acceleration bounds and curved trajectories need explicit treatment.
- Support hit positions, times, normals and stable component identifiers in four dimensions. Weapon and ship collision queries must use compatible authoritative geometry.
- Resolve contact through one momentum-response owner. Damage consumes collision results without applying the same impulse again.
- Preserve physically meaningful behavior across large mass ratios. An infinite-mass asteroid is a gameplay approximation, not interchangeable with a finite asteroid receiving recoil.
- Use suitable collision proxies for irregular asteroids and ships. Compound convex proxies are an initial option; general 4D boundary-mesh collision remains separate work.
- Treat hull accuracy, minimum relevant feature size, speed limits, penetration error and CCD failure handling as explicit configuration and acceptance criteria.

Acceptance evidence: grazing hits; high-speed projectiles; rotating long hulls; moving targets; simultaneous impacts; ramming; docking contacts; small fragments hitting massive asteroids. Cover exact adversarial inputs, not only random tests. The existing sphere-box GJK false negative is a blocking robustness issue for this acceptance gate.

Continuous detection can use a time-of-impact search rather than relying solely on smaller fixed steps. This is a general technique, not a claim that the current 4D port has complete rotational CCD. [Continuous collision reference](https://box2d.org/files/ErinCatto_ContinuousCollision_GDC2013.pdf).

### 5. Precision across large distances and velocities

Requirements:

- Separate global location from local contact calculations. A candidate representation is sector coordinates plus high-precision offsets and bounded local simulation frames.
- Preserve small relative motion when bodies share a large translational velocity. High-precision position alone does not solve velocity cancellation.
- Give frames explicit identifiers, origins, velocities and epochs. Initially prefer nonrotating frames; an accelerating or rotating frame requires the corresponding inertial-force terms.
- Make frame changes and cell migration atomic at declared simulation boundaries. Transform positions, velocities, field evaluations, collision caches and network history consistently.
- Avoid making the server's physics origin depend on one client's camera. Different views must observe the same authoritative interactions.

Acceptance evidence: repeat docking, collision and weapon tests near the origin and at extreme offsets/common velocities. Results must agree within the declared budget. Verify rebasing and frame migration do not change momentum, geometry or hit results.

Upstream Jolt documents that velocity storage remains float-based even with double-precision positions. The port's precision path needs its own audit and validation. [Jolt documentation](https://jrouwe.github.io/JoltPhysicsDocs/5.5.0/index.html).

### 6. Large object populations and populated environments

Requirements:

- Index ships, asteroids, stations, projectiles and debris spatially. Avoid all-pairs collision checks and unnecessary detailed contact work for separated objects.
- Distinguish resident objects, active dynamic bodies, potential collision pairs, actual contacts and replicated entities. Object count alone is not a useful capacity metric.
- Allow cheaper trajectory propagation for distant, noninteracting objects under a documented model and error bound. Analytic free drift is exact for force-free translation; arbitrary gravitational motion is not.
- Activate detailed simulation before an interaction becomes possible, using conservative bounds, speed, acceleration, weapon reach and the next scheduling interval. Camera distance alone is insufficient.
- Preserve state and physical consequences across activation, deactivation and region transfers. Never freeze moving asteroids simply because they are invisible.
- Decide which asteroids are finite dynamic bodies, prescribed trajectories or effectively static scenery. State the approximation; promote an object before an impact if gameplay requires it to react.
- Bound debris creation and resource use with explicit game rules. Cosmetic debris must not inconsistently block shots or apply authoritative impulses.

Acceptance evidence: sparse belts, dense local encounters, simultaneous battles, debris bursts and projectiles crossing region boundaries. Report active/resident counts, candidate pairs, contacts, memory, activation latency and approximation error.

### 7. Predictable server cost and reproducible PvP outcomes

Requirements:

- The server owns accepted inputs, physical state, resource consumption, damage and hit outcomes. Clients predict their local experience and reconcile with authoritative snapshots.
- Separate physics frequency, subsystem rates, network snapshots and rendering. All authoritative events still share a simulation timeline.
- Use stable object/event identifiers, controlled randomness, deterministic ordering where required, and versioned replay state. A repeated final hash is useful evidence but not proof of cross-platform determinism.
- If exact deterministic replay is not guaranteed, use state replication and corrections; do not depend on input-only lockstep.
- Define projectile authority, bounded lag compensation, rollback scope, late-input handling and overload behavior. Historical hit queries must not mutate the live world accidentally.
- Profile tail latency, memory and network bandwidth, including impacts, fragmentation, joins and region activation. Do not hide startup or combat bursts behind warmup exclusions.

Acceptance evidence: seeded command replays; thread-count comparisons; packet loss, delay and reordering; bounded correction error; simultaneous collisions and topology changes; overload tests. Validate the exact compiler/configuration combinations claimed to be deterministic.

## Performance and accuracy acceptance

The starting target is a **120 Hz authoritative motion step** (8.333 ms total step budget), with **150 Hz** (6.667 ms) as a stretch target if accuracy and workload permit. These are design goals, not current guarantees. Physics receives only part of the server budget: ship systems, combat, serialization and scheduling also cost time.

Define a versioned benchmark manifest before claiming capacity:

| Required field | What must be specified |
|---|---|
| Hardware and build | CPU, physical/logical cores, worker count, memory, OS, compiler, flags and commit |
| Population | Player ships, hull complexity, components, missiles, resident/active asteroids, active debris |
| Interactions | Contact counts, speed/size and mass-ratio ranges, peak shots/s, breakup sizes, simultaneous fights |
| Time and precision | Integration scheme, timestep, substeps, duration, spatial scale, coordinate/velocity precision |
| Accuracy | Absolute and relative momentum/energy errors, hit/penetration tolerances, trajectory and controller errors |
| Server/network | Snapshot rate, bandwidth, latency/loss envelope, prediction and reconciliation limits |

Player capacity, object budgets and numerical tolerances remain **TBD pending representative prototypes**. They must be assigned before a milestone passes; an unspecified tolerance or workload is not a pass. Do not invent a promise such as “10,000 fully active asteroids at 150 Hz” from the pyramid results.

Report median, p95, p99, maximum, deadline misses and backlog for entire runs, with event markers. Repeat fixed seeds and include long-duration conservation checks. A reliability claim applies only to its declared envelope and observed misses; average throughput is insufficient. Define explicit headroom and acceptable miss policy before production acceptance, and do not lower accuracy under load without declaring the changed model.

## Representative benchmark suite

| Scenario | Primary checks |
|---|---|
| Isolated ship | Free drift, asymmetric spin, propulsion, fuel use, long-run drift |
| Combat maneuver | Full 4D control, saturation, damaged thrusters, recoil |
| Asteroid engagement | Moving irregular obstacles, near misses, occlusion and high-speed hits |
| Docking / ramming | Low-speed precision, high mass ratios, coupled contacts and damage |
| Fragmentation burst | Momentum transfer, component ownership, debris and allocation spikes |
| Multiple battles | Spatial scaling, independent jobs, cross-region projectiles and activation |
| Network replay | Prediction, authoritative hits, loss/reordering, corrections and frame migration |

Keep the existing pyramid as a regression and worst-case contact stress scene. Do not replace its results with a lighter scene and describe that as an optimization.

## Development stages and gates

1. **Foundation audit:** reproduce and repair known collision failures; verify conventions, free-space settings and changing-mass behavior; establish analytical references. Gate: documented correctness envelope with no known failure inside it.
2. **Single-ship vertical slice:** distributed thrust, mass properties, resource limits, control and recorder. Gate: maneuver and conservation checks with bounded errors.
3. **Populated combat space:** irregular asteroid proxies, weapons, swept queries and local contacts. Gate: hit correctness across the speed/size envelope and measured combat cost.
4. **Damage and topology:** component loss, breakup, recoil and docking representation. Gate: transactional state changes with momentum/energy accounting and bounded burst cost.
5. **Large-world scale:** global coordinates, local frames, trajectory tiers and region handoff. Gate: frame-invariance and interaction-activation tests at intended world scale.
6. **Authoritative PvP:** input validation, snapshots, prediction, reconciliation and historical queries. Gate: declared network envelope and reproducible authoritative outcomes. Build serialization/replay hooks from stage 2 rather than retrofitting them here.
7. **Capacity and release:** combined workloads, long runs, overload and recovery. Gate: complete manifest, accuracy evidence, full-run deadline results and documented remaining limits.

These stages organize implementation dependencies. They do not authorize reducing physical fidelity to reach a headline tick rate. The [architecture](SpaceSimArchitecture.md) defines how to connect the stages without multiple systems competing to own the same state.
