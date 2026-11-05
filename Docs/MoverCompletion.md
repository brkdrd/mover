# Mover requirements and implementation boundary

**Status: bounded engine foundation implemented and independently reviewed; final publication checks are in progress. This is not a declaration that every product/specification gate passes.**

The user-approved public repository preserves the history of `brkdrd/4dJolt`; it is a standalone derivative, with the source retained as `upstream`. The implementation plan is [MoverImplementationPlan.md](MoverImplementationPlan.md). Executed results, exact configurations, unresolved findings and compatibility changes belong in [MoverValidation.md](MoverValidation.md), not in inferred feature claims.

The identifiers below enumerate all **42 requirement bullets** under the seven priorities in [SpaceSimSpecification.md](SpaceSimSpecification.md), in their original order. “Implemented” means a bounded implementation exists; it does **not** mean that the unspecified production workload or accuracy gate has passed. “Partial” identifies missing integration or application policy. The final-tree test matrix and outstanding reviews still govern acceptance.

## 1. Conservation and integration

| ID | Implementation / evidence boundary | Status |
|---|---|---|
| 1.1 | Fixed-step core; symmetric Hamiltonian gyro drift; independent free-spin momentum/energy and equal-duration refinement checks. Not a proof for arbitrary contacts or conditioning. | Implemented, bounded |
| 1.2 | External delivery bridge; component mass redistribution, expulsion, fragmentation and rigid merge; explicit recoil/exhaust ledgers. Game gravity is not chosen. | Partial |
| 1.3 | Energy/work/dissipation ledgers and long-spin measurements. Contact energy behavior is measured; a production contact-error budget is absent. | Partial |
| 1.4 | Authoritative integer ticks and fixed PhysicsSystem updates; documented rotational scheme and tick-end query-impact policy. No general within-tick motion-changing event loop. | Partial |
| 1.5 | Free-space fixtures disable damping/sleep; sticky clamp/nonfinite diagnostics and endpoint speed admission fail-stop rather than hide failure. Boundary-impulse and integration-preview repairs passed independent review. | Implemented, bounded |
| 1.6 | Independent Clifford/mechanical references and refinement over the same physical duration, not a comparison of unequal runs. | Implemented, bounded |

## 2. Four-dimensional mechanics and mass

| ID | Implementation / evidence boundary | Status |
|---|---|---|
| 2.1 | Four linear coordinates, six bivector components, proper SO(4)/Spin(4) attitudes and coupled rotation. | Implemented |
| 2.2 | Existing torque `F wedge r`, point velocity and virtual-work conventions retained and tested together. | Implemented |
| 2.3 | Physical inertia represented through a full central second-moment tensor and principal frame, yielding the coupled six-plane operator; gyro dynamics repaired. | Implemented, conditioned domain |
| 2.4 | Component aggregation and separate redistribution/expulsion/fragmentation/docking transactions; adapters install geometry COM and mass together. Float component-mass quantization can explicitly refuse a requested pulse. | Implemented, bounded |
| 2.5 | Recenter and redistribution distinguish coordinate changes from physical impulses. Separate central-spin closure and trap-safe checked mass export repaired the reviewed defects without tolerance relaxation. Subnormal-intermediate accuracy is excluded, and unsupported singular exports may reject. See MoverValidation.md. | Implemented, bounded; reviewed |

## 3. Propulsion and control

| ID | Implementation / evidence boundary | Status |
|---|---|---|
| 3.1 | Thruster mount/direction, force limit, response, fuel, power, heat and health; aggregation about current COM. | Implemented |
| 3.2 | Bounded actuator allocation produces achievable actuation/residuals, not an assigned ideal body velocity. | Implemented |
| 3.3 | Four translation and six rotation request components. No player UI is supplied. | Implemented |
| 3.4 | Saturation, failures, response memory and shared-resource/power limits; actual mass/COM participates in boundary installation. | Implemented, bounded |
| 3.5 | Actual integrated actuation charges resources; depletion cuts do not grant a full interval for free. The live fuel adapter uses a declared instantaneous electric point-reservoir pulse, not continuous variable-mass integration or a complete multi-rate thermal network. | Partial |
| 3.6 | Free-space and explicit external delivery are supported. Mutual gravity, an environmental field law and atmosphere/material models are application decisions, not invented defaults. | Open product policy |

## 4. Collisions and combat

| ID | Implementation / evidence boundary | Status |
|---|---|---|
| 4.1 | Finite-radius projectile queries against moving/rotating compound-box proxies; actual core LinearCast path exercised. General rotating-compound pair CCD is not implemented. | Partial |
| 4.2 | Conservative swept bounds and bounded advancement under the query's stated acceleration/spin model. Arbitrary curved/torque-driven collision trajectories are outside that model. | Partial |
| 4.3 | Four-dimensional positions/normals/times and stable component IDs; integration uses matching sharp-box collision assets. | Implemented, box-proxy domain |
| 4.4 | Physics owns contact response; query-owned impacts have a separate exactly-once path. Damage must not add the same impulse again. Clamp/publication boundary and real replay integration passed independent review. | Implemented, bounded |
| 4.5 | Finite bodies receive momentum; mass transactions retain explicit recoil. A general extreme-mass-ratio solver envelope is not established by these tests. | Partial |
| 4.6 | Core compound convex shapes and compound-box combat proxies. General 4D boundary-mesh collision and automatic proxy generation remain separate. | Partial |
| 4.7 | Explicit query work/precision limits, speed admission and unresolved-query results. The recorded predicate/projection defects are repaired and the unchanged 18,432-query gate passes in the parent's rebuild. Independent review and parent native/sanitizer suites passed; a production geometry/speed/workload envelope remains required. | Bounded gate verified |

## 5. Large-world precision

| ID | Implementation / evidence boundary | Status |
|---|---|---|
| 5.1 | Sector/offset global positions, bounded local solve frames, repaired four-double position/math/stream paths. | Implemented, bounded |
| 5.2 | Constant-velocity reference frames and relative arithmetic; transfers refuse nonzero discarded velocity residuals rather than lose future encounters. Arbitrary boosts are not guaranteed transferable. | Implemented, bounded |
| 5.3 | Frame IDs, origins, velocities and epochs; nonrotating constant-velocity frames only. | Implemented |
| 5.4 | Staged whole-solve-set migration, body poses/velocities, broad-phase refresh and contact invalidation. Pending weapons, external field/history caches and cross-region deployment require caller coordination. | Partial |
| 5.5 | Authoritative frames do not depend on client cameras. No renderer is an owner of physical state. | Implemented architecture boundary |

## 6. Populations and residency

| ID | Implementation / evidence boundary | Status |
|---|---|---|
| 6.1 | World interaction candidates and the engine broad phase avoid detailed solving for separated objects. Large distributed population capacity is not measured by the small core benchmark. | Partial |
| 6.2 | Resident/tier metadata and active native bodies are distinct; benchmark reports bodies, contacts and worker configuration. No replication-interest telemetry/backend yet. | Partial |
| 6.3 | Epoch-based force-free translation with stated precision/refusal bounds. Arbitrary gravitational trajectory approximation is not supplied. | Implemented, free-drift model |
| 6.4 | Swept scheduling bounds include motion, acceleration allowance and weapon reach; explicit horizon limits and preparation-before-promotion. | Implemented, bounded |
| 6.5 | Promotion/demotion guards and whole-set frame-transfer adapters preserve modeled state. General distributed region handoff and all external histories are not automatic. | Partial |
| 6.6 | Detailed/trajectory/distant representation is explicit; finite reacting bodies are used in the integrated scenario. Which world objects may be prescribed/static remains game policy. | Partial |
| 6.7 | Bounded queues/admission and reserved-body fragmentation examples. A game-wide debris lifetime/count/cosmetic-versus-authoritative policy is not chosen. | Partial |

## 7. Authority, replay and cost

| ID | Implementation / evidence boundary | Status |
|---|---|---|
| 7.1 | Authoritative command/event coordinator, snapshots and local reconciliation over real physics. No transport, authentication service, authoritative damage policy or deployed multiplayer server. | Partial |
| 7.2 | Shared integer timeline with independent scheduling/snapshot concepts. Complete multi-rate sensors/AI/power/rendering systems are not implemented. | Partial |
| 7.3 | Stable IDs, bounded ordered events, replay state/versioning and repeat/thread probes. This is not cross-platform bitwise lockstep certification. | Implemented, bounded |
| 7.4 | State-based reconciliation instead of relying on input-only lockstep; presentation smoothing is outside the engine layer. | Implemented boundary |
| 7.5 | Ownership/sequence/intended-time checks, bounded history and local prediction. Coupled-world rollback, solver-cache restoration and production lag-compensation/network envelope remain open. | Partial |
| 7.6 | Reproducible full-run core-Update timing including first calls and collision bursts, with setup/sampling separated. Network bandwidth, joins, distributed activation and full-server tail budgets are not measured. | Partial |

## Architecture boundaries that must not be inferred from helper APIs

- **One owner:** PhysicsSystem remains the integrator/contact-response owner. Mass helpers prepare state; they do not independently advance live bodies.
- **Failure semantics:** preparation rejects without publishing its output. An error after live physics mutation is fail-stop, **not rollback**. Exclusive quiescent world ownership is a precondition for multi-body installation/migration.
- **Event timing:** query impacts are currently handled at a declared tick boundary. General time-of-impact advancement, response, prediction invalidation and remainder integration are not complete.
- **Damage:** failure, detachment and rigid-merge mechanics exist; a general material/penetration/ricochet model and component-connectivity damage processor do not.
- **Fuel:** electrical work, heat and expelled momentum belong to an explicit pulse model. Fuel is not merely subtracted while an unrelated full thrust impulse is also applied.
- **Replay:** causal component/actuator/controller state must come from the image plus immutable assets, never from the current server scene. Local prediction does not restore all coupled solver caches.
- **Precision compatibility:** ObjectStream 2.00 rejects old records. Raw/state streams need application-level versioning; no automatic save migration is supplied.

## Release decisions still required

The specification deliberately leaves player/object budgets, hull complexity, world scale, accuracy tolerances, latency/loss envelope, headroom and acceptable deadline-miss policy undefined. They must be assigned and measured on combined workloads before production acceptance. The bounded benchmark does not replace the inherited pyramid/contact stress measurements or establish a speedup against the original 4D revision.

Remaining numerical defects and review findings must stay visible even when the registered unit suite is green. The broader GJK sweep's previously material failures were repaired without changing its budget; independent review and parent integrated sanitizer/double suites passed. The independently reproduced mass-conservation and checked-export overflow defects were repaired and independently approved within the declared numerical domain. Final publication must state which implementation is available and which gates remain unpassed.
