# Mover implementation and acceptance plan

Source baseline: `ecfb0c386e9497fe446905dfab5b2610d49f698e` (`brkdrd/4dJolt`, `4d`).

The two SpaceSim design documents are requirements, not evidence of implemented functionality. This work preserves the existing 4D engine and adds independently testable simulation support under `Jolt/SpaceSim`. The user approved a standalone public `brkdrd/mover` repository with preserved history because the authenticated account already owns the source fork. `upstream` remains the source repository.

## Invariants

- PhysicsSystem remains the only rigid-body integrator and contact-response owner.
- Four spatial coordinates, six bivector components, torque `F wedge r`, and orbital angular momentum `p wedge x` retain the existing conventions.
- No independently editable velocity/momentum duplicates.
- Fixed authoritative time, bounded queues, deterministic event order and explicit failure; never silently skip authoritative time or collisions under overload.
- Floating-point inputs, resource exhaustion, identifier reuse, allocation limits and transaction boundaries require tests.
- No claim of universal CCD, cross-platform determinism or production capacity from a subsystem test or a pyramid average.

## Workstreams and dependency order

1. **Baseline and foundation:** reproduce existing tests and the recorded adversarial sphere-box GJK failure. Repair its cause without weakening the accuracy threshold. Audit gyroscopic mechanics and independent references.
2. **Mass and topology:** aggregate component central second moments and center of mass; distinguish redistribution, expulsion, fragmentation and dissipative rigid merge. Compute validated results before commit.
3. **Propulsion:** bounded ten-axis allocation, response dynamics, actual integrated thrust and resource/heat ledgers; partial-tick shortages must not grant phantom impulse.
4. **World:** sector/offset coordinates, constant-velocity frames, conservative interaction bounds, spatial indexing and trajectory tiers with promotion/demotion guards.
5. **Authority:** validated bounded input, tick stages, ordered single-owner events, snapshots/history and versioned causal state with explicit restore semantics.
6. **Combat:** finite projectile queries against translating/rotating component proxies over an interval; stable IDs and an explicit unresolved-query status at resource limits.
7. **Integration:** connect the modules to real Jolt bodies, run a combined maneuver/combat/topology/frame/input scenario, and verify no competing owner or half-published state.
8. **Validation and publication:** independent review, optimized and assertion/sanitizer builds, existing regression suite, long-duration numerical checks, representative full-run timings and a requirement matrix with unpassed gates called out.

Independent implementation uses isolated Git worktrees. Numerical algorithms and data-integrity-sensitive subsystems receive critical-model agents. The coordinating agent integrates and independently reruns the delivered tests; subagent success claims alone are insufficient.

## Acceptance evidence

Each implementation reports its exact test command and observed result. New behavior follows failing-test / implementation / passing-test cycles. Preserve existing tests and accuracy bounds. Do not disable a failing test to obtain a green report.

| Requirement | Minimum evidence |
|---|---|
| Conservation | Free drift, off-center impulse, asymmetric spin, closed-system topology momentum and energy ledger |
| 4D mechanics | Existing Clifford references, virtual work, full second-moment coupling, COM geometry continuity |
| Propulsion | Symmetric/offset force, ten-axis authority, saturation/failure, partial resource exhaustion, timestep refinement |
| Collisions | Exact known GJK input, moving target, finite projectile, grazing/rotating hull, bounded conservative failure |
| Large world | Extreme sectors/common velocities, frame transforms, interaction-set migration and cross-region query epoch |
| Population | Indexed candidates, promotion horizon, hysteresis, demotion blockers, bounded debris/admission |
| Authority | Sequence/ownership/time validation, reorder/loss, one response per event, consistent bounded history and restore |
| Integrated behavior | A real PhysicsSystem scenario exercising subsystem boundaries together |
| Performance | Versioned workload/build manifest, full-run median/p95/p99/max/deadline misses, correctness metrics |

## Product decisions that remain open in the source specification

Player capacity, population budgets, ship complexity, global world scale, network latency/loss envelope, numerical accuracy budgets and acceptable deadline-miss/headroom policies are explicitly TBD. Prototype tests may use declared conservative configurations, but cannot convert these unknowns into production acceptance. Material damage/penetration rules, a game gravity law, presentation, transport/backend deployment and generalized 4D joints/boundary-mesh collision are not implied by helper APIs. The final report must distinguish implemented/tested capability, integration gaps, and product decisions still required.
