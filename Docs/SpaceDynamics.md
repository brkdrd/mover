# Opt-in free-space rotational dynamics

This documents the rigid-body core change, not complete space-sim or collision/CCD acceptance. The broader requirements remain in [SpaceSimSpecification.md](SpaceSimSpecification.md) and [SpaceSimArchitecture.md](SpaceSimArchitecture.md).

## Configuration and ownership

The existing `BodyCreationSettings::mApplyGyroscopicForce` / `Body::SetApplyGyroscopicForce` opt in to the repaired integrator. The default remains **false**, and the default contact workload keeps its old integration path. A free-space body should explicitly use:

```cpp
settings.mAllowedDOFs = EAllowedDOFs::All;
settings.mApplyGyroscopicForce = true;
settings.mLinearDamping = 0.0f;
settings.mAngularDamping = 0.0f;
settings.mAllowSleeping = false;
settings.mGravityFactor = 0.0f; // Or set PhysicsSystem gravity to zero.
// Choose explicit finite application limits with sufficient headroom.
// Conservation tests use 1e6 for both limits; no tested trajectory approaches them.
settings.mMaxLinearVelocity = 1.0e6f;
settings.mMaxAngularVelocity = 1.0e6f;
```

Use a fixed step, initially 1/120 s. Limits are safety policy, not numerical stabilization. An asymmetric body's angular speed can increase as its orientation changes while its energy stays constant. A limit hit changes physical momentum/energy; reject/report such a run rather than calling it a conservation pass. Do not add an application-side angular-momentum accumulator that can disagree with the body's authoritative angular velocity.

The free-top path requires all six rotational DOFs and finite positive principal inertias. Locked rotational DOFs or zero inverse inertias use the constant-angular-velocity orientation path even if the flag is set; they are **not** claimed to conserve isolated-top momentum. The caller must validate finite mass/inertia and normalized initial orientation. Sleeping, damping, gravity, driven kinematics and mass changes are different physical models.

## Integration scheme

Let `Q = body_rotation * inertia_rotation` map principal coordinates to world coordinates. The stored bivector convention is `[XY, XZ, XW, YZ, YW, ZW]`, point velocity is `Omega.Contract(r)`, and torque is **F wedge r**, not r wedge F.

For the torque-free rotational drift, compute a step-local world momentum from the current orientation and angular velocity. In principal coordinates, `H(L) = sum_k L_k^2 / (2 I_k)`. Each individual Hamiltonian `H_k` has an exact flow:

```text
angle = dt * L[k] / I[k]
Q' = Q exp(0.5 * angle * e_k)
L' = Ad(exp(-0.5 * angle * e_k)) L
```

`L[k]` is constant during its own plane flow. Compose planes 0..4 for half a step, plane 5 for a full step, and planes 4..0 for half a step. This **11-stage symmetric Hamiltonian (Strang/Lie-Poisson) splitting** is second order and time reversible in exact arithmetic. Unlike forward Euler plus length rescaling, it advances orientation and momentum on the same coadjoint orbit. Energy is not projected to its initial value: it has a bounded splitting error within the tested step/spin envelope. This fixed-work scheme avoids an implicit midpoint convergence/failure policy.

At the endpoint, derive world angular velocity from the original step-local world momentum and the **actual normalized endpoint body orientation**. No momentum is stored across updates. Float rotor adjoints are not exact inverses: repeated forward/inverse applications were measured to amplify momentum even for zero time. The opt-in drift therefore evaluates normalized SO(4) adjoints and endpoint momentum/velocity conversion in double precision, rounding the final angular velocity once. This corrects rotation representation error; it never rescales momentum or energy to conceal an integration error.

The existing pipeline has one owner per kick:

1. `ApplyForceTorqueAndDragInternal` applies external forces/torques and damping once at the starting orientation. The force kick retains its existing first-order sampling and safety limits.
2. Existing warm starting/contact velocity solves apply their impulses once.
3. `PhysicsSystem::JobIntegrateVelocity` calls `Body::IntegrateRotationInternal(dt)` once. Dynamic opt-in bodies call `MotionProperties::IntegrateGyroscopicRotationInternal`, which returns the endpoint rotor and updates angular velocity together. Other bodies retain the existing exponential step.
4. Existing CCD impulses and position corrections keep their existing owners. `AddRotationStep` / `SubRotationStep` remain geometric position-correction operations; they do not integrate free motion again.

There is no earlier standalone gyroscopic torque kick. The internal `ApplyGyroscopicForceInternal` entry point has been replaced, not kept as an independently callable second integrator. Applications should use `PhysicsSystem::Update`, not either internal integration entry point.

The free drift is second order; the entire externally driven/contact/CCD update is **not** thereby promoted to second order. Position projection can alter angular momentum; conservation across it is not promised. The opt-in contact regression disables position projection to isolate the actual velocity impulses. Existing floor/stack/pyramid regressions still exercise the default pipeline and its position correction.

## Prerequisite inertia-frame repair

`MassProperties::DecomposePrincipalMomentsOfInertia` sorted its eigenvectors and then used `sCross(c0,c1,c2).Dot(c3)` as a determinant. That expression has the **opposite sign** under this repository's 4D cross convention. It forced determinant -1 frames, which cannot be represented by a Spin(4) rotor. The correction uses `Mat44::GetDeterminant()` before optionally negating the final eigenvector. Reconstructing the tensor before conversion to a rotor was insufficient to catch this: a reflection also diagonalizes a symmetric tensor.

The regression checks proper orientation with an independent Leibniz determinant and reconstructs physical second-moment tensors **after the eigenframe -> rotor -> matrix round trip**, for `[1,2,4,8]` principal second moments and 64 seeded general rotations plus the diagonal case.

## Evidence and practical envelope

`UnitTests/Physics/SpaceDynamicsTests.cpp` drives the real `PhysicsSystem`, with zero gravity/damping, sleeping disabled, and distant inactive safety limits. Reference integration is a separate double-precision SO(4) matrix RK4 solver. It computes inertia action with matrix products, not production rotor adjoints, commutator formulas, exponential stepping or `PredictOrientation`. RK4 reference steps 1/3840 and 1/7680 s must agree to 1e-10 over the same two seconds.

Before production changes, the original implementation failed both free-spin cases: maximum relative world-momentum error reached 1.61--1.65 over 60 s at 60/120/150 Hz, and the two-second orientation error stayed near 2.994 under timestep refinement. The inertia-frame regression separately failed with determinant -1 and tensor reconstruction errors after rotor conversion. A further zero-time regression caught float-adjoint amplification (angular-velocity error 8.81e-4 over 1000 zero-time drifts) before the double-conversion correction.

The principal long-run sample is a box with half extents `(0.35,0.65,1.05,1.6)`, a general initial orientation, all six angular-velocity components `(1.1,-0.7,0.9,0.6,-1.2,0.8)` rad/s, and principal inertia ratio about 6.72. Acceptance bounds are maximum relative momentum error < 1e-4 and maximum relative rotational-energy error < 2e-4. The 120 Hz run lasts 600 s; 60 and 150 Hz runs last 60 s. The two-second refinement compares 30, 60 and 120 Hz against the same independent reference, requiring each halving to reduce orientation error by at least 2.5x. Eight further physical tensors with independently rotated `[1,2,4,8]` second moments exercise nontrivial inertia frames.

Additional checks cover the double adjoint against brute-force Clifford sandwiches, a commuting asymmetric double rotation against a Clifford power-series analytical solution, off-center force/impulse signs and one-time force clearing in **every** rotational plane, and actual opt-in asymmetric-body contact impulses.

Observed WSL Release results (GCC, profiler off, LTO off, `-j1`):

| Rate | Duration | Max relative world L error | Max absolute world L error | Max relative energy error | Max absolute energy error |
|---|---:|---:|---:|---:|---:|
| 60 Hz | 60 s | 2.01733e-6 | 0.0211643 | 9.87097e-5 | 1.02934 |
| 120 Hz | 600 s | 5.12417e-6 | 0.0537590 | 3.05582e-5 | 0.318659 |
| 150 Hz | 60 s | 1.72852e-6 | 0.0181343 | 1.67652e-5 | 0.174826 |

Absolute quantities use the test's consistent mass/length/time units. Relative errors use the nonzero initial angular-momentum norm and rotational energy. The two-second orientation Frobenius errors are **0.00227688, 0.000584069, 0.000143805** at 30/60/120 Hz; the angular-velocity errors are **0.00144988, 0.000376608, 0.0000886243** rad/s. These are measured independent-reference errors, not differences against the production integrator itself.

The focused suite passes **9 cases / 176,780 assertions**; the full regression passes **398 cases / 2,574,100 assertions**, and HelloWorld passes all checks. This includes the original 389 cases, not a filtered replacement suite. Native Windows, sanitizer, cross-architecture determinism and ship-population timing were not run in this workstream.

These are validated samples, **not** a universal guarantee for all spin rates, inertia ratios or durations. Very ill-conditioned inertia, large per-tick rotation, locked rotational DOFs, mass/topology changes, positional projection and rotational CCD need their own acceptance envelopes. Increase temporal resolution and repeat same-duration refinement rather than changing conservation tolerances or allowing clamps to hide errors.

## Build and cost

```sh
cmake -S Build -B /home/brkdrd/mover-build-dynamics \
  -DCMAKE_BUILD_TYPE=Release -DINTERPROCEDURAL_OPTIMIZATION=OFF \
  -DDEBUG_RENDERER_IN_DEBUG_AND_RELEASE=OFF \
  -DPROFILER_IN_DEBUG_AND_RELEASE=OFF \
  -DTARGET_UNIT_TESTS=ON -DTARGET_HELLO_WORLD=ON \
  -DTARGET_PERFORMANCE_TEST=OFF -DTARGET_SAMPLES=OFF -DTARGET_VIEWER=OFF
cmake --build /home/brkdrd/mover-build-dynamics --target UnitTests HelloWorld -j1
/home/brkdrd/mover-build-dynamics/UnitTests --test-suite=SpaceDynamicsTests
/home/brkdrd/mover-build-dynamics/UnitTests
/home/brkdrd/mover-build-dynamics/HelloWorld
```

The opt-in path costs eleven plane exponentials/adjoints plus two small double-precision adjoints, rather than one constant-velocity exponential. Both adjoints and world momentum are temporary; body layout, serialization and contact-row layout do not gain momentum state. Default gyro-off workloads do not pay this solver cost. No capacity or frame-time claim follows from these correctness tests; measure the intended ship workload before enabling the option broadly.
