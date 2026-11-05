# Mover validation record

This is an incremental engineering record, not a production-capacity certificate.
Source baseline: `ecfb0c386e9497fe446905dfab5b2610d49f698e`.

## Resolved review blockers and accuracy boundary

The redistribution central-spin bypass was repaired by `79626115` (integrated `ae9368e7`) and independently approved. Central-spin closure now has a separate gate using the original `1e-9` coefficient and fixed-COM scale; common boosts, orbital representation allowances and caller ledger budgets cannot excuse a failed solve. The original exact counterexample rejects atomically, and valid nonzero-spin/COM-rounding controls remain accepted. The 26-case mass suite passes 1,724,344 assertions, including 4,608 adversarial attempts.

The checked-export overflow was repaired by `dbe75f99` (integrated `d3ffa361`) and independently approved. The original mass-scaled public-Body input now rejects before overflowing diagnostic arithmetic, preserves output and FP control masks, and leaves ordinary solver diagnostics/tolerances intact. Independent probes covered 12,000 symmetric matrices, 26 scaled public-API bodies and a 5,120-export trapped sweep.

The final benchmark output-alias repair also passed independent review: 47 permanent regressions on each of Linux and NTFS through WSL, plus 35 independent probes. Dangling symlink destinations reject before workload execution or artifact mutation. Native Windows benchmark coverage and resistance to concurrent adversarial filesystem changes are not claimed.

**Accuracy limit:** MassModel's subnormal-intermediate exclusion remains explicit. It is a documented accuracy precondition, not comprehensive runtime rejection: a sufficiently tiny positive rounding product can still be reported as zero. The approvals do not establish conservative bounds for every finite input or production acceptance outside the declared numerical domain.

## Current integration status

The sections below headed **Early validation history** retain earlier results; unresolved items there may have been repaired subsequently. This section supersedes those status statements.

- **Latest parent integrated runs:** on `d3ffa361` plus the current benchmark files, freshly rebuilt native precise/FP-traps tests passed **586 cases / 4,566,822 assertions**; whole-engine assertions + ASan/UBSan/leak detection passed **586 / 4,566,725**; optimized double-position passed **588 / 4,566,739**. All three HelloWorld executables also exited zero. The native suite intentionally skips one subprocess-only entry and Linux skips two; their parent tests invoke the children separately. Logs: `out/native-mass-final-{build,tests,hello}.log` and `/home/brkdrd/mover-mass-final-{sanitized,double}-{build,tests,hello}.log`. The rebuilt standalone collision gate again passed all **18,432 unique query IDs** with zero misses, false hits, distance failures or oracle failures. Independent reviews subsequently approved the final central-spin, checked-export and benchmark artifact repairs; the stated accuracy boundaries still apply.
- **Bridge:** final independent review of the pre-integration revision passed on native Windows with FP traps and WSL, each 17 cases / 669 assertions. Combining the new gyro drift exposed an endpoint speed-limit violation after the clamp stage. The parent reproduced it and added endpoint validation; native 17 / 671 passed. Independent native review approved this endpoint follow-up; the query-boundary diagnostic follow-up also passed independent review and is committed as `ea90f3b4`. The separate geometry/reservoir preparation repair also passed independent review after storage recovery: 19 integrated cases / 1,380 assertions, full native suite, and all six original probe modes passed. Validation masking does not extend into PhysicsSystem::Update; numerical failure does not roll back live physics.
- **Dynamics:** independent 600-second/120 Hz checks measured maximum relative world-angular-momentum error `5.12417128e-6` and energy error `3.05582419e-5`. Two-second orientation refinement at 30/60/120 Hz gave `0.00227688`, `0.000584069`, `0.000143805`. Mathematical review passed; the separately discovered AVX2 default-allocator CCD crash was repaired and included in the parent sanitizer run. Independent allocator/endpoint review passed 30 native cases / 28,855 assertions, including real CCD and adversarial allocation probes.
- **Collision:** independent review approved the narrow double-precision pentachoron predicate repair without changing thresholds. Follow-up `099cf64c` repairs the dimensional inconsistency in tetrahedron projection and removes all 329 previously recorded scale-0.01 misses; independent review approved that narrow repair, including an 80-digit oracle and 6,000 projection/fallback assertions. The parent native full run also passed these registered regressions (its sole failure was the separately tracked integrated fuel assertion). The broader 18,432-query sweep historically had material misses, including 382/512 selected scale-0.01 cases. Lower-dimensional repair `b4e5ef3e` (integrated `81161f8b`) now resolves **all 2,154 material misses and 8,733 distance failures**, with zero new failed identities. The parent unchanged standalone gate exits 0 with no oracle failures. Independent review approved the repair without tolerance/termination/legacy-allowance changes, and parent native 577-case plus whole-engine sanitizer 576-case runs passed. This verifies the bounded envelope, not arbitrary collision/CCD correctness.
- **Mass core:** final checked-export follow-up `d3ffa361` passed independent review. The original finite mass-scaled state (`M = x*x^T + 4I`, `x=(1,2,3,4)`, initial mass 2, scaled to `ldexp(2.f,70)`) now rejects atomically without FP faults. Accepted-range response checks remain within the unchanged `64 * FLT_EPSILON` output policy. Ordinary eigensolver diagnostics are preserved; arbitrary oblique singular round trips remain intentionally unsupported.
- **Double precision:** independent review approved four-lane math/serialization/Morton repairs. Isolated double build: 383 cases / 2,473,442 assertions and HelloWorld passed. Independent assertions/FP-trap probes exercised SSE2, SSE4.1, AVX, AVX2, AVX512 and forced scalar. NEON/RVV were inspected, not executed. Inherited AVX512 signed-zero differences remain; unrestricted backend parity is not claimed. Legacy/intermediate-operation coverage has been restored. The parent native DVec4 suite passed 16 cases / 490 assertions after locally scoping a deliberate 0/0 fixture; the isolated double suite includes an additional precision-specific case.
- **Serialization compatibility:** ObjectStream is now **2.00** and rejects legacy records. Raw/state-stream layouts also change and require enclosing application versioning; no automatic migration is supplied.
- **Module review repairs integrated:** mass `4f1d71ef8999d00afb68ee7e27778b5ff9985de4`, propulsion `cf1f791cae739a049ab1a408501e1c0f1ff3a838`, world `d7b1ee5ec801ee02fdcaebaa635c69bc4abf2617`, authority `1171441305491934523099af6a4aaeb06ecc3ff1`. Authority and propulsion/world numerical semantics passed independent re-review. Final native World validation review passed 31 cases / 6,240 assertions plus a 50,196-pair independent oracle; scoped Propulsion guards were also approved. These trap-safe changes are committed as `b3e14da7`. The final mass redistribution follow-up also passed review as described above.

Packaging/configuration review passed and is committed as `e2a5e62e`. It verified all three Linux workflow modes, both standalone diagnostics, native Windows configuration, and byte-identical installation of all seven SpaceSim headers; the external double-position SDK consumer also ran successfully. `Build/CMakePresets.json` is explicitly tracked despite the inherited `/build/` ignore rule. This is local configuration and execution evidence, not a hosted GitHub Actions result.

Sanitizer configuration: GCC Release, precise math, LTO off, assertions on, `OVERRIDE_CXX_FLAGS=OFF`, `CMAKE_CXX_FLAGS_RELEASE="-O1 -g"`, `-fsanitize=address,undefined -fno-omit-frame-pointer`. Runtime: `ASAN_OPTIONS=detect_leaks=1:halt_on_error=1`, `UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`. The entire linked engine was instrumented, not just headers.

Native precise math avoids two inherited fast-math EPA failures. Presets/CI disable C4866/C4868 evaluation-order diagnostics and C5246 standard-conforming aggregate brace-elision warnings; numerical assertions/FP exceptions remain enabled. The build option `CROSS_PLATFORM_DETERMINISTIC` is not evidence of cross-platform simulation determinism.

The combined gameplay scenario, parent integrated matrix and scoped independent reviews passed. The benchmark was rebuilt from clean source `e503e2db` and all eight runs / 9,600 steps passed the unchanged checks; raw artifacts and dense-pile timings are archived in MoverPerformance.md. Hosted CI passed on implementation head `88a0925d`: both push run `41966332818` and PR run `41966456438` completed Linux Release, double-position, ASan/UBSan/assertions, and Windows MSVC jobs successfully. PR #6 was merged into `main` as `ca21b694`, preserving history. The final artifact audit independently reproduced the archived timing summaries, row counts and fingerprints without rerunning the benchmarks. These are bounded implementation gates, not production workload/network/accuracy acceptance; isolated branch counts must not be summed into a final result.

## Early validation history

### Baseline and discovered failures

A fresh WSL GCC 13.3 Release build passed 389 cases / 2,397,320 assertions and HelloWorld. Nevertheless, targeted independent probes exposed:

- The exact recorded sphere-box GJK false negative (~0.105954 overlap), due to cancellation in the float pentachoron containment predicate.
- An improper principal-inertia eigenframe (determinant -1) passed through a proper-rotation Rotor conversion.
- A wrong-sign gyroscopic update, with world angular-momentum error that did not vanish under timestep refinement.
- Incorrect solid-hypersphere second moments and incomplete mass/inertia export.
- Remaining three-dimensional assumptions in double-coordinate math/serialization and Morton indexing.

These findings illustrate why the pre-existing suite alone was not foundation acceptance. No tolerance relaxation is an acceptable substitute for repairing them.

## Reviewed standalone capabilities

### Continuous projectile/compound-box queries

Implementation: `8faf85625a5f513f50e32744e3c7bd843a9be1b6` (integrated as `ae27a6ed`).

The coordinating agent reran `UnitTests --test-suite=CombatTests`: **19 cases / 864 assertions passed**. An independent critical review found no blocking issue inside the documented contract. It additionally exercised an independent double Clifford/matrix reference, constructed collisions, extreme cases, and sanitizers. Review approval applies only to the standalone query, not full simulation event ownership.

The contract is constant world acceleration/spin over a bounded interval, authoritative sharp 4-box components, finite-radius projectiles, declared numerical allowances and bounded query work. `Indeterminate` is not a miss. General compound-pair rotational CCD and within-tick physical event integration are not supplied by this helper.

### Physics bridge

The coordinating agent reran the first numerical-fix revision: **15 cases / 487 assertions passed** under WSL Release. Independent review of that revision and native validation are still pending.

The initial four-test adapter was rejected in review: exact float conversion could overflow after validation, and gyroscopic evolution could trigger an unreported angular clamp. Permanent regressions reproduce both. Sticky core diagnostics and exact validated delivery storage were added. Numerical failure means fail-stop without tick/snapshot publication, not rollback of the live core.

Native MSVC 19.50.35721 builds with narrowly disabled inherited warnings C4866/C4868 (evaluation-order diagnostics on doctest expressions and independent braced constructor arguments). One ordinary bridge test passed. The complete bridge suite exits abnormally with floating-point exceptions enabled; this is an unresolved validation issue, not a pass.

## Integrated build matrix

Pending integration and independent execution:

- GCC optimized headless, complete existing and new suite, HelloWorld.
- GCC double-position build, complete suite.
- Assertions plus ASAN/UBSAN, complete linked implementation.
- Native MSVC, complete suite and HelloWorld.
- Combined ship/resource/combat/topology/frame/authority scenario.
- Fixed-seed whole-run timings with median/p95/p99/max, misses and numerical metrics.

The GitHub workflow is being expanded to exercise these compiler/configuration boundaries. A workflow definition is not evidence that CI passed.

## Explicit unpassed product gates

The source specification leaves production population/workload, numerical tolerances, network envelope, headroom and acceptable deadline-miss policy TBD. Prototype envelopes and observed test bounds must be stated separately. No 120/150 Hz production reliability, cross-platform deterministic lockstep, generalized collision envelope or complete authoritative game-server readiness is claimed here.
