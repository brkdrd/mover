# Pentachoron predicate precision repair

## Scope

Baseline: `ecfb0c386e9497fe446905dfab5b2610d49f698e`. The exact Windows-generated sphere–box pair has an analytic overlap of approximately 0.105954, yet raw GJK returned false. Its five-point support simplex strictly contains the origin: independently solved barycentrics are `[0.29662915, 0.30951689, 0.28358813, 0.09363796, 0.01662787]`. Float Gram/cofactor cancellation made the last numerator negative, incorrectly excluding the interior.

Only `ClosestPoint::OriginOutsideOfPentachoronPlanes` changes in production. Dot products promote coordinates **before** multiplication/accumulation; Gram entries, cofactors, determinant, numerators and sign comparisons remain double. Edge subtraction and input/output vectors stay float. The `1.0e-6f` Hadamard degeneracy threshold and `FLT_EPSILON` barycentric slop are unchanged, as are all other collision tolerances, GJK convergence rules, tetrahedron/triangle algorithms, and existing statistical accuracy bounds.

This is a bounded cancellation repair, **not** an exact-arithmetic predicate or a general robust-simplex redesign.

## Permanent regressions

- `TestGJKSphereVsBoxKnownFalseNegative`: exact input bits, independent double analytic positive-overlap check, raw intersection and closest-distance checks at the existing `1e-4` tolerance; both argument orders and zero/X/W/mixed initial axes.
- `TestPentachoronKnownInteriorCancellation`: exact support bits, all 120 vertex permutations at scales 0.01, 1 and 100. Require mask zero, full feature set 31, and exactly zero closest point for both unrestricted and `MustIncludeE` paths.
- `TestPentachoronSliversAgainstKnownBarycentrics`: 48 seeded rotations/weight choices, all 120 permutations at five scales (0.001–1000), in two subcases. Known positive/negative barycentrics test interior and every outside-face bit; near-face weights are ±1/32768. The unrotated worst anchor has Gram/Hadamard ratio about `3.64e-6`, above the unchanged `1e-6` threshold. There are 57,600 direct mask checks, plus interior feature/closest-point checks. The oracle is prescribed barycentric coordinates, not a second Gram solve. The existing explicitly degenerate fixture remains unchanged.
- `TestGJKSeededRotationsAndScalesAgainstAnalyticIntersections`: 64 seeded configurations, scales 0.1/1/100, four initial axes, both argument orders, both overlap/separation, sphere–box and sphere–sphere. All 6,144 raw intersection queries must agree; no disagreement allowance. Prescribed corner distance and center distance give analytic references, with signed margin `0.125 * scale`; query tolerance remains `1e-4 * scale`.

Random components are drawn sequentially from explicit `mt19937` bits, avoiding CRT distribution and constructor-argument evaluation-order dependencies in the new seeded tests.

## TDD and executed validation

The exact GJK and support-simplex regressions were written and run before changing the production header. Original code returned mask 1 rather than 0, feature 15 or fallback vertex 1 rather than 31, and GJK squared distance `0.000126710147` rather than zero.

```text
RED (original production code):
[doctest] test cases:    2 |    0 passed |  2 failed | 389 skipped
[doctest] assertions: 1821 | 1787 passed | 34 failed |
[doctest] Status: FAILURE!

GREEN (same exact regressions after the repair):
[doctest] test cases:    2 |    2 passed | 0 failed | 389 skipped
[doctest] assertions: 1821 | 1821 passed | 0 failed |
[doctest] Status: SUCCESS!

FINAL FULL SUITE (both Release and assertions + ASAN/UBSAN):
[doctest] test cases:     393 |     393 passed | 0 failed | 0 skipped
[doctest] assertions: 2482950 | 2482950 passed | 0 failed |
[doctest] Status: SUCCESS!
```

Additional standalone header comparisons run the final exact + seeded simplex tests against a pinned original-header overlay: 2,842 assertions fail with the original header; all 79,083 pass with the repair. The original exact reproducer now reports `Intersects=1`, closest squared distance zero, mask zero and feature 31. Release and sanitizer HelloWorld both report `4D HelloWorld: all checks passed`, sleeping after 63 steps. Builds use `-j1`; no compiler warnings/errors in the final build logs.

Environment: WSL Ubuntu GCC 13.3.0, x86-64 AVX2, `-ffp-contract=off`, headless and profiler-disabled, **LTO disabled**. Release uses `-O3`. The separate sanitizer build uses `USE_ASSERTS=ON`, `-O1 -fsanitize=address,undefined -fno-omit-frame-pointer`, leak checking and halt-on-UB. Native MSVC validation is still required; this report does not claim it was run.

Build location: `/home/brkdrd/mover-build-collision` (sanitizers in its `sanitized` subdirectory).

```bash
cmake -S /mnt/c/Users/goliy/mover-worktrees/collision/Build \
  -B /home/brkdrd/mover-build-collision -DCMAKE_BUILD_TYPE=Release \
  -DDEBUG_RENDERER_IN_DEBUG_AND_RELEASE=OFF -DPROFILER_IN_DEBUG_AND_RELEASE=OFF \
  -DINTERPROCEDURAL_OPTIMIZATION=OFF -DTARGET_SAMPLES=OFF -DTARGET_VIEWER=OFF \
  -DTARGET_UNIT_TESTS=ON -DTARGET_HELLO_WORLD=ON
cmake --build /home/brkdrd/mover-build-collision --target UnitTests HelloWorld -j1
/home/brkdrd/mover-build-collision/UnitTests
/home/brkdrd/mover-build-collision/HelloWorld
```

## Remaining small-scale failure: not hidden by the green suite

An exploratory version of the seeded GJK test also included scale **0.01**. It produced **329 sphere–box false negatives**, all for overlapping pairs with prescribed overlap `0.00125` and query tolerance `1e-6`. There were no failures at 0.1/1/100, no sphere–sphere failures, and no false positives in that probe.

The same 329 failed identities occur with **both the pinned original header and repaired header**. These failures are pre-existing and not repaired by the pentachoron change. Their detailed root cause has not been isolated here; existing lower-dimensional float projection/absolute degeneracy handling warrants a separate investigation. The default permanent seeded GJK test explicitly covers the passing 0.1/1/100 envelope, with a source comment pointing here. Its tolerance was not increased to excuse the 0.01 failures. The direct pentachoron tests still cover 0.001 and 0.01.

Preserved executable evidence under the build's `artifacts` directory:

- `GJKTests-small-scale-diagnostic.cpp`, `UnitTests-small-scale-diagnostic`, `small-scale-diagnostic.log`;
- standalone `gjk-diagnostic-original` / `gjk-diagnostic-fixed` and corresponding logs;
- `collision-validation.py` with compiler commands and baseline pinned by SHA;
- `additional-verification.json`, exact repro source/log, red/green/final suite logs, sanitizer logs.

Run either standalone diagnostic with `--test-case=TestGJKSeededRotationsAndScalesAgainstAnalyticIntersections` to reproduce its expected nonzero exit. To regenerate the diagnostic on a different machine, copy the permanent GJK test source to a scratch runner and add `0.01f` to that test's scale list; do not weaken the assertions.

A separate exploratory variant using double cofactors but already-rounded float dot products also passes the current regression corpus. That does **not** validate this shortcut over a wider envelope; the retained implementation evaluates the dot products themselves in double. The tests are evidence for the measured cases, not proof that every lower-precision variant is excluded.

Consequently this change closes the recorded 0.105954 false negative, but **does not close the SpaceSim collision-robustness acceptance gate** or establish general small-scale GJK correctness.

## Performance sanity, not a server-capacity claim

Same compiler/flags, original versus fixed header, AMD Ryzen 5 9600X, profiler and LTO off. An external benchmark uses 4,096 seeded inputs, non-inlined entry wrappers, checksum-consumed results, one million predicate calls and 200,000 sphere–box queries per sample. ABBA ordering yields 14 samples per variant:

| Median ns/call | Original | Fixed |
|---|---:|---:|
| Direct pentachoron predicate | 30.57 | 32.78 |
| Mixed random sphere–box GJK | 118.96 | 110.68 |

Predicate cost is about 1.07× in this probe. GJK timing ranges overlap substantially (original 104.32–148.70 ns; fixed 103.89–130.31 ns); concurrent host activity means the apparent GJK speedup is **not** a reliable optimization claim. This only rules out an obvious order-of-magnitude regression in the tested workload. It is not whole-scene, tail-latency, or 120 Hz acceptance evidence. Source, compiler commands and every measurement are preserved in `predicate-benchmark.cpp`, `collision-benchmark.py`, `benchmark.json` and raw logs.
