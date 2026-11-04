# Collision predicate and tetrahedron projection precision repairs

## Scope

Baseline: `ecfb0c386e9497fe446905dfab5b2610d49f698e`. The exact Windows-generated sphere–box pair has an analytic overlap of approximately 0.105954, yet raw GJK returned false. Its five-point support simplex strictly contains the origin: independently solved barycentrics are `[0.29662915, 0.30951689, 0.28358813, 0.09363796, 0.01662787]`. Float Gram/cofactor cancellation made the last numerator negative, incorrectly excluding the interior.

The initial repair (`b2084e89`) changes only `ClosestPoint::OriginOutsideOfPentachoronPlanes` in production. Dot products promote coordinates **before** multiplication/accumulation; Gram entries, cofactors, determinant, numerators and sign comparisons remain double. Edge subtraction and input/output vectors stay float. The `1.0e-6f` Hadamard degeneracy threshold and `FLT_EPSILON` barycentric slop are unchanged, as are all other collision tolerances, GJK convergence rules, tetrahedron/triangle algorithms, and existing statistical accuracy bounds.

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

## Small-scale follow-up: recorded 329 failures repaired, broader gate still fails

Baseline for this follow-up: `b2084e89`. Restoring scale **0.01** to the permanent seeded test reproduces all **329 sphere–box false negatives**, with prescribed overlap `0.00125` and query tolerance `1e-6`. These are material overlaps, not grazing cases within the declared tolerance.

### Isolated cause and bounded production change

For sample 0, initial axis W, sphere first, the last tetrahedron has float Gram determinant `3.13343101466e-13` and determinant/Hadamard ratio `0.00084322324`. Its face predicate accepts this interior with the existing relative `1e-6` conditioning threshold, but the four-point `GetBaryCentricCoordinates` overload rejects it against an **absolute `1e-12` determinant cutoff**. This determinant has units length^6. Falling back to a triangle returns a squared distance `1.42379761314e-7`, slightly worse than the preceding `1.42379732893e-7`, so GJK exits as no-progress/false. The independent 65-digit projection instead has squared distance `2.93892382265e-8` and strictly positive barycentrics `[0.00832342377, 0.37596215505, 0.49389370775, 0.12182071343]`.

Only that **four-point barycentric overload** changes in this follow-up:

- Replace the dimensionful absolute cutoff with `det <= 1.0e-6f * d00 * d11 * d22`, matching the existing tetrahedron face predicate's scale-invariant conditioning policy. Rank-deficient tetrahedra still fall back; no geometric/query tolerance is increased.
- Promote coordinates before dot multiplication/accumulation, and retain double through cofactors, determinant and the solve; explicitly round only the returned float weights. The relative-only change fixes the 329 GJK failures, but still fails 48 exact-support projection assertions under vertex permutation/scaling at the fixed regression precision. The double solve independently fixes those cancellation errors.

No change to `GJKClosestPoint.h`, support mappings, triangle/line algorithms, tetrahedron face predicates, pentachoron predicates, GJK termination rules, `cCompareTolerance`, or `cMaxDisagreements`. This is not a shape-specific GJK workaround or a general robust-simplex rewrite.

### Permanent regressions and RED/GREEN evidence

- `TestTetrahedronSmallScaleProjection`: prescribed centroid weights 1/4, all 24 permutations, scales 0.001/0.01/0.1/1/100, both closest-point paths. Independent analytical projection is `(0, 0, 0, 0.125 * scale)`.
- `TestTetrahedronSmallScaleKnownSupportProjection`: exact support bits from the trace, independent 65-digit reference, all 24 permutations and scale multipliers 0.1/1/10/100. Both closest-point paths must match within `1e-7 * multiplier`; this is tighter than the original query tolerance, not an increased allowance.
- `TestGJKSeededRotationsAndScalesAgainstAnalyticIntersections`: scale 0.01 is now permanently included alongside 0.1/1/100. All **8,192 raw intersection queries** agree with the prescribed analytic references, without any disagreement allowance. Query tolerance stays `1e-4 * scale`.

The tests ran on the unchanged follow-up baseline before the repair: the seeded test had 329 failed assertions; the two direct tetrahedron tests had 336 failed assertions. With only the relative cutoff, the seeded test passes and the direct tests retain 48 failures. With the final double solve, all 8,704 seeded-test assertions and all 1,776 direct-test assertions pass.

Final rebuilt **Release and assertions + ASAN/UBSAN**, each with `-j1`, both report:

```text
[doctest] test cases:     395 |     395 passed | 0 failed | 0 skipped
[doctest] assertions: 2486897 | 2486897 passed | 0 failed |
[doctest] Status: SUCCESS!
```

Both rebuilt HelloWorld executables pass all checks and sleep after 63 steps. The build configurations and compiler remain those documented above. **Native MSVC/precise validation is not claimed.**

### Bounded remaining-envelope investigation: material failures remain

A separate fixed sweep checks 64 deterministic rotations, four initial axes, both argument orders, six scales, and signed margins `±0.125`, `±0.001`, `±0.00005` times scale: **18,432 sphere–box queries per variant**, baseline versus final. The double reference uses the actual float sphere center/radius and the transformed box corner, and verifies that the residual lies in that corner's outward normal cone. Thus it does not silently assume exact rotor/matrix orthogonality.

For each scale, 1,024 of 3,072 queries lie within `abs(signed overlap) <= query tolerance`; these are reported separately and are **not counted as material errors**. Each positive material margin has 512 queries. Final false negatives are:

| Scale | Overlap `0.125 * scale` | Overlap `0.001 * scale` (about 10 query tolerances) |
|---|---:|---:|
| 0.001 | 450 | 502 |
| 0.003 | 194 | 499 |
| 0.01 | **0** (baseline 329) | **382** (baseline 480) |
| 0.03 | 0 | 125 |
| 0.1 | 0 | 2 |
| 1 | 0 | 0 |

There are no material false positives in this sweep. Comparing exact query identities finds **zero newly introduced material intersection or closest-distance failures**, and 1,552 resolved material false negatives, across the two variants. Closest-distance accuracy also remains outside the declared tolerance on smaller-scale inputs; counts and all exact query identities are retained in the CSV/JSON evidence, rather than hidden behind the intersection-only green suite. For example, sample 0 at scale 0.01, margin 0.001, axis W, sphere first still returns false for analytic overlap `1.00004253752e-5` with tolerance `9.99999997e-7`, and closest squared distance `3.27894333907e-8`.

The remaining failures are **not** excused as near-grazing within tolerance. Further lower-dimensional absolute-degeneracy/projection handling requires separate work; this bounded follow-up deliberately does not relax those limits or add an unproved blanket-hit fallback. **The SpaceSim collision-robustness acceptance gate remains FAIL/unsafe to certify.** Raw GJK no-progress is not a robust separation proof over the wider measured envelope.

Follow-up evidence is preserved under `/home/brkdrd/mover-build-collision/artifacts/small-scale-followup/`: trace source/overlay/dependency file and log; exact oracle JSON; RED/relative-only/GREEN logs; final Release and sanitizer build/full-suite/HelloWorld logs; baseline header overlay; envelope source, runner, binaries, CSVs and summary JSON; and benchmark source/scripts/results. The original diagnostic artifacts in the parent directory remain intact.

The follow-up small-scale mixed sphere–box ABBA benchmark uses 14 samples per variant with the same compiler/flags and 4,096 inputs. Median baseline/final GJK cost is 110.39/109.30 ns per call, with overlapping ranges 103.29–120.45 and 104.30–125.73 ns. This excludes an obvious order-of-magnitude regression in that workload; it is not a speedup or server-capacity claim.

A separate exploratory variant of the initial repair using double cofactors but already-rounded float dot products also passes the initial regression corpus. That does **not** validate that shortcut over a wider envelope; both retained double solves promote coordinates before multiplication/accumulation.

## Performance sanity, not a server-capacity claim

Same compiler/flags, original versus fixed header, AMD Ryzen 5 9600X, profiler and LTO off. An external benchmark uses 4,096 seeded inputs, non-inlined entry wrappers, checksum-consumed results, one million predicate calls and 200,000 sphere–box queries per sample. ABBA ordering yields 14 samples per variant:

| Median ns/call | Original | Fixed |
|---|---:|---:|
| Direct pentachoron predicate | 30.57 | 32.78 |
| Mixed random sphere–box GJK | 118.96 | 110.68 |

Predicate cost is about 1.07× in this probe. GJK timing ranges overlap substantially (original 104.32–148.70 ns; fixed 103.89–130.31 ns); concurrent host activity means the apparent GJK speedup is **not** a reliable optimization claim. This only rules out an obvious order-of-magnitude regression in the tested workload. It is not whole-scene, tail-latency, or 120 Hz acceptance evidence. Source, compiler commands and every measurement are preserved in `predicate-benchmark.cpp`, `collision-benchmark.py`, `benchmark.json` and raw logs.
