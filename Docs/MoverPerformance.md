# Mover performance: bounded space-dynamics throughput

Status: **measurement record**, 2026-11-05. This reports what one bounded workload
did on one machine. It is **not** a production capacity, frame-budget, reliability or
speedup claim. No production tolerance or population budget exists yet, and none is
established here.

Harness: [`Tools/SpaceSimBenchmark`](../Tools/SpaceSimBenchmark) (new, self-contained
CMake subproject; no production source or CMake file was modified). Raw per-run JSON,
per-step CSV and logs are in `Tools/SpaceSimBenchmark/results/`; a condensed copy of
the numbers below is in
[`Docs/benchmark-results/2026-11-05-spacesim-throughput.json`](benchmark-results/2026-11-05-spacesim-throughput.json).

## Harness review and final-source acceptance

The final harness passed independent review after repairs to exact time serialization, compile-command provenance validation and artifact destination aliasing. All 47 permanent regressions passed on Linux and NTFS through WSL; independent dangling-link, filesystem-error and valid-path probes also passed. Numerical workloads, reference calculations and the `1e-4` separated-scene gates are unchanged. No protection against concurrent adversarial filesystem mutations is claimed.

The historical measurements below are retained with their original provenance. Final-source capture is recorded in the next section; it supersedes older measurements for the current source without rewriting their provenance.

## Final-source capture (`e503e2db`)

The parent reconfigured and rebuilt the approved harness from clean commit `e503e2db3d653072e0e45e12e29f44f48ca73609` (zero dirty entries). All **8 runs / 9,600 measured steps** passed the unchanged gates; both repeats and the separated 256-body thread-count comparison have matching state hashes. CTest passed the 47-regression suite again. Raw JSON, CSV, logs and environment/compile-command metadata are archived under `Tools/SpaceSimBenchmark/results/final-e503e2db/`.

Binary SHA-256: `71069ca796123d9b5fb1965b8011eb547d16118c2904c83d55ef83465e6660ed`. Core source fingerprint: `ffdd3318db5698795ea1722e46eac16e028ea7cabed3eaecc48b6f338b8ca1ba`. GCC 13.3, float positions, no fast-math, library LTO on and cross-platform-deterministic option off; exact per-target commands are retained. This measurement configuration differs from the precise validation presets.

| Run | p50 ms | p95 ms | p99 ms | Max ms | Misses at 120 Hz |
|---|---:|---:|---:|---:|---:|
| interacting-128.rep1 | 0.147 | 0.223 | 0.447 | 0.533 | 0 |
| interacting-128.rep2 | 0.146 | 0.189 | 0.220 | 0.411 | 0 |
| separated-128.rep1 | 0.137 | 0.167 | 0.416 | 0.448 | 0 |
| separated-128.rep2 | 0.137 | 0.153 | 0.193 | 0.309 | 0 |
| separated-256-4t.rep1 | 0.456 | 0.513 | 0.738 | 0.879 | 0 |
| separated-256-4t.rep2 | 0.448 | 0.529 | 0.621 | 0.907 | 0 |
| separated-256.rep1 | 0.267 | 0.374 | 0.555 | 0.699 | 0 |
| separated-256.rep2 | 0.267 | 0.320 | 0.553 | 0.623 | 0 |

Zero deadline misses in these small workloads do not establish population capacity or server-level headroom. The interacting runs measured total energy drift `0.003217738`; that field is reported but not gated for contacts, as specified before the run. Separated-scene energy/spin gates remain `1e-4`.

### Retained dense-contact stress workload

The native precise/FP-traps `PerformanceTest` was also rebuilt from the same engine source and ran the existing Pyramid scene at an actual 120 Hz for 1,200 steps, with no warmup samples discarded. Both thread counts completed with identical final hash `0x1a5f81a7674216e6`. Raw timings and manifest are in `Docs/benchmark-results/pyramid-e503e2db/`; a preliminary six-thread run made during concurrent review is also retained separately.

| Total threads (including main) | Mean ms | p95 ms | p99 ms | Max ms | Misses at 120 Hz |
|---|---:|---:|---:|---:|---:|
| 6 | 8.210 | 13.595 | 20.463 | 163.386 | 164 / 1,200 |
| 12 | 6.661 | 10.468 | 15.273 | 85.183 | 108 / 1,200 |

**The dense-contact workload does not meet a hard 120 Hz deadline.** One run per configuration is not a reliability estimate; these are not matched upstream speedup measurements. The much lighter space-dynamics results must not replace or conceal these stress results.

## What was measured

`PhysicsSystem::Update` called in a loop on a real `PhysicsSystem` at a fixed
`dt = 1/120 s`, **1200 steps per run = 10.0 s simulated for every body count**. Every
step is timed with `std::chrono::steady_clock`, **including the first**; there is no
warmup phase, no discarded startup and no early exit. Sleeping is disabled both per
body and in `PhysicsSettings`, so no step is skipped by the sleep heuristic. That is a
deliberate choice for a throughput workload and it is recorded in the manifest.

**separated** (primary): 128 or 256 dynamic tesseracts on a 4D grid with 8 m spacing
against a 1.83 m bounding radius, so the measured minimum bounding-surface gap is 4.344 m and no
contact is possible. Half extents 0.40 / 0.70 / 1.00 / 1.30 (1 : 1.75 : 2.5 : 3.25
anisotropy, distinct principal moments), convex radius 0. Zero gravity, zero linear
and angular damping, common translational velocity (0.35, -0.20, 0.15, 0.25) m/s so
mutual separation is analytically invariant, six nonzero angular-velocity components per
body, `mApplyGyroscopicForce = true`.

**interacting** (secondary): 128 dynamic tesseracts in two 64-body lattice clusters at
x = ±16, closing at 4 m/s each (spacing 5 m, so no initial overlap), restitution 1,
friction 0, same spin and gyroscopic setting. Contacts occur in the mid-run collision
phase; the run recorded 376 `OnContactAdded`, 217 `OnContactPersisted` and 376
`OnContactRemoved` **callback invocations** (callback counts, not distinct body pairs).

## Manifest

| Field | Value |
|---|---|
| Hardware | AMD Ryzen 5 9600X, 6 cores / 12 threads, 32 MiB L3; WSL2 Ubuntu on Windows 11 |
| Build | GCC 13.3.0, CMake 3.28.3, Release `-O3`; **no `-ffast-math`**; `-ffp-contract=off` on the library TUs; LTO/IPO on for the library (`-flto=auto`); `JPH_USE_AVX2`/SSE4.2/FMA/BMI as selected by the Jolt target; profiler **off**, debug renderer **off**, asserts **off**, float positions. The harness TU is built separately: no `-ffp-contract=off`, no `-flto`, `-Wall -Wextra` (TU flags differ; see below) |
| Revision | `6a32e275b38c096f5ef886d381ad032314c1564e`, working tree dirty (14 entries, all other workstreams; no core physics edit in this task). The dirty-entry count is not an identity of the tree: **no source fingerprint exists for this capture** (see Provenance) |
| Binary | SHA-256 `d75bd07e5c07f5a9e94e69cfefa0a0095ab5535518ca88a23ef0b72fc916caed` |
| Population | 128 / 256 dynamic tesseracts, one shape instance shared, no static bodies |
| Time/precision | 120 Hz, one collision step per update, 1200 steps, float positions |
| Solver | `mNumVelocitySteps = 10`, `mNumPositionSteps = 2` (engine defaults), deterministic simulation on |
| Worker threads | 0 (single-threaded job system) for the primary rows; 4 for the scale check |
| Instrumentation | one `steady_clock` pair per step, measured at 16.4–18.1 ns; numeric sampling every 60 steps **between** timed calls (0.23–0.53 ms per run total), i.e. **20 sampled states** for a 1200-step run |

The two translation units are not compiled with the same flags, and the manifest records
both. The library gets `-ffp-contract=off`, `-flto=auto`, `-fno-rtti`, `-fno-exceptions`
and `-Wall -Werror` from the `Build/` directory scope. The harness TU gets
`-O3 -std=c++17 -Wall -Wextra`, the Jolt PUBLIC ISA options, and the compiler's default
FP contraction (GCC may contract; there is no `-ffast-math` anywhere). IPO is enabled for
the library; the harness TU and the link command add no `-flto` (GCC's linker plugin
still processes the library's LTO objects). Exact command lines are recorded in
`environment_*.txt` from `compile_commands.json`.

## Results — full-run latency distribution

`steps/s` is 1200 / (sum of all step times). Percentiles are nearest-rank over all
1200 steps. `setup` covers allocation, `PhysicsSystem::Init`, body creation and
`OptimizeBroadPhase`; it is excluded from the step distribution and reported separately.

| run | repeats | p50 ms | p95 ms | p99 ms | max ms | mean ms | steps/s | setup ms | wall ms | misses > 8.333 ms |
|---|---|---|---|---|---|---|---|---|---|---|
| separated 128 | 1 | 0.142 | 0.175 | 0.417 | 0.495 | 0.151 | 6622 | 0.122 | 181.5 | 0 |
| separated 128 | 2 | 0.142 | 0.173 | 0.303 | 0.455 | 0.150 | 6687 | 0.128 | 179.7 | 0 |
| separated 256 | 1 | 0.277 | 0.388 | 0.583 | 0.667 | 0.295 | 3385 | 0.181 | 355.1 | 0 |
| separated 256 | 2 | 0.274 | 0.311 | 0.552 | 0.594 | 0.284 | 3522 | 0.190 | 341.2 | 0 |
| interacting 128 | 1 | 0.157 | 0.204 | 0.241 | 0.643 | 0.164 | 6106 | 0.131 | 196.8 | 0 |
| interacting 128 | 2 | 0.158 | 0.210 | 0.395 | 0.583 | 0.168 | 5954 | 0.138 | 201.9 | 0 |
| separated 256, 4 threads | 1 | 0.471 | 0.540 | 0.793 | 0.935 | 0.447 | 2238 | 0.619 | 536.8 | 0 |
| separated 256, 4 threads | 2 | 0.469 | 0.548 | 0.752 | 1.017 | 0.445 | 2247 | 0.511 | 534.6 | 0 |

No step exceeded 6.667 ms (the 150 Hz stretch deadline) in any run either.

The first step is included and is visible in the data: 0.162 ms against a 0.142 ms
median for separated 128, and 0.292/0.322 ms against a 0.277/0.274 ms median for
separated 256. Nothing was excluded to make the distribution look better.

Per-body cost is essentially flat across the doubled population: 1.11 µs per body per
step at 128 bodies, 1.08 µs at 256 (medians), i.e. this workload scales linearly in
body count up to 256 on one core.

## Numerical checks (independent double-precision reference)

The harness recomputes each body's spin angular momentum `I(Ω) = Q diag(I) Qᵀ Ω` and
kinetic energy `½m|v|² + ½Ω·L` from public state with 4×4 skew matrices for the six angular components — no engine
inertia operator, gyroscopic integrator or energy helper is used. Angular momentum
about the origin is `I(Ω) + m(v ∧ x)`. The reference itself is checked by
`--selftest` against analytic fixtures (skew round-trip, `L = IΩ` and
`E = ½ΣIω²` at identity orientation, rotor → matrix orthogonality and determinant +1,
symmetry/positive-definiteness of the rotated inertia operator, energy invariance under
rigid reorientation, nearest-rank percentiles on 1..100). All self-tests pass.

Drift is reported relative to `Σ|E₀|`, `Σ|L₀|` and `Σm|v₀|` — well-conditioned scales;
normalising by a nearly cancelling total (the interacting clusters have almost zero net
momentum) would be numerically meaningless.

| run | energy drift (per body / total) | spin L drift (per body) | total L drift | total p drift | min pair gap |
|---|---|---|---|---|---|
| separated 128 | 5.92e-5 / 1.66e-5 | 2.62e-6 | 2.63e-5 | 0 | 4.344 m |
| separated 256 | 6.21e-5 / 1.66e-5 | 2.62e-6 | 2.38e-5 | 0 | 4.344 m |
| interacting 128 | not gated / 4.81e-3 | not gated | 1.88e-5 | 1.08e-8 | n/a |

All separated-workload gates pass. State stays finite, all rotations stay normalized,
every body stays active, no velocity clamp fires and no `PhysicsUpdateError` was
returned in any run.

These state checks are evaluated at the initial state and at the **20 sampled states**
(every 60 steps plus the final step), not at all 1200 steps. The clamp flags and
`MotionProperties::HadNonFiniteVelocity()` are sticky, so a clamp or a non-finite
velocity is still seen at the next sample; the original capture did not consult the
non-finite-velocity diagnostic at all, and that check was added afterwards (it appears
in the validated rerun below, not in the historical runs above).

The reference now forms the translational energy and the `Σm|v₀|` scale from double
components; the original capture used the engine's float `Vec4::LengthSq()`, which fed
a rounded float into the reference. That affects reported reference values only, at a
level far below the drift gate; the historical numbers above are left exactly as
captured.

Declared gates (in `src/main.cpp`, fixed before the runs, not tuned afterwards):
per-body energy and spin-L drift ≤ 1e-4 for free rotation, total L ≤ 1e-3, total p
≤ 1e-5, positive pair gap, contacts > 0 for the interacting workload, plus the
finiteness/activity/clamp/non-finite-velocity checks. They are prototype float32-arithmetic
bounds, not production tolerances. The measured 5.9–6.2e-5 per-body energy drift is
within the declared 1e-4 gate; the gate was not widened after seeing it.

The interacting workload's total energy changes by 4.81e-3 (0.48 %) across the
collision phase. It is **reported, not gated**: no accuracy budget exists for the
contact/position solver, and this run does not establish one.

## Reproducibility

Each configuration was run twice with the fixed seed; the two end-state hashes are
identical in every case. The 256-body separated run produces the **same** final state
hash (`85c794c0eeb3435f`) with the single-threaded job system and with 4 worker
threads, i.e. this workload is bit-identical across those thread counts on this build.

## Provenance and the validated rerun

The numbers above are the **original capture** of 2026-11-05. That binary is identified
by SHA-256 `d75bd07e5c07f5a9e94e69cfefa0a0095ab5535518ca88a23ef0b72fc916caed`, built
from revision `6a32e275` with a dirty working tree of 14 entries. The exact tree is no
longer reconstructible, so **no build manifest or source fingerprint is claimed for
it** — none was recorded at the time, and none is invented here. The capture remains
valid evidence for that core revision and its recorded numbers are not rewritten.

The core has since changed (`aee62ba5`, rank-aware mass creation; the principal
frame/eigenvalue ordering differs, which can change the gyroscopic splitting and the
numerical trajectory of these fixtures). A rerun was therefore built from the current
tree into a **fresh, separate build directory** and written to a separate output
directory, [`Tools/SpaceSimBenchmark/results/verified-current/`](../Tools/SpaceSimBenchmark/results/verified-current/),
with its own fingerprints:

| field (verified rerun) | value |
|---|---|
| revision | `aee62ba525d089ff45a993fa507de9ee32e2d680`; 21 dirty entries recorded at configure time and at run time (other workstreams were editing the tree; the count is not an identity, the fingerprints are) |
| binary | SHA-256 `624df6ee240532c6041136fc26c43cc589bf9d785f02ed2665686d7203b05991` |
| harness source fingerprint | `4d4ca3c35601dc92623ab30fcac294b4e6b965aeb12902fcae39818d0979563b` |
| core source fingerprint | `8a2d06bc01286f1b98f210738854b536bee50348fb08c4f7bdda7b2bc58a4c87` |
| build fingerprint | `e186d3da38175222782abaae44bae9b134487c1aa2589f65b4b1c66003ef42a3` |
| build flags | GNU 13.3.0 Release `-O3`; no `-ffast-math`; library `-ffp-contract=off` and `-flto=auto`; harness TU neither |

Results (repeat 1 of each configuration; both repeats, per-step CSVs and the exact
compile commands are in that directory):

| run | p50 ms | p95 ms | p99 ms | max ms | steps/s | energy drift (per body / total) | total L | total p | state hash |
|---|---|---|---|---|---|---|---|---|---|
| separated 128 | 0.141 | 0.158 | 0.181 | 0.278 | 6915.6 | 1.58e-5 / 2.98e-6 | 2.63e-5 | 0 | `4e94ad1aa6486f57` |
| separated 256 | 0.275 | 0.316 | 0.366 | 0.401 | 3529.5 | 1.62e-5 / 3.03e-6 | 2.38e-5 | 0 | `4eae817888d2cf89` |
| interacting 128 | 0.147 | 0.186 | 0.219 | 0.310 | 6455.2 | not gated / 3.22e-3 | 6.01e-7 | 1.36e-8 | `f86af8530ee92cba` |
| separated 256, 4 threads | 0.460 | 0.509 | 0.568 | 0.720 | 2292.1 | 1.62e-5 / 3.03e-6 | 2.38e-5 | 0 | `4eae817888d2cf89` |

All gated checks pass, no step exceeded 8.333 ms, no update error occurred, and the
non-finite-velocity diagnostic (new in this harness revision) reports no event. The
rerun reproduced its own repeats and the 0-vs-4-thread hash (`4eae817888d2cf89`) for
the 256-body separated workload. Contact callbacks in this trajectory: 342 added, 201
persisted, 341 removed.

The working tree kept moving during this work (`Jolt/SpaceSim/Simulation.h` and
`Jolt/SpaceSim/PhysicsBridge.h` were edited at 18:40Z and 18:44Z, after the library was
built at 18:39Z). The fingerprints above identify the **configure-time content the
binary was built from**, not the later tree; they are not a claim about the frozen
sources. The parent's subsequent source-freeze rerun is recorded in the final-source section above.

The state hashes and drifts differ from the historical capture (`85c794c0eeb3435f` /
`01b2945a285796cd`, per-body energy drift 5.9–6.2e-5): this is a **different core
trajectory**, not a reproduction, and the difference is not presented as a regression
or an improvement. The 1e-4 gate was not widened or otherwise changed; the rerun's
1.6e-5 per-body energy drift passes it as measured. This rerun is a snapshot taken
before the parent source freeze, not a final number.

## Observations (bounded by the data above)

- On one core, 1200 steps of this workload take 181 ms (128 bodies) and 355 ms
  (256 bodies): p50 step times 0.142 ms and 0.277 ms against an 8.333 ms step budget,
  with the first step included. This is a throughput headroom observation for **this**
  workload, not a statement that a game server meets 120 Hz.
- Running the same 256-body workload on 4 worker threads is **slower** (p50 0.471 ms
  vs 0.277 ms). With zero contacts the per-update parallel work is too small to repay
  barrier and scheduling cost on this machine. Reported as measured; not generalised
  to contact-heavy scenes.
- The interacting workload's tail is higher in the collision phase (p99 0.24–0.40 ms
  against p50 0.16 ms); the largest single step (0.643 ms) falls at step 524 of 1200,
  i.e. during the collision, not at startup.
- WSL2 scheduling noise is visible in the separated tails (p99 up to 0.42 ms with a
  0.142 ms median). These are full-run numbers including that noise, not a quiet-host
  best case.

## Not established here

- **No baseline comparison.** No before/after upstream 4D Jolt comparison of this workload was produced, so no
  speedup ratio is claimed — only absolute latencies for this configuration.
- No claim beyond 256 bodies, one process, one machine, one seed, one build.
- No production frame-budget, headroom or acceptable-miss policy. Zero misses in these
  short runs is not a reliability certificate; the specification still requires a
  declared envelope and miss policy.
- The contact solver's energy behaviour is measured, not certified, and no accuracy
  budget is implied.

## Reproduce

```bash
# from the repository root, in WSL/Linux (use a fresh build and output directory)
cmake -S Tools/SpaceSimBenchmark -B /home/<user>/mover-perf-build-verified2 -DCMAKE_BUILD_TYPE=Release
cmake --build /home/<user>/mover-perf-build-verified2 --target SpaceSimBenchmark -j2
Tools/SpaceSimBenchmark/run_benchmark.sh /home/<user>/mover-perf-build-verified2/SpaceSimBenchmark \
    Tools/SpaceSimBenchmark/results/verified-current-2
```

The script records the environment (CPU, compiler, `git rev-parse HEAD`, dirty-entry
count **and dirty-path list**, binary SHA-256, the harness source fingerprint, and the
actual library/benchmark compile commands from `compile_commands.json`) into
`environment_*.txt`. It runs each configuration twice, requires 16-hex-digit
end-state hashes and identical repeats, requires the separated-256 hash to agree
between 0 and 4 worker threads, parses every JSON strictly (rejecting `nan`/`inf`),
checks each manifest against the requested configuration and against the checked-out
harness sources, requires one CSV row per step, runs the harness self-test, and exits
non-zero if any check fails. It writes only into the output directory it is given,
refuses to write into a directory that already contains per-run JSON (unless
`SPACESIM_BENCH_ALLOW_OVERWRITE=1`), and does not modify the original `results/`
capture. Nothing was committed to git for this task.

The CLI/artifact regressions can be re-run on any built binary:

```bash
python3 Tools/SpaceSimBenchmark/tests/test_cli_and_artifacts.py /home/<user>/mover-perf-build-verified2/SpaceSimBenchmark
```
