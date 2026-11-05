# SpaceSimBenchmark — bounded space-dynamics throughput harness

A small, reproducible throughput and numerical-integrity measurement of the real
`PhysicsSystem` loop on the 4D engine. It is **measurement code**: no production code
is changed and no production API is added.

## Artifact destination safety

Artifact destinations are validated before workload execution. Colliding JSON/CSV paths, dangling symlinks (including chains and parent components), symlink loops and filesystem inspection failures reject with exit 2 without creating or modifying artifacts. Resolving symlinks and distinct ordinary files remain supported. Independent review passed the 47-regression suite on Linux and NTFS through WSL, plus 35 additional probes. This is CLI correctness validation, not protection against concurrent adversarial filesystem mutations or a claim of native Windows/MSVC benchmark coverage.

## What it measures

`PhysicsSystem::Update` at a fixed 120 Hz timestep over a fixed number of steps (1200
by default = 10.0 s of simulated time, identical for every body count). Each step is
timed individually with `std::chrono::steady_clock`, including the first one. There is
**no warmup, no discarded startup and no early exit**. Reported:

* per-step min / p50 / p95 / p99 / max / mean (nearest-rank percentiles),
* total wall time, sum of step times, and setup time separately,
* deadline misses against 8.333 ms (120 Hz) and 6.667 ms (150 Hz),
* update errors, non-finite state, non-normalized rotation, velocity clamps,
* energy and angular-momentum drift against an independent double-precision reference,
* a state hash so repeated runs can be checked for bit-identical determinism.

Timing overhead of the instrument is measured and printed (a `steady_clock` pair costs
tens of ns on the target machine, against a ~ms step). Numerical sampling happens
**between** timed calls and its cost is reported separately; it never inflates a step
time. Sleeping is disabled by construction (both the body setting and the system
setting) so that no step is skipped by the sleep heuristic; this is a throughput
workload, not a resting-scene workload.

### What the numerical checks actually cover

State checks run on the initial state and on every sampled state (every 60 steps plus
the last step: 20 in-run samples for the 1200-step runs), **not on all 1200 steps**.
The velocity-clamp flags and `MotionProperties::HadNonFiniteVelocity()` are sticky, so
a clamp or non-finite velocity during a step is still observed at the next sample; a
repaired transient between samples is only caught by those core diagnostics, not by
the state scan.

## Workloads

* `--scenario separated` — N dynamic **anisotropic** tesseracts (half extents
  0.40/0.70/1.00/1.30, no convex radius) on a 4D grid with spacing 8 m, which is more
  than twice the bounding radius, so no two bodies can ever touch. Zero gravity, zero
  damping, a **common** translational velocity (so mutual separation is analytically
  invariant), and a six-component angular velocity per body with the gyroscopic term
  enabled (`BodyCreationSettings::mApplyGyroscopicForce`). The six principal plane moments
  differ, so the gyroscopic term is actually exercised. Each body's spin angular
  momentum and kinetic energy are constant analytically; both are checked.
* `--scenario interacting` — N dynamic tesseracts in two lattice clusters closing on
  each other in zero gravity, restitution 1, friction 0. Contact callbacks are counted
  so the run proves it really interacted. Total linear and angular momentum are checked
  against the declared prototype bounds (contact impulses are internal); total energy
  and individual spin drift are not gated. These bounds are not a production
  contact-solver accuracy budget.

## Independent numerical reference

The harness recomputes, in double precision, from public body state only:

* spin angular momentum `L = Q diag(I) Q^T (Omega)` via 6 skew matrices (plane
  ordering `e12,e13,e14,e23,e24,e34`) — no engine inverse-inertia operator,
* rotational energy `0.5 * Omega . L`, plus translational `0.5 m |v|^2`,
* angular momentum about the origin `I(Omega) + m (v ^ x)`,
* total linear momentum `sum m v`.

The translational terms are formed from the double components, not from the engine's
float `Vec4::LengthSq()`; otherwise a rounded float input would leak into the
"independent" reference.

The reference's own correctness is checked by `--selftest` against hand-computable
fixtures (skew round-trip, identity-orientation `L = I Omega` and `E = 0.5 sum I omega^2`,
rotor→matrix orthogonality and determinant +1, symmetry and positive definiteness of
the rotated inertia operator, rotation invariance of energy, all four spatial
components in the translational energy, and nearest-rank percentiles on 1..100). A
failing self-test fails the run script.

## Command line

```
SpaceSimBenchmark --scenario separated|interacting --bodies N [--steps 1200] [--threads 0]
                  [--rate 120] [--seed 0x...] [--json out.json] [--csv steps.csv]
                  [--selftest] [--help]
```

`--bodies` must be an integer in [2, 4096], `--steps` in [1, 10000000], `--threads` in
[0, 64], and `--seed` an unsigned integer (base 0). Every numeric token must be fully
consumed: `128x`, `-1`, `1e3` for an integer option and `120junk` for `--rate` are
rejected rather than silently truncated.

`--rate` must be a **finite positive** number, and the resulting float timestep and the
simulated duration (`steps * dt`) must both be positive and finite. `--rate inf`,
`nan`, `-120`, `0`, `1e300` (dt underflows to 0) and `1e-300` (dt overflows) are
rejected before any physics runs: they do not describe a measurement.

Exit codes:

| code | meaning |
|---|---|
| 0 | every gated check passed **and** every requested artifact was written completely |
| 1 | failed check, non-finite reported value, or artifact I/O failure (open/write/flush/close) |
| 2 | invalid command line |

A JSON/CSV that cannot be written (missing parent directory, path is a directory,
`/dev/full`, short write, failed flush or close) makes the process exit 1 and print
`RESULT: FAIL`; a failed artifact is never reported as a successful run. JSON is never
emitted with `nan`/`inf` tokens: a non-finite reported value fails the run and
suppresses the artifact.

## Build

Portable CMake subproject. It uses the repository's own `Build/` project to compile the
Jolt library with the engine's exact compiler and ISA definitions (AVX2/SSE4.2/FMA/BMI
come from the `Jolt` target's PUBLIC interface — these are not hand-picked flags), and
links one executable. The measurement configuration is pinned: Release, profiler and
debug renderer off, asserts off, float precision.

```bash
# from the repository root, in WSL/Linux
cmake -S Tools/SpaceSimBenchmark -B /home/<user>/mover-perf-build -DCMAKE_BUILD_TYPE=Release
cmake --build /home/<user>/mover-perf-build --target SpaceSimBenchmark -j2
```

The build writes `compile_commands.json` and a generated
`SpaceSimBenchmarkBuildInfo.h` that is compiled into the binary. The library and the
benchmark translation unit do **not** share flags: the library gets
`-ffp-contract=off -flto=auto -fno-rtti -fno-exceptions -Wall -Werror` from the
`Build/` directory scope, while the harness TU gets `-O3 -std=c++17 -Wall -Wextra` with
the Jolt PUBLIC ISA options and the compiler's default FP contraction. Neither uses
`-ffast-math`. The per-run JSON records both flag strings, the definitions, the
relevant cache values (`INTERPROCEDURAL_OPTIMIZATION=ON`,
`CROSS_PLATFORM_DETERMINISTIC=OFF`, ...), and the fingerprints below; the exact command
lines are in `environment_*.txt` from `compile_commands.json`.

## Provenance

Each JSON manifest carries:

* `harness_source_sha256` — content hash of `src/main.cpp` + `CMakeLists.txt`;
  recomputed by the run script and required to match, so a stale binary is rejected,
* `core_source_sha256` — content hash over every existing source file of the `Jolt`
  target (the library actually linked),
* `build_sha256` — harness sources + core hash + definitions + both flag strings +
  cache values,
* `git_head_at_configure` and `git_dirty_entries_at_configure` — configure-time git
  state. The dirty **count is not an identity** of the tree: the fingerprints are.
  The run script also records the full dirty-path list at run time; it can differ from
  the configure-time count when other workstreams are editing the tree.

## Run

```bash
Tools/SpaceSimBenchmark/run_benchmark.sh /home/<user>/mover-perf-build/SpaceSimBenchmark Tools/SpaceSimBenchmark/results
```

The script runs each configuration twice with the same seed, requires 16-hex-digit
end-state hashes and identical repeats, verifies that `separated-256` has the same hash
with 0 and 4 worker threads, runs `--selftest`, writes per-step CSVs, per-run JSON,
logs, an environment manifest (CPU, compiler, git HEAD, dirty scope, binary SHA-256,
actual compile commands) and `summary.json`. Every JSON is parsed strictly (nan/inf
tokens rejected), its manifest is checked against the requested configuration, its
harness fingerprint against the checked-out sources, and every CSV must have exactly
one row per step. It exits non-zero if any check fails; failures are reported, never
suppressed. Use a fresh output directory for a rerun; the script refuses to write into a
directory that already contains per-run JSON unless `SPACESIM_BENCH_ALLOW_OVERWRITE=1`
is set, so the original capture cannot be overwritten by accident.

`results/` holds the **original capture** (2026-11-05, core revision `6a32e275`,
binary SHA-256 `d75bd07e...`) taken before build fingerprints existed.
`results/verified-current/` holds a separate, validated rerun on a later core
revision; see `Docs/MoverPerformance.md` for both and for why they are not
interchangeable.

## Regression tests

`tests/test_cli_and_artifacts.py` drives the real binary as a subprocess and covers the
CLI/artifact contract above: invalid rates and integers, missing artifact directories,
directory paths, `/dev/full` write failures, strict-JSON content, finite numbers,
manifest provenance, CSV row counts, repeat determinism and 0-vs-4-thread hash
equality.

```bash
python3 Tools/SpaceSimBenchmark/tests/test_cli_and_artifacts.py <path-to-SpaceSimBenchmark>
ctest --test-dir <build-dir> -R SpaceSimBenchmark
```

## What this does not establish

No production capacity, frame-budget or 120/150 Hz reliability claim. No population
budget. The acceptance gates in `src/main.cpp` are prototype float32-arithmetic gates,
declared in the manifest, not production tolerances. The interacting workload's
contact-solver energy behaviour is measured and reported, not certified.
