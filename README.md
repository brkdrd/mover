# 4D Jolt

**Rigid-body physics in four spatial dimensions, built from Jolt Physics.**

4D Jolt extends [Jolt Physics](https://github.com/jrouwe/JoltPhysics) to a world with four spatial coordinates: **X, Y, Z and W**. The fourth coordinate is another direction in space, not time. The project explores collision detection, contact response and rotational dynamics beyond three dimensions.

The intended application is a **4D multiplayer space combat simulator**: physically simulated ships, distributed thrusters, damage, projectiles, debris and asteroid environments. The repository currently provides an experimental headless physics foundation. Ship systems, large-world management and multiplayer are proposed work, not implemented game features.

[Build and run](#build-and-run) · [Space-sim specification](Docs/SpaceSimSpecification.md) · [System architecture](Docs/SpaceSimArchitecture.md) · [Performance investigation](Docs/PerformanceOptimization.md)

## Project status

The headless v1 milestone builds and runs a 4D simulation. This is an independent, incomplete port of **Jolt 5.5.1**, not an official Jolt release or a drop-in replacement for the 3D library.

| Area | Current scope |
|---|---|
| Mathematics | Four-component spatial vectors, six-component angular bivectors, Spin(4) rotors and 4D rigid transforms |
| Geometry | 4D GJK/EPA, convex hulls, bounds, clipping and ray-query foundations |
| Shapes | Hyperspheres, tesseract boxes, capsules, convex hulls, hyperplanes, compound and decorated shapes; some operations remain incomplete |
| Simulation | Rigid bodies, broad phase, contact manifolds, normal response and three tangent friction directions |
| Examples and checks | Headless `HelloWorld`, unit tests and `Pyramid` / `MaxBodies` performance scenes |
| Space-sim layer | Design specification only; no production ship controller, weapon system or networking stack yet |

### Known limits

- **Collision robustness is unfinished.** A reproducible sphere-box GJK false negative exceeds the existing accuracy bound. It occurs in both the v1 baseline and the optimized code. See the [investigation and exact-input record](Docs/PerformanceOptimization.md#validation-and-an-existing-accuracy-gap).
- **Rendering is deferred.** Build with the debug renderer disabled. The inherited 3D Samples and JoltViewer are not a supported 4D viewer.
- General joint constraints, character controllers, upstream vehicle systems, soft bodies, ragdolls, hair, and mesh/height-field collision are excluded from the current v1 build. Source files remaining in the tree do not imply 4D support.
- Large-world precision, all upstream platforms and configurations, and cross-platform determinism have not been comprehensively validated for this port.
- A reliable full-workload 120 Hz deadline has **not** been established. Higher average throughput in some measurements does not remove collision-burst latency.

Recorded validation on 2026-11-04: WSL/GCC Release and assertions-enabled ASAN/UBSAN passed **389 cases / 2,397,320 assertions**. Native MSVC passed **389 of 390 cases**, with the known GJK failure above. These counts describe the enabled 4D suite, not all upstream Jolt features. Full evidence and conditions are in the [performance report](Docs/PerformanceOptimization.md) and [machine-readable results](Docs/PerformanceOptimizationResults.json).

## What changes in four dimensions?

A rigid body has four translational and six rotational degrees of freedom. Rotation takes place in planes: XY, XZ, XW, YZ, YW and ZW. The port uses `Rotor` for orientation and `Bivec` for angular quantities; a 3D quaternion or three-component angular velocity is insufficient.

Collision geometry changes too: a full GJK simplex can contain five vertices, convex boundary cells are three-dimensional, and a contact has three independent tangent directions. This requires changes to geometry and contact solving as well as vector types.

The repository's torque convention is **`Bivec::sWedge(force, offset)`**, consistent with its rotor integration and point-velocity contraction. See the [coordinate and dynamics contract](Docs/SpaceSimArchitecture.md#coordinate-and-dynamics-contract) before integrating custom forces.

## Build and run

Use CMake 3.20 or newer and a C++17 compiler. The currently documented validation environments are x86-64 Linux/WSL with GCC 13.3 and native Windows with MSVC 19.50. The default x86 options enable AVX2; other CPU configurations require their own validation.

### Linux / WSL

From the repository root:

```sh
cmake -S Build -B out/headless -DCMAKE_BUILD_TYPE=Release \
  -DDEBUG_RENDERER_IN_DEBUG_AND_RELEASE=OFF \
  -DPROFILER_IN_DEBUG_AND_RELEASE=OFF \
  -DTARGET_SAMPLES=OFF -DTARGET_VIEWER=OFF \
  -DTARGET_UNIT_TESTS=ON -DTARGET_HELLO_WORLD=ON \
  -DTARGET_PERFORMANCE_TEST=ON
cmake --build out/headless --target UnitTests HelloWorld PerformanceTest --parallel 6
./out/headless/UnitTests
./out/headless/HelloWorld
```

`HelloWorld` is a headless simulation and smoke test; a nonzero exit code indicates a failed check. Start an application integration with [HelloWorld.cpp](HelloWorld/HelloWorld.cpp).

### Windows

Use a Visual Studio developer shell and the [tested native build recipe](Docs/PerformanceOptimization.md#reproducing-the-native-build). It documents the compiler options used for the existing measurements and the known test failure. The inherited [build manual](Build/README.md) remains useful for general options, but its platform and feature lists describe upstream Jolt and are not a 4D support guarantee.

### Performance tests

```sh
./out/headless/PerformanceTest -s=Pyramid -q=Discrete -t=6 -i=600 -hz=60
```

The 4D pyramid contains 1,296 dynamic boxes with sleeping disabled. It stresses sustained contact solving; it is not a ship-combat benchmark. `-hz` sets the simulated timestep, while reported steps/second measures execution throughput. Do not equate these quantities.

See [benchmark options and the comparison runner](Docs/PerformanceTest.md). The `MaxBodies` stress scene has very large memory requirements; inspect its setup before running it.

## Direction: a 4D space combat simulator

The proposed design combines detailed rigid-body ship dynamics with scalable collision and world management. It prioritizes:

1. Conservation laws and stable numerical integration.
2. Correct 4D rotational mechanics and changing mass properties.
3. Component-based forces, propulsion and flight control.
4. Reliable collisions at combat speeds.
5. Precision across large distances and velocities.
6. Large object populations, including asteroid belts and debris.
7. Predictable server cost and reproducible PvP outcomes.

The [specification](Docs/SpaceSimSpecification.md) defines these requirements, validation gates and development stages. The [architecture](Docs/SpaceSimArchitecture.md) explains state ownership, scheduling, collision response, damage, reference frames and networking. Both distinguish proposals from existing capabilities.

## Repository guide

| Path | Purpose |
|---|---|
| [Jolt/](Jolt/) | Physics library and 4D mathematics |
| [UnitTests/](UnitTests/) | Enabled 4D checks and retained upstream tests awaiting porting |
| [HelloWorld/](HelloWorld/) | Minimal working 4D simulation |
| [PerformanceTest/](PerformanceTest/) | Headless benchmarks and repeatable comparison runner |
| [Docs/](Docs/) | Design, measurement reports and inherited documentation |
| [Build/](Build/) | CMake configuration and platform build scripts |
| [CLAUDE.md](CLAUDE.md) | Detailed porting history, mathematical conventions and development notes |

The `4d` branch contains the port. `master` preserves the upstream 3D reference. Upstream documentation, release notes, integrations and bindings retained in this tree should not be assumed compatible with the 4D API.

## Contributing

Open an [issue](https://github.com/brkdrd/4dJolt/issues) with a reproducible problem or discuss a substantial change before implementing it. Target contributions at `4d`. Include relevant correctness tests and disclose build settings when reporting performance. Do not weaken accuracy bounds or reduce physical workload to claim an optimization.

Useful contributions include robust 4D collision algorithms, independently checked rotational mathematics, representative space-combat benchmarks and visualization. The historical upstream contributor agreement and automation are retained; they do not establish a separate CLA process for this fork.

## Credits and license

Based on **Jolt Physics by Jorrit Rouwe and contributors**. The upstream architecture and implementation make this experiment possible. The 4D adaptation is developed in [brkdrd/4dJolt](https://github.com/brkdrd/4dJolt).

The library is distributed under the [MIT license](LICENSE). Preserve existing copyright and license notices. Some examples and assets carry their own notices; see their files and [Assets/LICENSE](Assets/LICENSE).
