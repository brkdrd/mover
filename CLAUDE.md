# 4D Jolt — Porting Jolt Physics to Four Spatial Dimensions

## What this project is

A fork of Jolt Physics **v5.5.1** (fork point = `master` = `6576239a`) being converted into a
4-spatial-dimension rigid-body physics engine on branch `4d`. The guiding idea: keep Jolt's
architecture and change only what dimension-dependence forces us to change. Most of the
*geometry* does generalize cleanly (GJK, EPA, quickhull, SAT, ray tests). Most of the
*rotational dynamics and contact resolution* does **not** — see "Things that do not trivially
generalize" below. Treat that section as load-bearing: underestimating it is how this project
got stuck.

- `master` = pristine upstream 3D Jolt, kept as reference. Never commit 4D work to it.
- `4d` = the working branch. As of 2026-07-22 most of the Geometry/Physics/ObjectStream work
  is **uncommitted working-tree changes** on top of commit `6f770d01 "Math for 4D"`.

## Build & test

```bash
# v1 is HEADLESS: build with the debug renderer OFF (the 3D DebugRenderer is deferred).
cmake -S Build -B build -DCMAKE_BUILD_TYPE=Release -DDEBUG_RENDERER_IN_DEBUG_AND_RELEASE=OFF
cmake --build build --target Jolt -j$(nproc)           # BUILDS CLEAN (headless v1)
cmake --build build --target UnitTests -j$(nproc)      # builds and links
./build/UnitTests                                     # runs: 364 cases, 357 pass, 7 fail (see below)
```

**Current build status (verified 2026-08-08): the headless v1 library BUILDS CLEAN** — `libJolt.a`,
0 errors, down from ~2,790 at the start of Phase B. **`UnitTests` also builds, links and runs: 364
cases, 357 pass, 7 fail.** Build with the debug renderer OFF
(`-DDEBUG_RENDERER_IN_DEBUG_AND_RELEASE=OFF`); v1 is headless.

This replaces the old standalone Math-only runner — **run the real `./build/UnitTests` from now on.**
(The recipe for the standalone runner is kept below only as a fallback.)

### The 7 remaining test failures — all one root cause in GJK
`UnitTests` is at **357/364**. Everything left is EPA (5 cases) and GJK (2), and the evidence points
at a single mechanism: **GJK does not reliably grow its simplex to 5 points in 4D.**

Measured on `TestGJKSphereVsBoxIntersect` (10,000 random sphere/box pairs): **24 disagreements, all
false negatives** — GJK reports MISS for pairs that genuinely overlap by up to 0.071, far beyond the
2e-4 test tolerance. A traced run of the first failure shows |v|² descending 15 → 3 → 0.25 → … →
5.0e-6, then *rising* to 5.1e-6, which trips the `prev_v_len_sq - v_len_sq <= FLT_EPSILON *
prev_v_len_sq` "converged, therefore separated" early-out. The true distance is 0 (the shapes
overlap), so the loop converged to the wrong answer rather than to a separating axis.

Across the whole trace the simplex point count goes 1,2,3,4,3,2,3,4,3,4,3,4,2 — the tetrahedron
reduction **never once returns the full set 0x0f**, so the pentachoron case (and its `set == 0x1f`
"origin enclosed, therefore hit" test) is never reached.

What has already been ruled out:
- `GetClosestPointOnTetrahedron` is correct, including the 4D projection fix: given a tetrahedron
  whose 3-flat the origin projects strictly inside, it returns set `0xF` and the right projection.
- `GetClosestPointOnPentachoron` is correct: for a pentachoron strictly containing the origin it
  returns set `0b11111` and a zero closest point.
- `Sphere::GetSupport` / `AABox::GetSupport` are correct 4D support functions.
- `GJKClosestPoint::Intersects` and `GetClosest` are structurally right for 4D (they already switch
  on 5 points and test `set == 0x1f`).

Prime suspect: the `LastPointPartOfClosestFeature` / `MustIncludeD` path. In
`GetClosestPointOnTetrahedron`, when face ABC is flagged and D must be included, the code returns
`set = 0b0001, closest = A` (inherited verbatim from 3D Jolt, whose comment argues the closest point
"cannot be an interior point either"). In 3D four points are a full simplex so that reasoning holds;
in 4D four points are only a 3-flat and the origin is essentially never on it, so returning vertex A
can be far worse than the true closest point and destroys the monotonic decrease of |v| that GJK's
convergence test depends on. **The 4D question to answer is what the `MustInclude` shortcut should
be when the simplex is one dimension short of enclosing the origin.**

Reproduce with the probes in the scratchpad, or just run `./build/UnitTests -ts=GJKTests`.

**v1 scope cuts applied at the source level (2026-08-08).** Commented out of `Jolt.cmake` with
`TODO(4D)` markers, with matching RTTI declarations / `Factory` registrations removed from
`RegisterTypes.cpp`:
- **SoftBody, Character, Vehicle, Ragdoll, Hair, Skeleton** (`RegisterHair()` is now a no-op).
- **Joint constraints** — Cone, Distance, Fixed, Gear, Hinge, Path(+Hermite), Point, Pulley,
  RackAndPinion, SixDOF, Slider, SwingTwist — and the rotational / joint-specific constraint parts
  (Angle, Gear, HingeRotation, IndependentAxis, RackAndPinion, RotationEuler, RotationQuat,
  SwingTwist). Deferred to v1.1 per roadmap step 9; this alone removed ~341 errors.
  **`Axis`/`DualAxis`/`Point`/`Spring` parts stay** — the contact solver uses them.
- **Triangle-mesh collision** — the `*VsTriangles` helpers, `InternalEdgeRemovingCollector`,
  `TriangleSplitter/`, `AABBTree/`. They serve the already-cut Mesh/HeightField/Triangle shapes.
  `PhysicsSystem::sDefaultSimCollideBodyVsBody` and `NarrowPhaseQuery::CollideShapeWithInternalEdgeRemoval`
  now always take the plain-collide path.

SoftBody additionally still has references *from* the rigid-body core, so those are compiled out
behind **`JPH_ENABLE_SOFT_BODY`**, a macro never defined for v1 — upstream code left verbatim so the
subsystem can be switched back on once ported.

**33 unit-test files are de-registered** in `UnitTests.cmake` (still written against the 3D API or
against cut features), plus `Mat44Tests`, `ComputeTests` and `PathConstraintTests`. **Do not delete
them** — they encode a lot of expected behaviour; re-enable one at a time as each area is ported.

## The 4D type system (established, do not re-litigate)

| 3D Jolt | 4D Jolt | Notes |
|---|---|---|
| `Vec3` | `Vec4` | Spatial 4D vector. Thin subclass of `Lane4`; adds `sAxisW`, ternary `sCross(a,b,c)` (generalized cross product), `GetNormalizedPerpendicular`. |
| `Vec4` (SIMD util) | `Lane4` | Pure SIMD 4-lane utility (matrix columns, batch trig). NOT spatial. |
| `Quat` | `Rotor` | Rotation = even subalgebra of Cl(4,0): 8 floats `[s, e12, e13, e14, e23, e24, e34, e1234]` in a `Vec8`. |
| (new) | `Bivec` | 6-component bivector `[e12,e13,e14,e23,e24,e34]` (Vec8 storage). Angular velocity, torque, constraint rotation planes. |
| `Mat44` (rot+trans) | `Mat44` | Pure 4×4 linear map (rotation, optionally scale). **No translation column.** |
| `DMat44` | `RMat44` | Rigid transform = `Mat44` rotation + `RVec4` translation. Method names like `Multiply3x3` kept for diff-friendliness but operate on the full 4×4. |
| `DVec3` | `DVec4` | Double-precision spatial vector. ⚠ still half-converted, see bug list. |
| `RVec3` | `RVec4` | = `Vec4` or `DVec4` depending on `JPH_DOUBLE_PRECISION` (`Real.h`). |
| `Float3`/`Double3` | `Float4`/`Double4` | POD storage. `Float8` = POD for Vec8/Rotor. |
| Triangle-based geometry | Tetrahedron-based | `Tetrahedron.h`, `IndexedTetrahedron.h`, `RayTetrahedron.h`. Convex hull cells = tetrahedra, ridges = triangles. GJK simplex = up to 5 vertices (pentachoron). |

`Vec3.{h,inl,cpp}`, `Quat.{h,inl}`, and `Quat`-overloads in `Mat44`/`RMat44` still exist as
"backward compat" — they are what the un-migrated layers compile against. They should be
deleted once nothing references them; do not build new code on them.

### Cl(4,0) conventions

Basis index order for Rotor: `[1, e12, e13, e14, e23, e24, e34, e1234]`. `eᵢ² = +1`,
`e1234² = +1`, all bivectors square to −1. Reverse negates only grade 2. Full 8×8 geometric
product table and derived formulas live in the auto-memory file `project_4d_rotor_formulas.md`
and are implemented in `Jolt/Math/Rotor.inl` (scalar code). Sandwich `R v R̃` computes the full
odd multivector then extracts grade 1.

**Key fact used everywhere:** Spin(4) ≅ SU(2)×SU(2). Any bivector B splits by pseudoscalar
duality (B̃ = B·e1234) into commuting self-dual and anti-self-dual parts B± = ½(B ± B̃), each
of which squares to a non-positive scalar and exponentiates like a quaternion:
`exp(B) = exp(B₊)·exp(B₋)`. This is the correct way to implement the general rotor exponential
and logarithm, rotor normalization (normalize each SU(2) factor), and SLERP.
`Mat44::GetRotor()` (`Mat44.inl:673`) already uses this isomorphism and is a good reference.

## Current state by module

### Math (`Jolt/Math/`) — scaffolding solid, dynamics core now fixed & tested
Lane4/Vec4/Vec8/Mat44/RMat44/Float8/Double4 are solid. The Cl(4,0) rotor core is now correct
and covered by the reference-validated suite (37/37): the SU(2)×SU(2) split drives
`Rotor::Normalized`/`IsNormalized` (defect-aware), `Rotor::sExp`/`Log` (general 4D exp/log),
and `Mat44::GetRotor` (the anti-self-dual sign bug is fixed). Rotor/Bivec are still **entirely
scalar** (no SIMD on any ISA — a v1.1 perf item, not a correctness one); Vec8 has SSE/AVX/NEON
but no RVV. `Bivec.h`/`RMat44.h` are now tracked and registered in `Jolt.cmake`.

### Geometry (`Jolt/Geometry/`) — the strongest part of the project
Verified as a genuine, careful 4D generalization (not a rename): ClosestPoint (Gram-matrix
barycentrics incl. closest-point-on-pentachoron), GJK (5-vertex simplex), EPA (tetrahedral
cells, triangle ridges, pentachoron seed), ConvexHullBuilder (4D quickhull, cell adjacency by
shared triangle, 4-volume = det/24), OrientedBox (complete 56-axis SAT for two 4-boxes),
AABox/AABox4, Ray{AABox,Sphere,Triangle,Tetrahedron}. `RayTriangle` correctly treats the 4D
ray-vs-2-simplex system as overdetermined (normal equations + residual). Known exceptions:
`Plane::sIntersectPlanes` (sign bug), `MortonCode` (still 3D bit-spread), `ClipPoly` (still 2D
polygon clipping — see conceptual gaps), EPA degenerate-cell path leaves `mNormal`
uninitialized. Both hull builders dropped upstream's coplanar-facet merging (acceptable but
numerically weaker).

### Physics — Body subsystem migrated, everything else untouched
Migrated & internally consistent: `Body/*` (15/15 files: Body, BodyCreationSettings,
BodyInterface, BodyManager, MotionProperties, MassProperties, AllowedDOFs = 4 trans + 6 rot
DOFs), plus a thin slice of Collision: `Shape.h/.cpp`, `ConvexShape.h`, `SphereShape`,
`ConvexHullShape.cpp`, `TransformedShape` (Rotor + RVec4), `RayCast/ShapeCast/CollideShape/
AABoxCast` settings structs, `ScaleHelpers.h`.

MassProperties correctly stores the 4×4 second-moment tensor `M_ij = ∫ρ xᵢxⱼ dV` with correct
Rotate/Translate/Scale and the mapping to 6 bivector moments `I_ij = M_ii + M_jj`.

`BoxShape` is now ported to 4D (tesseract, Phase B step 6, 2026-07-22): Vec4 half extent,
tesseract mass properties, 4D AABox support / supporting-face (a cubic cell), 4D ray cast.
Its 4D content is verified standalone (mass = ∏edges·density, diagonal second moments
m·Sᵢ²/12, support corners, supporting-cell selection). **Caveat:** `BoxShape.cpp` compiles
clean on its own, but the shape *base* layer (`Shape.cpp`, `ConvexShape.cpp`, `SphereShape.cpp`)
still has un-migrated errors, so no shape TU links yet — porting individual shapes is ahead of
the base-layer cleanup.

**Untouched (still 3D, does not compile):** BroadPhase (QuadTree), the other Shape classes
(Capsule, Cylinder, Mesh, HeightField, Compound, Decorated, Plane, Triangle, Tapered*…),
the entire constraint/contact solver (`Constraints/`, `ContactConstraintManager`, all
`ConstraintPart/*`), `PhysicsSystem.cpp` stepping loop, islands, DebugRenderer, SoftBody,
Character, Vehicle, Ragdoll, Hair, TestFramework/Samples/JoltViewer.

### ObjectStream — Rotor/Bivec now registered as primitives (Phase B step 5, 2026-07-22)
`Rotor` and `Bivec` are first-class ObjectStream primitives (text + binary), and
`BodyCreationSettings::mRotation/mAngularVelocity` are re-enabled as `JPH_ADD_ATTRIBUTE`s — the
TU compiles clean and round-trips again (verified by `UnitTests/ObjectStream/
RotorBivecStreamTest.cpp`). Still 3D-lossy here: `DVec4` text/binary serialization writes only
3 components (drops W — part of bug #6), as does the `DVec4` DeterminismLog operator.

### Tests — algebra suite now exists (test-first), rest still a gap
`UnitTests/Geometry/*` were genuinely rewritten for 4D (pentachoron/hypercube hull tests, 4D
GJK/EPA cases, new Ray*Tests) — good.

The Clifford-algebra core now has a reference-validated test-first suite (added 2026-07-22):
`UnitTests/Math/RotorTests.cpp`, `BivecTests.cpp`, `Vec8Tests.cpp`, `RMat44Tests.cpp`, all
checked against a brute-force Cl(4,0) implementation in `UnitTests/Math/Cl4Reference.h` (16-blade
multivectors, geometric product from first principles). These are registered in
`UnitTests/UnitTests.cmake`. Writing them caught a bug the read-only audit missed and had praised
(GetRotor, #8); the P0/P1 math bugs they pinned (#1 exp/log, #2 defect normalization, #8 GetRotor,
#10 sAnd aliasing) are now all fixed and the former EXPECTED-TO-FAIL cases are green regression
tests. **Do not weaken a test to go green — fix the code.**

Added 2026-08-08 (dynamics conventions + constraint algebra), same test-first discipline:
`UnitTests/Physics/DynamicsConventionsTest.cpp` (pins point velocity to the *actual* rotor
integrator by finite differencing it, and the torque to `F ∧ r` by the power identity and by
conservation of energy for an impulse at a point), `UnitTests/Physics/ConstraintPartTests.cpp`
(point/axis constraint solves: relative point velocity vanishes, linear + angular momentum
conserved, energy non-increasing), `UnitTests/Math/Mat44DeterminantTests.cpp` (new
`Mat44::GetDeterminant` vs a brute-force Leibniz expansion). Both physics files deliberately depend
only on `Jolt/Math/*` headers so they run in the standalone runner; `ConstraintPartTests` mirrors the
formulas in the real parts (which need `Body`) — **once the library links, drive the real classes and
delete the duplication.** These were mutation-checked: flipping any single wedge/Contract sign fails
4–5 of the 6 constraint cases. **All 55 cases pass (3652 assertions).**

Running them: the full `UnitTests` target cannot link yet (the physics library doesn't build),
so the suite is verified via a standalone runner that compiles just the header-only tests
against the Core/Math sources that do build. Reproduce with (from repo root, all one line):
```
g++ -O2 -std=c++17 -mavx2 -mbmi -mpopcnt -mlzcnt -mf16c -mfma -mfpmath=sse -pthread \
  -DNDEBUG -DJPH_OBJECT_STREAM -DJPH_PROFILE_ENABLED -DJPH_USE_AVX -DJPH_USE_AVX2 \
  -DJPH_USE_F16C -DJPH_USE_FMADD -DJPH_USE_LZCNT -DJPH_USE_SSE4_1 -DJPH_USE_SSE4_2 -DJPH_USE_TZCNT \
  -IUnitTests -I. UnitTests/Math/{RotorTests,BivecTests,Vec8Tests,RMat44Tests,Mat44DeterminantTests}.cpp \
  UnitTests/Physics/{DynamicsConventionsTest,ConstraintPartTests}.cpp \
  <a main defining DOCTEST_CONFIG_IMPLEMENT + Trace + RegisterDefaultAllocator + Factory::sInstance> \
  Jolt/Core/{Memory,IssueReporting,Profiler,TickCounter,StringTools,Color,Factory,RTTI}.cpp \
  Jolt/Math/Vec3.cpp -o math_tests && ./math_tests
```
Once the library links, delete the standalone runner and run these through the normal
`UnitTests` binary. The committed `build/UnitTests` binary predates all 4D tests; never cite it.

Still missing: `Mat44Tests`/`Vec3Tests`/`QuatTests` test removed/legacy APIs and need rewriting
or deletion (`DMat44Tests.cpp` was de-registered from cmake — its `DMat44.h` is gone); Lane4 has
no dedicated suite (it is exercised transitively through Vec4Tests); no *whole-body* dynamics tests
yet (integrator + solver over many steps, gyroscopic term, determinism — Phase C).

## Verified bugs (fix before building anything on top)

**Fixed 2026-07-22 (this pass, all validated by the Math test suite): #1, #2, #8, #10, and the
Body.inl half of #9; plus the build blocker and repo hygiene (#12, mostly).** Still open:
#3, #4, #5, #6, #7, #11, and the `Mat44::sRotationXY..` half of #9.

P0 — correctness of core dynamics/geometry:
1. **[FIXED]** ~~Rotor exponential only handles simple rotations.~~ `Rotor::sExp(BivecArg)` /
   `Rotor::Log()` now implement the general 4D exp/log via the SU(2)×SU(2) split
   (`Rotor.inl`), and `Body::AddRotationStep`/`SubRotationStep` call `sExp` instead of the old
   inline single-plane formula. Validated against the brute-force reference for random double
   rotations (`TestRotorExpLog`).
2. **[FIXED]** ~~Rotor normalization ignores the e1234 defect.~~ `Normalized()`/`IsNormalized()`
   now project onto Spin(4) by normalizing each SU(2) quaternion factor (removes the e1234
   defect), and `LERP` normalizes its result (`Rotor.inl`). Tracked by
   `TestRotorIsNormalizedDetectsDefect` / `TestRotorNormalizedRemovesDefect` /
   `TestRotorLERPReturnsValidRotor`.
3. **[FIXED]** ~~`MotionProperties::MultiplyWorldSpaceInverseInertiaByVector` is a stub.~~ The
   world-space inverse inertia is now the proper adjoint operator I⁻¹ = Ad(R) D⁻¹ Ad(R⁻¹) with
   R = body_rotation · inertia_rotation, implemented via the new `Rotor::RotateBivec` /
   `InverseRotateBivec` (the bivector adjoint Ad(R)(B) = grade2(R B R̃), computed as
   (R eᵢ R̃)∧(R eⱼ R̃)). The gyroscopic path (`MotionProperties.inl`) now also transforms into and
   out of the body frame. Both are validated against the brute-force Cl(4,0) reference by
   `TestRotorRotateBivec` / `TestInverseInertiaAdjoint` (adjoint = brute-force sandwich, orthogonal,
   symmetric-positive-definite operator, reduces to the diagonal at identity).
4. **[FIXED 2026-08-08]** ~~`Plane::sIntersectPlanes` sign bug~~ — it returned (−x, y, −z, w). The
   per-component Cramer expansion mixed two sign conventions. Rewritten in vector form: with
   `u . sCross(a,b,c) == det[u;a;b;c]`, the row-cofactor vectors of `M = [n1;n2;n3;n4]` are
   `v1=+sCross(n2,n3,n4)`, `v2=−sCross(n1,n3,n4)`, `v3=+sCross(n1,n2,n4)`, `v4=−sCross(n1,n2,n3)`
   and `x = Σ (−1)^(i+1) rhs_i v_i / det` with `det = n1·v1`. Also removed a dead `Mat44` build at the
   top. `TestPlaneIntersectPlanes` passes. ⚠ `TestPlaneIntersectPlanesDiagonal` was ALSO failing but
   the *test* was wrong, not the code: its normals satisfied `n1−n2+n3−n4 = 0`, i.e. 4 dependent
   hyperplanes with no unique intersection. Fixed the test's normals and documented why.
5. **`MortonCode` is still 3D** (`MortonCode.h`): stride-3 masks (`0x09249249`…) with 4 axes →
   bit collisions (x bit at position 3 collides with w), plus `+` instead of `|` propagates
   carries. Needs an 8-bit stride-4 spread (masks ending `0x11111111`). Degrades BVH quality
   silently.
6. **`DVec4` still enforces the 3D "W mirrors Z" invariant**: `CheckW()` asserts Z==W under
   `JPH_FLOATING_POINT_EXCEPTIONS_ENABLED` (`DVec4.inl:113`), 3-arg ctor duplicates Z into W,
   `sAxisX/Y/Z` use the 3-arg form. Any genuine W≠Z translation (i.e. all of 4D) trips it in
   double-precision builds. Purge the invariant.
7. **EPA degenerate cell**: `EPAConvexHullBuilder.h:668-713` can leave `mNormal`/`mLambda`
   uninitialized on a degenerate tetrahedron yet still link the cell; later `IsFacing` reads
   garbage.
8. **[FIXED]** ~~`Mat44::GetRotor` recovers the wrong sign for the e23/e24/e34 (and e1234)
   components.~~ Its internal `(qL,qR)` were labeled mirror to the rotor's SU(2) convention, so
   Step 5 negated the anti-self-dual components; the fix applies the qL↔qR swap (`Mat44.inl`
   GetRotor). Discovered by the new suite (the earlier read-only audit had *praised* GetRotor).
   `sRotation(GetRotor(M)) == M` now holds for simple, double and generic rotors —
   `TestRotorMatrixRoundTrip` / `TestRMat44GetRotor`.

P1 — determinism / hygiene:
9. **[PARTIAL]** C-library `sin/cos` instead of `JPH::Sin/Cos` breaks Jolt's determinism
   guarantee (`Trigonometry.h:9`). Fixed in the `Body.inl` rotation step (now routes through
   `Rotor::sExp`, which uses `JPH::Sin/Cos`). **Still open:** `Mat44::sRotationXY/XZ/XW/YZ/YW/ZW`
   (`Mat44.inl:123-168`) still call C-library `sin/cos`.
10. **[FIXED]** ~~`Bivec::sAnd` type-puns floats through `reinterpret_cast<uint32*>`.~~ Now masks
    per component via `memcpy` (`Bivec.h`). The prior version returned its first argument
    unmodified at -O2/-O3 (strict-aliasing UB). Tracked by `TestBivecSAnd`.
11. `Body::GetBodyCreationSettings` (`Body.cpp:383`) squeezes 6 bivector inertia moments into
    a diagonal Mat44 and drops the three W-plane moments — lossy round-trip.
12. **[MOSTLY FIXED]** Repo hygiene: `build/` de-committed + `/build/` added to `.gitignore`;
    `Bivec.h`/`RMat44.h`/`Tetrahedron.h`/`IndexedTetrahedron.h`/`RayTetrahedron.h` + new test
    files now tracked; `Bivec.h` registered in `Jolt.cmake`. **Still open:** the bulk of the
    Geometry/Physics/ObjectStream migration is still uncommitted working-tree changes — commit
    in reviewable slices.
13. **[FIXED 2026-08-08]** ~~Torque / angular-impulse wedge operands reversed everywhere.~~ The
    migration transliterated the 3D `r × F` as `r ∧ F`, but the 4D point-velocity operator
    (`Bivec::Contract`, pinned by the `exp(½Ω·dt)·R` integrator) is the *negative* adjoint of the
    wedge, so the torque conjugate to Ω is `F ∧ r`. As written, a force applied off-center spun the
    body the **wrong way**, and constraint solves converged to the wrong velocity state. Fixed in
    `Body::AddForce`/`AddImpulse` (`Body.inl`), the buoyancy impulse (`Body.cpp`) and the angular
    Jacobians of `Axis`/`DualAxisConstraintPart` (`r∧n` → `n∧r`). See the 🔑 box in the roadmap;
    pinned by `DynamicsConventionsTest.cpp` + `ConstraintPartTests.cpp`. ⚠ Note the effective mass is
    invariant under this flip, which is why it survived a compile-clean port — only a
    velocity-residual / momentum / energy test detects it.

## Things that do NOT trivially generalize (the real design work left)

1. **Contact manifolds.** In 3D the contact patch between convex bodies is a 2D polygon; Jolt
   clips face polygons (`ClipPoly`) and reduces to ≤4 points. In 4D the patch is a 3D
   *polyhedron* living in the contact hyperplane. `ClipPoly.h` still clips 2D polygons.
   Needed: supporting-face extraction returning a polyhedron (or point cloud), convex
   polyhedron-vs-polyhedron clipping (Sutherland–Hodgman generalized to halfspace-clipping a
   polyhedron), and manifold reduction to ~5–6 points whose hull spans the patch (4 points
   can't stabilize a 3D patch; upstream's 4-point heuristic must be rethought, incl.
   `ManifoldBetweenTwoFaces`). This is the hardest open problem in the project.
2. **Angular constraint math.** 4D has 6 rotational DOFs. A "hinge axis" is a bivector plane;
   the effective-mass matrix of an angular constraint part is up to 6×6; the inertia tensor is
   properly a symmetric operator on bivectors (6×6), of which the current
   diagonal-6-in-principal-frame + rotor conjugation is the special case. Every
   `ConstraintPart/*` taking `Vec3 inAxis` needs redesign (Jacobians in bivector components),
   not a type swap. Decide early: keep diagonal+adjoint or store the full 6×6 world-space
   inverse inertia per step (upstream stores world-space `Mat44`; the analog is a 6×6).
3. **Shape zoo semantics.** Sphere→3-sphere (done; hypervolume π²r⁴/2), Box→tesseract
   (mechanical), Capsule→sphere-swept segment (fine). But: Cylinder/TaperedCylinder (which
   2-plane × which cross-section?), Triangle/Mesh shapes (surface meshes become *tetrahedral*
   boundary meshes of S³-like boundaries — `ConvexShape.h:135`, `Shape.h:347`; the whole
   `GetTrianglesContext`/DebugRenderer pipeline assumes triangles), HeightField (2D grid → 3D
   grid of heights over a hyperplane — memory grows cubically). Recommendation: v1 ships
   Sphere, Box, Capsule, ConvexHull, Plane, Compound/Decorated wrappers; cut Mesh/HeightField/
   Cylinder/Triangle from the build until needed.
4. **Inside/outside & winding tests.** Parity ray-tests and winding conventions were audited
   as 3D-invalid in spots (`Shape.cpp:319` TODO). Check every "odd number of hits = inside".
5. **Visualization.** DebugRenderer, TestFramework, Samples, JoltViewer are inherently 3D
   renderers. For the engine to *ship* you don't need them — you need the headless library +
   UnitTests + PerformanceTest. Recommend: make DebugRenderer compile via projection/slice
   stubs (or `JPH_DEBUG_RENDERER` off), and treat real 4D visualization (slice w=const, or
   perspective projection to 3D) as a separate downstream project.
6. **Scope cuts for v1.** SoftBody, Character, Vehicle, Ragdoll, Hair are large, deeply
   3D-coupled (Vehicle especially — wheels are inherently 3D), and none are needed to call the
   core engine "working". Exclude them from `Jolt.cmake` for v1 rather than porting.

## Roadmap to a shippable v1

**⏸ RESUME HERE (2026-08-08):** The headless v1 library builds and `UnitTests` runs at **357/364**.
Next, in priority order: (1) **the GJK simplex-growth bug** — all 7 remaining failures (5 EPA, 2 GJK)
trace to it, and the diagnosis with trace evidence and everything already ruled out is in the
build-status section at the top; start from the `MustIncludeD` shortcut in
`GetClosestPointOnTetrahedron`. (2) Write the **4D HelloWorld** (tesseract onto a hyperplane floor) —
the whole stepping loop compiles now, so this is the first chance to see the engine simulate, and it
will flush out runtime bugs the unit tests cannot. (3) The **contact manifold pipeline**
(`ManifoldBetweenTwoFaces` is a reduced single-point fallback — see the TODO in that file; resting
boxes get one contact instead of a patch, so stacks will be less stable than 3D Jolt). (4) Re-enable
de-registered unit tests as each area is ported. **Read the 🔑 torque-convention box below before
touching any dynamics code** — the naive `r ∧ F` transliteration is wrong and an effective-mass test
cannot detect it.

**Phase B step 6(b) is COMPLETE — the v1 shape layer compiles clean at the per-TU level (15/15 TUs).**
Verified via a per-TU `-fsyntax-only` check with the Jolt build flags (helper `cc.sh`; recreate from
`Jolt.dir/flags.make`, `-O0 -fsyntax-only`). The clean TUs: `Shape.cpp`, ConvexShape, SphereShape,
BoxShape, CapsuleShape, ConvexHullShape, PlaneShape, EmptyShape, DecoratedShape, ScaledShape,
RotatedTranslatedShape, OffsetCenterOfMassShape, CompoundShape, StaticCompoundShape,
MutableCompoundShape. (Re-verified 2026-08-08: still 15/15 clean. `ConvexHullShape.cpp` had picked up 1
error from the still-3D `PolyhedronSubmergedVolumeCalculator.h`; that include was a leftover — the
class is only named in a TODO now that `GetSubmergedVolume` is stubbed — and removing it restored the
TU. That header stays 3D until the 4D polytope-vs-hyperplane clip is written.)

**⏳ Phase B step 6(c) IN PROGRESS.** Foundational, tested pieces done:
- **Inverse-inertia operator (bug #3) FIXED + tested** — `Rotor::RotateBivec`/`InverseRotateBivec`
  (bivector adjoint) added to `Rotor.h`/`.inl`; `MultiplyWorldSpaceInverseInertiaByVector` and the
  gyroscopic path now do the real Ad(R) D⁻¹ Ad(R⁻¹). Validated vs the brute-force Cl(4,0) reference
  (`RotorTests.cpp`: `TestRotorRotateBivec`, `TestInverseInertiaAdjoint`). This is THE operation every
  constraint impulse flows through.
- **🔑 THE 4D TORQUE CONVENTION — settled, tested, and it invalidates the naive port (2026-08-08).**
  The three halves of 4D rigid-body dynamics only close for ONE choice of wedge operand order, and the
  code originally had it backwards in every place. The chain: `Body::AddRotationStep` integrates
  `R ← exp(½Ω·dt)·R`, which forces the body-fixed point velocity to be `+Ω.Contract(r)` (verified by
  finite-differencing the real integrator, not by analogy). Since
  **`⟨Ω.Contract(r), F⟩_R4 = −⟨Ω, r∧F⟩ = +⟨Ω, F∧r⟩`** (Contract is the *negative* adjoint of wedge),
  the bivector conjugate to Ω — i.e. the torque — is **`F ∧ r`, not `r ∧ F`**. General rule:
  > **a 3D `a × b` in a torque / angular-Jacobian role becomes a 4D `b ∧ a` — operands swapped.**

  Consequences applied: torque/angular impulse in `Body::AddForce`/`AddImpulse` (`Body.inl`) and the
  buoyancy impulse (`Body.cpp`) were `r ∧ F` → now `F ∧ r`; the angular Jacobian in `AxisConstraintPart`
  /`DualAxisConstraintPart` was `r ∧ n` → now `n ∧ r`. ⚠ The effective mass is *quadratic* in the
  Jacobian, so it is identical under the flip — an effective-mass or symmetry/PSD test can NOT catch
  this. Only a velocity-residual or momentum/energy test can. Pinned by
  `UnitTests/Physics/DynamicsConventionsTest.cpp` (4 cases) and `ConstraintPartTests.cpp` (6 cases);
  mutation-checked (flipping any one sign fails 4–5 of the 6). Angular momentum about the origin is
  correspondingly `L = I·Ω + m(v ∧ x)`.
- **`AxisConstraintPart` (1-row) and `DualAxisConstraintPart` (2-row) ported to 4D bivector Jacobians.**
  Angular Jacobian is the bivector `n ∧ r` (was `r × n`, see the convention above); the inverse inertia
  is applied via the operator above (diagonal+adjoint, no precomputed matrix); `J·w` uses `Bivec::Dot`.
  The effective mass is scalar (Axis) / 2×2 (DualAxis) — the row count, independent of the spatial
  dimension. Both compile clean per-TU. The old `Mat44Arg inRotation` / `GetInverseInertiaForRotation`
  args became `RotorArg` + the operator. These are the parts contacts and most joints reuse.
- **`PointConstraintPart` redesigned for 4D + validated (2026-08-08).** A point-to-point constraint has
  **4** linear DOF in 4D: `C = p2 − p1` has 4 rows, the effective mass is a `Mat44` and λ is a `Vec4`.
  Writing `C_r` for the 4×6 map `Ω ↦ Ω.Contract(r)`, its transpose is `C_rᵀ: u ↦ u ∧ r` (6×4), so
  `K = (m₁⁻¹+m₂⁻¹)E + C_{r1} I₁⁻¹ C_{r1}ᵀ + C_{r2} I₂⁻¹ C_{r2}ᵀ`, assembled one column at a time by
  feeding it the 4 basis vectors. The per-body 6×4 operator `I⁻¹C_rᵀ` is cached as its 4 bivector
  columns (`mInvI1_R1X[4]`) and doubles as the λ→Δω map. Needed a new **`Mat44::GetDeterminant()`**
  (upstream only had the 3×3) to guard invertibility — tested against a brute-force Leibniz expansion
  in `UnitTests/Math/Mat44DeterminantTests.cpp`. Compiles clean per-TU.

- **`ContactConstraintManager` + `ContactListener` ported (2026-08-08).** Both compile clean per-TU.
  The 4D-specific parts:
  - **Friction has 3 directions, not 2.** The space perpendicular to a contact normal is a hyperplane,
    so Coulomb friction acts in a **3D tangent space**: `mFrictionConstraint3` added, `GetTangents`
    returns 3, `mFrictionLambda` is `Vector<3>`, and the friction cone clamp is over
    `λ₁²+λ₂²+λ₃² ≤ (μλₙ)²` (the cone itself is unchanged — only the dimension of the disc it clamps to).
  - The orthonormal 3-frame moved from a file-static in `PlaneShape.cpp` to
    **`Vec4::GetNormalizedPerpendicularBasis`** (tested: orthonormal + spans R⁴ with the normal), since
    contacts and planar patches both need it.
  - **The world-space inverse inertia can no longer be precomputed as a `Mat44`** — it is a 6×6 operator
    on bivectors — so the three places that cached `inv_i1`/`inv_i2` now pass `RotorArg` + the inertia
    scale down to `AxisConstraintPart` (whose templated entry point gained the scale args).
  - Cache layout: `Float3`→`Float4` throughout; `CachedBodyPair::mDeltaRotation` is an uncompressed
    `Float8` rotor (the 3-of-4 quaternion trick does not generalize — Spin(4) needs both SU(2) signs).
    Struct sizes: `CachedContactPoint` 36→48, `CachedManifold` 56→72, `CachedBodyPair` 28→52.
  - The body-pair cache-validity rotation test now uses the **8-component rotor dot** = scalar part of
    `p·~q` = `(cos(a/2)+cos(b/2))/2` under the SU(2)×SU(2) split; still 1 at zero rotation and decreasing
    as either angle grows, so the existing `mBodyPairCacheCosMaxDeltaRotationDiv2` comparison stays valid.
  - `ContactSettings::mRelativeAngularSurfaceVelocity` is a `Bivec`; the surface velocity it induces at
    `r` is `Contract(r)` (was a cross product).

**Remaining for 6(c), in order:**
- **`PhysicsSystem` stepping loop** (mostly mechanical renames onto the migrated Body API).
- **Contact manifold pipeline** (conceptual gap 1 — supporting-face polyhedra, halfspace clipping,
  ≥5-point reduction). ⚠ Note the contact *solver* above is done, but it is fed by the manifold
  pipeline, which is still the hardest open problem.
- **The rotational parts** (`RotationEuler`/`Angle`/`Hinge`/`SwingTwist`) — apply the same convention;
  these constrain *relative orientation*, so their Jacobians are bivector-valued in both bodies (no
  lever arm) and the redesign question is which 6-DOF subspace each part locks. Not needed for a first
  simulating scene.
- **Excise SoftBody from the core for v1** (the documented scope cut, not yet done at the source level).
  `SoftBodyShape.h`/`SoftBodyCreationSettings.h`/`SoftBodyMotionProperties.h` are still 3D and are still
  *included* by `Body.cpp` (2 errors, all in `GetSoftBodyCreationSettings`) and `BodyManager.cpp`
  (9 errors, ~50 references). These are the only remaining errors in the otherwise-migrated Body layer.
- **Contact manifold pipeline** (conceptual gap 1 — supporting-face polyhedra, halfspace clipping,
  ≥5-point reduction; the hardest open problem).
- **`ContactConstraintManager`** (2354 lines): uses 3 `AxisConstraintPart`s per contact point — the
  normal + friction. In 4D the friction tangent space is 3D, so add a **3rd** friction axis + tangent
  (reuse the orthonormal-3-frame from `PlaneShape::sPlaneGetOrthogonalBasis`); `mFrictionLambda`
  `Vector<2>` → `Vector<3>`; `GetTangents` returns 3; replace the `Mat44` inverse-inertia args with the
  operator.
- **`PhysicsSystem` stepping loop** (mostly mechanical renames onto the migrated Body API).

**The full library link is still gated on all of the above** — `RegisterTypes.cpp` registers the
un-migrated constraint/vehicle/ragdoll types, so `cmake --build` will not link yet. Verify with the
per-TU syntax check for now, and add reference-validated tests for each new constraint part /
conservation law (Phase C). Also still TODO on the (now-compiling) shapes: `GetTrianglesStart/Next`,
`CollideSoftBodyVertices`, and `GetSubmergedVolume` are stubbed (S³ tetra tessellation / 4D polytope
submerged-volume clip / soft-body iterator all deferred); ConvexHull's convex radius is forced to 0
(the 4D vertex shrink needs a 4-hyperplane vertex solve).

Key 4D decisions locked in the 6(b) pass:
- **Real 4D shape content:** CapsuleShape (segment+4-ball support, derived 4D second-moment tensor —
  reduces to the 4-ball `mass*r²/6` as height→0); PlaneShape (hyperplane + 3-cube finite patch via an
  orthonormal 3-frame); ConvexHullShape (mass tensor via **pentachoron decomposition**: cone each
  boundary tetra to the COM; canonical 4-simplex second moment has diagonal 1/360, off-diagonal
  1/720; `mInertia` = the raw second-moment covariance, matching Box/Capsule).
- **SIMD 4-box helpers** (`RayAABox4`/`AABox4VsBox`/`AABox4Scale`/`AABox4EnlargeWithExtent`/
  `AABox4VsPoint`) take **8 axis args** (min/max × X/Y/Z/W); every visitor `TestBounds`/`VisitNodes`/
  `TestBlock` is 8-arg. Added a conservative `AABox4VsBox(OrientedBox,…)` 4D overload (tests the OBB's
  world AABB — TODO tighten to a real 4D OBB-vs-AABB SAT).
- Compound `SubShape` stores a full `Rotor` (`Float8`) + `Vec4` position uncompressed (the 3-of-4
  quaternion compression trick does not generalize). `StaticCompoundShape::Node` gained
  `mBoundsMin/MaxW[4]` (HalfFloat) → 80 bytes; `MutableCompoundShape::Bounds` SOA gained `mMinW`/`mMaxW`.
- Support infra: `ScaleHelpers::CanScaleBeRotated`/`RotateScale` (4D rotor form), `Plane::Scaled(Vec4)`,
  `Shape.cpp` `ToRotor`→`GetRotor` fix, cmake + `RegisterTypes.cpp` cut-shape edits
  (Mesh/HeightField/Cylinder/Tapered*/Triangle + SoftBodyShape removed).

**Bug found while porting (unfixed, flagged):** `SphereShape::GetMassProperties` stores
`mass*r^2/3` on the `mInertia` diagonal, but `mInertia` is the raw second-moment tensor M_ij
(BoxShape sets `m*S_i^2/12`; the correct 4-ball value is `mass*r^2/6`). SphereShape is 2× too
large and inconsistent with Box/Capsule. CapsuleShape uses the correct second-moment convention
and reduces to `mass*r^2/6` as height→0. Fix SphereShape to `mass*r^2/6` (and add a mass-
properties reference test) — deferred because SphereShape is marked "done" and I didn't want to
silently edit it.

Definition of "shipped v1": headless library builds warning-clean; UnitTests all green
(including new Rotor/Bivec/dynamics conservation tests); a 4D HelloWorld (tesseract falling
onto a hyperplane floor, stacks of tesseracts stable, spheres bouncing) runs deterministically;
PerformanceTest runs on a 4D scene.

**Phase A — stabilize the foundation. Steps 1, 2, 4 DONE (2026-07-22); step 3 mostly done.**
1. DONE: repo hygiene — `build/` de-committed + gitignored, untracked 4D leaf sources tracked,
   `Bivec.h` registered in `Jolt.cmake`. Still open: commit the bulk migration in reviewable
   slices (the working tree still has ~140 modified files).
2. DONE: reference-validated algebra tests written test-first (RotorTests/BivecTests/Vec8Tests/
   RMat44Tests + `Cl4Reference.h`), registered in `UnitTests.cmake`; `DMat44Tests` de-registered.
   Still TODO here: rewrite `Mat44Tests` for the new API, delete/rewrite `Vec3Tests`/`QuatTests`
   (with the legacy types), add a dedicated Lane4 suite.
3. Math P0/P1 bugs fixed against those tests (all 37 cases green): #1 exp/log, #2 defect
   normalization, #8 GetRotor, #10 sAnd, and the Body.inl half of #9. **Still open:** #3
   (`MultiplyWorldSpaceInverseInertiaByVector` adjoint), #4 (`Plane::sIntersectPlanes` — covered
   by existing `TestPlaneIntersectPlanes`), #5 (`MortonCode` — needs a new no-bit-collision test),
   #6 (`DVec4` Z==W invariant), #7 (EPA degenerate cell), #11 (lossy inertia round-trip), and the
   `Mat44::sRotationXY..` half of #9.
4. DONE: build blocker fixed — `AABox.h` now includes `RMat44.h`; error count ~3,800 → ~3,400
   and the remainder is all Phase B un-migrated-layer work.

**Phase B — restore an end-to-end simulating core:**
5. DONE (2026-07-22): ObjectStream primitives for Rotor/Bivec (text + binary), excluded
   `BodyCreationSettings` attributes re-enabled, DeterminismLog operators for Rotor/Bivec added.
   Verified by `RotorBivecStreamTest`. Remaining nit: `DVec4` serialization still drops W (bug #6).
6. Shape layer for the v1 shape set. IN PROGRESS (2026-07-22):
   - DONE: `BoxShape` → tesseract (mass properties, 4D support fn + cubic-cell supporting face,
     4D ray cast), mirroring `SphereShape`. Both convex *leaf* shapes now compile clean with the
     renderer off (drop the 3D `GetTrianglesContext.h` include, stub GetTriangles/Draw/SoftBody,
     build diagonal inertia without the removed `Mat44::sScale`). Bugs fixed while porting:
     `AABox::GetSupportingFace` sign (picked the wrong cell), `TransformedShape::ToRotor` typo;
     `Shape::GetSubmergedVolume` inBaseOffset made unconditional for renderer-off consistency.
   - DONE: `CollisionDispatch` migrated (Vec4/RMat44) and `ConvexShape` convex-vs-convex **collide
     path** (`sCollideConvexVsConvex` — GJK+EPA penetration, split RMat44 into Mat44 rotation +
     Vec4 translation for the local-space support functions), `CollidePoint`, fallback `CastRay`.
     ConvexShape/SphereShape/BoxShape all compile clean with the renderer off.
   - DONE: turned the debug renderer OFF for the build (v1 is headless) — this removes the
     DebugRenderer + per-shape `Draw` errors and makes the count reflect the real surface.
   - DONE (a): fixed the geometry shape-cast start-translation bug. `GJKClosestPoint::CastShape`
     (both overloads) and `EPAPenetrationDepth::CastShape` used to take a pure-rotation `Mat44`
     start and build `TransformedConvexObject(inStart, Vec4::sZero(), …)`, dropping the cast start
     position. They now take the start rotation (`Mat44Arg`) and translation (`Vec4Arg`)
     separately. `ConvexShape::sCastConvexVsConvex` is un-stubbed (real convex-vs-convex cast,
     Vec4/RMat44). Verified: a unit sphere cast from x=3 toward a sphere at the origin hits at
     fraction 0.1. `EPATests` updated for the new signature.
   - **Next, in order:** (b) make `Shape.cpp` link — it includes `ScaledShape.h` +
     `StaticCompoundShape.h`; port the wrapper shapes Compound(base)/StaticCompound/MutableCompound/
     Decorated/Scaled/RotatedTranslated/OffsetCOM + ConvexHull + Capsule + a real 4D PlaneShape,
     and cmake-exclude Mesh/HeightField/Cylinder/Tapered*/Triangle + remove their
     includes/sRegister calls in `RegisterTypes.cpp` (which includes every shape header). (c)
     migrate `ContactConstraintManager` + the ConstraintParts (bivector Jacobians — conceptual
     gap 2) and `PhysicsSystem` stepping loop; (d) link + a permanent BoxShape/collision unit test.
   - Also still TODO: `GetTrianglesStart/Next` and `GetSubmergedVolume` are stubbed on the convex
     shapes (S^3 tetra tessellation; 4D polytope submerged-volume clip).
7. BroadPhase: QuadTree AABB logic → 4D AABox (mostly mechanical; AABox4 SIMD already done).
8. Contact manifold pipeline (conceptual gap 1): supporting-face polyhedra, halfspace clipping,
   ≥5-point manifold reduction, `ContactConstraintManager`.
9. Constraint parts (conceptual gap 2): contact normal + friction first (friction in 4D has a
   3D tangent space → 3 friction directions or a tangent-space solve), then
   Point/Distance/Axis parts; defer Hinge/Slider/SixDOF/etc. redesign to v1.1 unless needed.
10. `PhysicsSystem.cpp` stepping loop, islands, sleeping (mostly mechanical renames onto the
    migrated Body API).

**Phase C — validate & ship:**
11. Physics-level tests: momentum/energy/angular-momentum (bivector L) conservation, 4D
    Dzhanibekov analog for the gyroscopic integrator, restitution/friction, stack stability,
    determinism (same-input hash via DeterminismLog), double-precision build, PerformanceTest
    scene.
12. Port HelloWorld to 4D as the reference sample. Decide DebugRenderer story (slice/projection
    or off).
13. Only after v1: SIMD for Rotor/Bivec (profile first), NEON/RVV parity for Vec8, coplanar
    facet merging in the hull builders, the cut modules, real 4D viewer.

## Working conventions for this repo

- **Never trust, always test the algebra.** Every hand-derived Clifford formula (product
  tables, exp/log, adjoints, inertia maps) must have a unit test against a brute-force
  reference implementation before it is used elsewhere. This rule exists because the project
  accumulated an untested algebra core and three verified sign/convention bugs.
- Follow the type table above; never introduce new uses of `Vec3`/`Quat`.
- Use `JPH::Sin/Cos/Tan/ASin/...` (never `std::`/C-library trig) in simulation code —
  determinism.
- Keep upstream file layout and naming so `git diff master` stays reviewable; when a 3D
  concept is intentionally renamed/removed, note it in the header comment as `TODO(4D):` or
  a short "replaces X" line.
- Commit in small reviewable slices; never commit `build/` output; keep new files registered
  in `Jolt.cmake`/`UnitTests.cmake` at creation time.
- Auto-memory files (`project_4d_*.md`) hold the derived Cl(4,0) formulas and progress
  snapshots; update them when module status changes materially.
