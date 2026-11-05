# Four-dimensional double positions and Morton keys

## Contracts

- `DVec4` has four independent doubles on every backend. XYZ-only compatibility
  constructors embed with W = 0. Explicit `Vec3`/`Double3` projections are still
  deliberately three-dimensional; spatial arithmetic uses `Lane4`/`Vec4` instead.
- Equality, comparisons, truth masks (four bits), NaN detection, reductions,
  hashes, indexed writes and float conversions include W. The historical
  `ToVec3RoundDown/Up` names are aliases returning **Vec4**, not XYZ projections.
- `RMat44` still stores a float linear map. Applying it to double positions,
  composing translations, inversion and post-scaling no longer narrow positions
  to floats before doing the calculation. This does not make its rotation matrix
  double precision, or change the float velocity representation.
- Core binary streams store four doubles per `DVec4`, both scalar and array.
  They are unversioned: applications must version their enclosing state format.
- Text and binary ObjectStream store four doubles, including arrays. The text
  token is `dvec4`; scalar precision conversions are `Vec4` <-> `DVec4`.
  ObjectStream is now **2.00** and rejects 1.00, rather than reading a legacy
  three-double record's following field as W. There is no legacy stream importer.
- Determinism logging writes the bit patterns of all four doubles, without
  strict-aliasing type punning.

## Morton range and ordering

`MortonCode` encodes a **32-bit key with 8 bits per axis**. Each coordinate is
normalized within its bounds and quantized to the nearest integer in [0, 255].
Bit `4*b + 3/2/1/0` contains X/Y/Z/W bit `b`, respectively. The bit sets do not
intersect; they are combined with OR, not addition. Zero-extent axes map to zero.
Distinct quantized four-tuples have distinct keys. Different floating points
inside the same quantization cell intentionally share a key.

The tests cover all 256 spread values, an exhaustive 16^4 grid, 10,000 seeded
full-width tuples against a separate bit-by-bit encoder, the old X/W collisions,
maximum values, translated bounds and flat bounds. The helper and tests are in
the normal CMake targets. No new world-residency system or BVH speedup is claimed.

## Verified build

WSL Ubuntu, GCC 13.3, headless Release, profiler and LTO disabled, build jobs = 1:

```sh
cmake -S Build -B /home/brkdrd/mover-build-precision \
  -DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON \
  -DDEBUG_RENDERER_IN_DEBUG_AND_RELEASE=OFF \
  -DPROFILER_IN_DEBUG_AND_RELEASE=OFF \
  -DINTERPROCEDURAL_OPTIMIZATION=OFF -DUSE_ASSERTS=OFF \
  -DTARGET_SAMPLES=OFF -DTARGET_VIEWER=OFF \
  -DTARGET_PERFORMANCE_TEST=OFF -DTARGET_HELLO_WORLD=ON
cmake --build /home/brkdrd/mover-build-precision \
  --target UnitTests HelloWorld -j1
/home/brkdrd/mover-build-precision/UnitTests
```

The real double-precision `UnitTests` executable passes **383 cases / 2,473,442
assertions**. Separate focused AVX2, SSE (no AVX) and scalar-backend probes also exercise
four-lane math, serialization and Morton keys (21 cases / 76,245 assertions per
backend) with assertions and invalid,
divide-by-zero and overflow traps enabled. NEON and RVV were not executed.

Regression tests were observed failing before the repairs: independent W writes
and reductions, four-double wire size and sentinel alignment, missing fourth
stream component, narrowing of large double transform coordinates, translated
ray-direction cancellation, Morton bit collisions, precision conversion and
legacy ObjectStream rejection, and the fourth determinism-log field.

An assertions-enabled full GCC Release build encountered an existing optimizer
array-bounds warning in `BoxShape::GetSurfaceNormal` / `AABox::GetSupportingFace`
(via Lane4 component-index helpers); it was not suppressed or changed here. The
full build above uses the project's default `USE_ASSERTS=OFF`; the focused probes
separately test the modified code with assertions on.

## Sleeping repair and remaining boundaries

The initial double HelloWorld ran but exited 1: its tesseract at W approximately
10 did not sleep. A direct probe of `MotionProperties::ResetSleepTestSpheres`
with the first point `(1,2,3,10)` read back `(1,2,3,0)`, proving loss in its
`Double3` offset storage. A separate unit regression reproduced the same loss
with sub-float-resolution coordinates at large offsets.

The authorized storage-only repair changes `mSleepTestOffset` to `Double4`,
loads it with `DVec4(mSleepTestOffset)` and stores it with `StoreDouble4`.
Existing generic SaveState/RestoreState calls now carry the fourth double;
applications must version this changed state record size. After the repair,
double HelloWorld passes every check and sleeps after **63 steps**.

The gyroscopic/clamp code and `Body.inl` are untouched. Independently of this
storage repair, the existing three sleep sample points (COM plus two independent
offsets) leave a rotational plane unobserved in 4D. Correcting that sampling needs
coordinated changes to all point/sphere counts and rotation-specific tests; it is
not part of this patch. Large-world rebasing/velocity precision and cross-platform
determinism are also not established by these tests.
