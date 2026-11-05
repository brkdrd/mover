# Body mass export: supported domain and failure policy

`MassProperties::mInertia` is the symmetric **second-moment** matrix
`M = integral x x^T dm`, not the rotational inertia. In a principal frame,
`I_ij = d_i + d_j`. `MotionProperties` retains six float reciprocals and a
principal-frame rotor; it does **not** retain an independently editable copy of M.
A zero physical planar moment has **zero inverse** by the existing particle
convention. It must not become a small positive moment with a huge reciprocal.

## Checked export, not a universal lossless serialization claim

Use `bool Body::TryGetBodyCreationSettings(BodyCreationSettings &out)` when
export may fail. Failure leaves **all of `out` and the body unchanged**. The
legacy `GetBodyCreationSettings()` calls the checked path and **terminates the
process** on failure, including in release builds. Its return type has no error
channel; returning plausible but corrupted inertia or silently locking DOFs is
not an acceptable fallback. Existing well-conditioned successful callers need
no change. Callers migrating arbitrary bodies between worlds should adopt the
checked API and keep the original body if it fails.

The supported mass domain is finite physical mass state, with normal binary32
moments/reciprocals representable by the existing creation path, satisfying the
following **additional export checks**:

1. Unlocked translation requires positive finite inverse mass with a finite
   representable reciprocal. Fully locked translation exports canonical mass 1.
2. Principal planar inverses must be nonnegative and finite. Reconstructing
   principal second moments must produce nonnegative finite representable values;
   an inconsistent inverse-setter override need not have this representation.
3. **If any planar inverse is exactly zero, the encoded body-local M must be
   exactly diagonal.** Off-diagonal entries are tested against zero, not a
   tolerance. This is a deliberately conservative supported subset, not a claim
   that every other singular float matrix is mathematically unrepresentable.
   Coordinate-aligned rank-one and rank-two tensors are supported. An oblique
   null space is rejected even when a particular rounded matrix happens to work.
   Arbitrary *world body orientation* is independent of this body-local condition.
4. Before the real `SetMassProperties` dry-run, a checked preflight runs the
   same float Jacobi eigensolver, with its **unchanged** normalization/eigenpair
   tolerances, and the same sorted matrix reconstruction test. Checked Jacobi
   first requires every entry finite and `abs(M_ij) <= sqrt(FLT_MAX)/(16*n*n)`
   (`n=4` for mass decomposition). This deliberately conservative intermediate-
   range restriction leaves headroom for float Jacobi arithmetic and squared
   residuals; it is **not** a rank cutoff or relaxed residual tolerance. It also
   checks the double-evaluated theta ratio before the ordinary float division
   and square. Unsupported ranges return false **before** overflow, including
   with traps enabled. Both eigenpair and reconstruction diagnostics explicitly
   require finite squared norms: infinity must not become a passing tolerance.
   Numerical diagnostic or convergence failures return false **in every build**,
   without invoking the assertion callback. Exactly diagonal tensors use the
   decoder's existing diagonal shortcut (and do not need the Jacobi range cap).
   The temporary `MotionProperties` then decodes the proposed settings through
   the real path; zero-inverse count must be unchanged and every decoded inverse
   must remain finite and nonnegative. Checked success retains the ordinary
   solver's float operations and reconstruction criterion for that rerun.
   This adds no exception masking: the solver's existing scoped invalid-operation
   guard is unchanged, and the caller's trap masks are restored. The ordinary
   creation decoder's historical diagnostics/tolerances remain unchanged.
5. Let A and B be the original and decoded **unmasked body-local 6x6 inverse
   response operators**, evaluated on all six coordinate bivectors using the
   actual float rotor/adjoint path. Export requires
   `||B - A||_F <= (64 * FLT_EPSILON) ||A||_F`. Differences and squared norms are
   accumulated in double. Checking the unmasked operator also protects latent
   inertia behind partial world-space DOF locks and static/kinematic motion types.

The last bound is an explicit *output accuracy policy*, `2^-17` relative
Frobenius error, not an eigenvalue cutoff or a mathematically derived error
bound for the eigensolver or the matrix/rotor/adjoint chain. The coefficient
64 is a chosen policy; acceptance is decided by the **measured decoded
response**. It is a normwise response contract, not a relative error guarantee
for every impulse in or arbitrarily close to a null space. No conditioning
threshold turns a small positive moment into zero. If cancellation or
conditioning defeats these checks, export rejects; it does not regularize M
or fabricate sphere inertia.
Existing tighter physics regression tolerances are unchanged.

Static bodies without motion properties and completely locked rotational
sectors export the existing finite canonical, unobservable mass sectors. DOF
masks, pose, COM/shape relationship and velocities retain their existing export
semantics; rejection does not change any of them.

## Why the oblique singular case cannot be promised generally

The review fixture is two unit particles at `+/-(0,0,1,1)`, mass 2, with
`M22 = M23 = M32 = M33 = 2` and all other entries zero. Its principal second
moments are `[4,0,0,0]`, and its planar inverse spectrum is
`[1/4,1/4,1/4,0,0,0]`. The principal frame must pass through a float rotor.
Reconstructing a body-local float matrix rounds its entries independently and
can destroy its exact null space. The former export/recreate path produced
two inverse moments of `8388608` instead of zero. The permanent regression
first reproduced this through actual body creation/export/recreation, without
inverse setters, and checked the analytical world response.

More double-precision arithmetic cannot recover rank information absent from
the final float matrix. Likewise a blanket `epsilon * largest_eigenvalue`
clamp cannot distinguish that lost null space from a genuine skinny positive
mass distribution. General preservation needs a representation carrying
principal rank/factors explicitly through creation and serialization, rather
than guessing rank from rounded matrix entries. That representation/API redesign
is not included here. The checked exporter rejects the fixture instead.

## Exact rank-one creation, without epsilon rank decisions

The original float Jacobi eigenpair assertion divides its residual scale by an
eigenvalue-sized quantity; a null eigenpair has no useful relative scale.
For an **exact rank-one submitted M**, mass decomposition instead certifies
`M_ij * M_pp == M_ip * M_pj` for every i,j, with a positive maximal diagonal
pivot p and exact symmetry. Both binary32 products are evaluated in binary64:
each uses at most 48 significand bits, fitting its 53 bits, with no float-range
underflow/overflow. The test is therefore exact for the supported normal float
inputs, not a numerical rank threshold. A skinny positive diagonal fails the
identity, however small its nonzero entry is.

After certification, the occupied column is assigned **directly from the
normalized pivot column**, preserving structural zero components exactly. A
stable double Householder construction supplies only the other three columns,
with eigenvalues `[trace(M),0,0,0]`. Evaluating the occupied column via the
Householder subtraction instead can leave `2.220446049250313e-16` in a required
zero component: for `M = 2^40 x x^T`, `x=(0,-2,10,-10)`, it previously reconstructed
`M10` as approximately `0.00697405` instead of zero. The existing proper-orientation
and reconstruction assertions are retained unchanged. This fixes the review
fixtures' *initial creation* assertions without weakening the generic eigensolver
or its tolerances. It does not claim to recover an ideal rank-one particle
distribution from an already rounded matrix whose exact rank is different.

`SetMassProperties` also keeps identity as the principal frame when M is
already exactly diagonal. Sorting/permuting it and converting that permutation
to a float rotor is unnecessary and perturbs coordinate-aligned null modes.
The former absolute `IsNearZero` gate is removed: a successful decomposition
of zero M keeps six exact zero inverses, and a successful small positive M
keeps its actual reciprocals rather than acquiring unit-sphere inertia.
For export, an exact zero `I_ab` plus PSD implies `d_a = d_b = 0`; other d's are
recovered from incident moments instead of subtracting almost equal
reciprocals. These are algebraic zero constraints, not clamping rules.

The ordinary non-diagonal, non-rank-one `SetMassProperties` creation path keeps
its historical assertion diagnostics and **decomposition-failure** sphere fallback.
Only export opts into recoverable numerical checks; it never changes an installed
body or uses the sphere fallback to justify acceptance. For mass 2 and the exactly
represented SPD tensor `M = x x^T + 2^-8 Identity`, `x=(-4,-3,3,-3)`, creation can
succeed while the rounded export fails the unchanged eigenpair accuracy test.
The checked API now rejects that case with zero assertion callbacks.
Arbitrary rank-two/rank-three ill-conditioned creation, underflow-scale input,
and arbitrary inverse-setter states are **not newly certified** by this change.
In particular, successful initial creation is not a promise that checked export
must succeed.

## Permanent coverage and verification

`UnitTests/Physics/MassRoundTripTests.cpp` covers:

- the exact review fixture, its analytical `I^+ = P/4` world response and atomic
  checked-export rejection;
- exact rank-one particle directions, multiple scales and world orientations,
  with a velocity-Gram analytical pseudoinverse and unchanged exact nullity;
- `s x x^T`, `x=(0,-2,10,-10)`, for `s = 2^-40`, `1`, `2^40`,
  and every choice of structurally zero axis: direct decomposition and actual
  body creation have zero assertion callbacks, exact zero reconstruction entries,
  orthonormal completion, proper orientation and three exact zero planar inverses;
- mass-2 `x x^T + 2^-8 Identity`, `x=(-4,-3,3,-3)`: successful creation followed
  by atomic checked-export rejection, with zero assertion callbacks throughout;
- mass-2 `x x^T + 4 Identity`, `x=(1,2,3,4)`, created normally and then scaled by
  public `MotionProperties::ScaleToMass(ldexp(2.f,70))`: all stored mass/inertia
  inverses remain normal; checked export rejects atomically, raises no overflow,
  invalid-operation or divide-by-zero flags, and preserves the caller's trap
  masks. Linux and native Windows self-subprocesses enable real overflow traps
  and require a zero exit within a bounded wait (20 seconds);
- the same physical mass-scaling fixture at exponents `-100,0,40,50,54,70,100`:
  accepted exports recreate through the ordinary Body API with zero assertion
  callbacks and the unchanged `64*FLT_EPSILON` response budget; the larger
  unsupported non-diagonal tensors reject without flags;
- a Linux self-subprocess legacy death test requiring **SIGABRT** plus a pipe
  marker emitted only by the legacy API's specific failure trace. A sanitizer
  or other abort before that trace cannot pass. An assertion callback exits
  with a distinct non-signal status. Its helper case is skipped except when
  explicitly spawned by the parent; even `--no-skip` without the parent's
  marker cannot accidentally trigger the death. The wait is bounded;
- every coordinate null-plane choice of diagonal rank-two M;
- positive skinny diagonals with entries down to `2^-80`: creation keeps every
  rotational inverse positive; export either preserves them or explicitly
  rejects unrecoverable reciprocal cancellation;
- zero M and uniformly small positive M at scales `2^-40` and `2^-60`, which
  must retain the particle/reciprocal convention instead of sphere fallback;
- the existing full-rank particle tensors, mass scales, static/kinematic latent
  state, shape/COM offsets, DOF masks and sphere quadrature.

### Earlier structural-zero verification (before the overflow repair)

The isolated assertions-enabled GCC 13 Release build used `-ffp-contract=off`
and `-j2` at `/home/brkdrd/mover-rank-mass-build`. The two assertion-path
regressions were each executed RED before their respective fixes. Evidence and
final verification are in `review-artifacts/cycle2-*.log`:

- mass suite: **13 cases / 3,956 assertions pass**, including the self-subprocess
  SIGABRT check; the child-only helper is skipped by normal test selection;
- full suite: **540 of 541 executed cases pass**, with one child helper skipped.
  The sole failure remains the pre-existing
  `GyroscopicClampFailStopsEvenWithoutAnyDelivery` (12 assertions), independently
  reproduced by the original-object baseline already in `review-artifacts/`;
- HelloWorld passes;
- the independent 2,048-rank-one scan passes, with maximum orthogonality error
  `9.22872507e-8` and normalized reconstruction error `5.52172781e-8`;
- all 512 well-conditioned independent SPD exports pass; maximum world-response
  error remains `2.02206906e-6`;
- the independent conditioning scan exercises 1,664 successfully created SPD
  tensors with **zero export assertion callbacks**. Its other 2,432 inputs
  already assert during ordinary creation and are outside this export repair;
- an isolated native MSVC `/O2 /fp:precise`, assertions/traps-enabled standalone
  executes **2,060 rank-one decomposition/normal-creation cases**, with zero
  assertion callbacks; its standalone checked SPD decode also rejects with zero
  callbacks. The same targeted MassProperties/MotionProperties harness passes
  GCC ASAN+UBSAN, and a separate assertions-disabled GCC Release run still rejects
  the checked SPD decode. These are not native Body API or whole-engine sanitizer
  certification; no shared/main native build was rebuilt.

Native probe sources, command and log are
`C:/Users/goliy/AppData/Local/hermes/cache/scratch/rank-cycle2-native.*`.
The historical normal-creation eigenpair assertion on unsupported ill-conditioned
non-rank-one inputs remains intentional; for example the independent skinny
non-diagonal probe at exponent 16 still diagnoses insufficient accuracy.

### Overflow-trap repair verification

The exact public-Body regression was executed RED before the range guard:
Linux returned true while setting `FE_OVERFLOW` with exceptions masked, and
terminated with **SIGFPE** under overflow traps. A native full-Body replay
compiled against the original solver header (confirmed by `/showIncludes`)
likewise set the overflow flag and terminated with **0xc0000091** when trapped.
The repair does not mask that fault or accept an infinite comparison tolerance:
it rejects the unsupported range before entering the float diagnostics.

Verification reused the same isolated GCC build, verified its CMake source
path, and rebuilt with **`-j1`**. Results:

- registered mass suite: **16 cases / 4,259 assertions pass**, plus the
  overflow-trapped child (**21 assertions**) and the marked legacy death child;
- full suite: **543 of 544 executed cases pass**, two child-only cases skipped.
  The only failure remains `GyroscopicClampFailStopsEvenWithoutAnyDelivery`
  with its same 12 assertions; HelloWorld passes;
- an additional trapped public-Body sweep of **512 SPD fixtures over 10 mass
  scales** accepts 2,653 exports and rejects 2,467, with zero assertion callbacks
  or overflow/invalid/divide-by-zero flags. Accepted exports recreate normally;
  the largest measured relative response error is `7.56616691e-6`, below the
  unchanged `64*FLT_EPSILON` budget;
- the earlier independent rank-one/SPD/conditioning probes still pass: all
  512 well-conditioned exports remain accepted, and the 1,664 successfully
  created conditioning cases produce zero export assertion callbacks;
- a negative control replacing the legacy call with an unrelated `abort()`
  produces SIGABRT but **fails** the hardened death test for its missing pipe
  marker. The unchanged ordinary eigensolver still diagnoses the unsupported
  skinny non-diagonal exponent-16 input;
- native MSVC **`/O2 /fp:precise` full public Body API**: **15 cases / 4,234
  assertions pass**, plus the Windows trapped child (**19 assertions**).
  Only Body/MassProperties/MotionProperties and this test translation unit
  were rebuilt into scratch artifacts, linked against the existing native-
  precise library/framework with its matching defines (assertions disabled,
  cross-platform deterministic, FP exceptions enabled). The native runner's
  trap masks are restored exactly; no shared native build was modified.

Logs and reproducible standalone commands/scripts are under
`C:/Users/goliy/AppData/Local/hermes/cache/scratch/mass-overflow-*`.
The original independent `final-mass-overflow-review-{masked,traps}.log` now
also confirms false with **zero** `FE_OVERFLOW`, `FE_INVALID` and `FE_DIVBYZERO`.
This repair does not add whole-engine sanitizer or double-position certification;
those shared configurations are verified independently by the integrator.
