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
4. A temporary `MotionProperties` decodes the proposed settings using the real
   `SetMassProperties` path. The count of exactly zero planar inverses must be
   unchanged; every decoded inverse must remain finite and nonnegative.
5. Let A and B be the original and decoded **unmasked body-local 6x6 inverse
   response operators**, evaluated on all six coordinate bivectors using the
   actual float rotor/adjoint path. Export requires
   `||B - A||_F <= (64 * FLT_EPSILON) ||A||_F`. Differences and squared norms are
   accumulated in double. Checking the unmasked operator also protects latent
   inertia behind partial world-space DOF locks and static/kinematic motion types.

The last bound is an explicit *output accuracy budget*, `2^-17` relative
Frobenius error, not an eigenvalue cutoff or a claimed a-priori error theorem
for the eigensolver. Six bits of headroom over binary32 epsilon are reserved for
the matrix/rotor/adjoint encode-decode chain, but acceptance is decided by the
**measured decoded response**, not by that operation-count heuristic. It is a
normwise response contract, not a relative error guarantee for every impulse in
or arbitrarily close to a null space. No conditioning threshold turns a small
positive moment into zero. If cancellation or conditioning defeats these
checks, export rejects; it does not regularize M or fabricate sphere inertia.
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

After certification, a stable double Householder construction supplies the
occupied direction and an orthogonal completion, with eigenvalues
`[trace(M),0,0,0]`. The existing proper-orientation and reconstruction assertions
are retained unchanged. This fixes the review fixture's *initial creation*
assertion without weakening the generic eigensolver or its tolerances. It does
not claim to recover an ideal rank-one particle distribution from an already
rounded matrix whose exact rank is different.

`SetMassProperties` also keeps identity as the principal frame when M is
already exactly diagonal. Sorting/permuting it and converting that permutation
to a float rotor is unnecessary and perturbs coordinate-aligned null modes.
The former absolute `IsNearZero` gate is removed: a successful decomposition
of zero M keeps six exact zero inverses, and a successful small positive M
keeps its actual reciprocals rather than acquiring unit-sphere inertia.
For export, an exact zero `I_ab` plus PSD implies `d_a = d_b = 0`; other d's are
recovered from incident moments instead of subtracting almost equal
reciprocals. These are algebraic zero constraints, not clamping rules.

The generic non-diagonal, non-rank-one decomposition and its historical
**decomposition-failure** fallback remain unchanged. Arbitrary rank-two/rank-three
ill-conditioned creation, underflow-scale input, and arbitrary inverse-setter
states are **not newly certified** by this change. In particular, successful
initial creation is not a promise that checked export must succeed.

## Permanent coverage and verification

`UnitTests/Physics/MassRoundTripTests.cpp` covers:

- the exact review fixture, its analytical `I^+ = P/4` world response and atomic
  checked-export rejection;
- exact rank-one particle directions, multiple scales and world orientations,
  with a velocity-Gram analytical pseudoinverse and unchanged exact nullity;
- every coordinate null-plane choice of diagonal rank-two M;
- positive skinny diagonals with entries down to `2^-80`: creation keeps every
  rotational inverse positive; export either preserves them or explicitly
  rejects unrecoverable reciprocal cancellation;
- zero M and uniformly small positive M at scales `2^-40` and `2^-60`, which
  must retain the particle/reciprocal convention instead of sphere fallback;
- the existing full-rank particle tensors, mass scales, static/kinematic latent
  state, shape/COM offsets, DOF masks and sphere quadrature.

The isolated assertions-enabled GCC 13 Release build uses `-ffp-contract=off`
and `-j2` at `/home/brkdrd/mover-rank-mass-build`. RED evidence is
`/home/brkdrd/mover-rank-mass-red.log`; build/test and fail-stop subprocess logs
are in that build and its `review-artifacts/` directory. The final mass suite
passes **10 cases / 3,568 assertions**, and HelloWorld passes. The full suite
passes **537 of 538 cases**; its sole failure is the pre-existing
`GyroscopicClampFailStopsEvenWithoutAnyDelivery` (12 assertions). A separate
baseline executable relinked the original HEAD Body/MassProperties/
MotionProperties objects and reproduced that same failure, without this repair.
Native Windows and sanitizer verification are separate acceptance work, not
implied by WSL results.
