// Space-dynamics throughput benchmark for the 4D Jolt physics core.
//
// This is measurement code, not a production API. It builds a bounded, deterministic
// PhysicsSystem workload out of 4D tesseracts and reports the *full-run* latency
// distribution of PhysicsSystem::Update at a fixed 120 Hz timestep. Every step is
// measured, including the first one; there is no warmup phase, no discarded startup
// and no early exit. Sleeping is disabled so that the throughput figure describes
// the declared workload rather than the sleep heuristic.
//
// Two workloads are provided:
//
//   separated   N dynamic anisotropic tesseracts on a 4D grid, zero gravity, a
//               nonzero six-plane angular velocity each (gyroscopic term enabled),
//               and one common translational velocity. The common velocity keeps
//               the mutual separation of the bodies exactly constant, so this is a
//               pure free-rotation / integration / broadphase workload with zero
//               contacts. Each body's spin angular momentum and kinetic energy are
//               analytically constant; they are checked against an independent
//               double-precision 6x6 reference built from the inertia diagonal and
//               the body's rotation (no engine inertia/gyro code is used for it).
//
//   interacting N dynamic tesseracts in two lattice clusters closing on each other
//               in zero gravity, restitution 1, friction 0. Contacts occur mid-run.
//               Linear momentum must still be conserved (contact impulses are
//               internal); energy and angular momentum are reported, not gated,
//               because the contact/position solver is not an exact elastic map and
//               no accuracy budget for it has been established.
//
// The reported per-step cost includes the per-call timing overhead, which is
// measured and printed separately, and includes all first-call allocations inside
// the physics system. Numerical sampling is done between timed calls and its cost
// is measured and reported separately as well.
//
// Exit code is 0 only if every declared check passes AND every requested
// artifact (JSON/CSV) was written completely. 1 is a failed check, a non-finite
// reported value, or an artifact I/O failure; 2 is a command-line error. Invalid
// numbers are rejected before the workload runs: the step rate must be a finite
// positive number and the resulting float timestep and simulated duration must be
// positive and finite, because a run with rate=inf/NaN/-1 is not a defined
// measurement. JSON output is never emitted with nan/inf tokens.
//
// Build: see Tools/SpaceSimBenchmark/README.md (CMake subproject over Build/).

#include <Jolt/Jolt.h>

#include <Jolt/RegisterTypes.h>
#include <Jolt/Core/Factory.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Core/JobSystemSingleThreaded.h>
#include <Jolt/Core/JobSystemThreadPool.h>
#include <Jolt/Physics/PhysicsSettings.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/Physics/EPhysicsUpdateError.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Body/BodyLockInterface.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/ContactListener.h>
#include <Jolt/Physics/Collision/BroadPhase/BroadPhaseLayerInterfaceTable.h>
#include <Jolt/Physics/Collision/BroadPhase/ObjectVsBroadPhaseLayerFilterTable.h>
#include <Jolt/Physics/Collision/ObjectLayerPairFilterTable.h>

#include <atomic>
#include <cerrno>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "SpaceSimBenchmarkBuildInfo.h"

JPH_SUPPRESS_WARNINGS
using namespace JPH;

namespace
{
	// ---------------------------------------------------------------------
	// Manifest constants. These are part of the benchmark definition; they are
	// printed into the result file so a run can be reproduced exactly.
	// ---------------------------------------------------------------------
	constexpr float cDeltaTime = 1.0f / 120.0f;					// 120 Hz authoritative step
	constexpr double cDeadlineNs = 1.0e9 / 120.0;				// 8.333... ms
	constexpr uint cDefaultSteps = 1200;						// 10.0 s simulated, identical for every body count
	constexpr int cSampleInterval = 60;							// numerical sampling every 0.5 s, outside the timed call
	constexpr uint64 cSeed = 0x5EED4D0000000001ull;				// fixed seed, no randomness outside the harness

	constexpr ObjectLayer cLayerNonMoving = 0;
	constexpr ObjectLayer cLayerMoving = 1;
	constexpr uint cNumObjectLayers = 2;
	constexpr uint cNumBroadPhaseLayers = 2;

	// 1 : 1.75 : 2.5 : 3.25 anisotropy. Distinct principal inertias are what make
	// the gyroscopic term non-trivial; a hypercube would have six equal moments.
	constexpr float cHalfX = 0.40f, cHalfY = 0.70f, cHalfZ = 1.00f, cHalfW = 1.30f;
	constexpr float cGridSpacing = 8.0f;						// > 2 * bounding radius (3.66 m), so bodies never touch
	constexpr float cCommonVelocityX = 0.35f;					// identical translation for every body: separation is invariant
	constexpr float cCommonVelocityY = -0.20f;
	constexpr float cCommonVelocityZ = 0.15f;
	constexpr float cCommonVelocityW = 0.25f;

	constexpr float cClusterSpacing = 5.0f;					// > 2 * bounding radius, so no initial overlap
	constexpr float cClusterCenterOffset = 16.0f;			// interacting: clusters start at +-16 on X
	constexpr float cClusterClosingSpeed = 4.0f;			// each cluster moves inward at this speed

	// Declared prototype acceptance gates. They are chosen from float32 arithmetic
	// (6e-8 relative per operation, ~1e3 operations per body-step) with margin for
	// the second-order splitting error of the free-rotation integrator. They are
	// NOT production tolerances: no production accuracy budget exists in the
	// specification yet. Reported measurements are printed next to each gate.
	constexpr double cEnergyRelMax = 1.0e-4;					// separated: kinetic energy is constant analytically
	constexpr double cAngularMomentumRelMax = 1.0e-4;			// separated: spin L is constant analytically
	constexpr double cLinearMomentumRelMax = 1.0e-5;			// both: no external force acts on the system
	constexpr double cTotalAngularMomentumRelMax = 1.0e-3;		// interacting: report gate for L about the origin
	constexpr float cRotationNormalizedTol = 1.0e-4f;

	// ---------------------------------------------------------------------
	// Deterministic pseudo-random source (splitmix64). Only used at setup.
	// ---------------------------------------------------------------------
	struct Rng
	{
		uint64 mState;
		explicit Rng(uint64 inSeed) : mState(inSeed) {}
		uint64 Next()
		{
			uint64 z = (mState += 0x9E3779B97F4A7C15ull);
			z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
			z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
			return z ^ (z >> 31);
		}
		float Uniform(float inLo, float inHi)
		{
			return inLo + (inHi - inLo) * float(double(Next() >> 40) / double(1ull << 24));
		}
	};

	// ---------------------------------------------------------------------
	// Thread-safe counting contact listener. Contact callbacks run on job system
	// threads, so the counters are atomic.
	// ---------------------------------------------------------------------
	class CountingContactListener final : public ContactListener
	{
	public:
		std::atomic<uint64> mAdded { 0 };
		std::atomic<uint64> mPersisted { 0 };
		std::atomic<uint64> mRemoved { 0 };

		virtual void OnContactAdded(const Body &, const Body &, const ContactManifold &, ContactSettings &) override
		{
			mAdded.fetch_add(1, std::memory_order_relaxed);
		}
		virtual void OnContactPersisted(const Body &, const Body &, const ContactManifold &, ContactSettings &) override
		{
			mPersisted.fetch_add(1, std::memory_order_relaxed);
		}
		virtual void OnContactRemoved(const SubShapeIDPair &) override
		{
			mRemoved.fetch_add(1, std::memory_order_relaxed);
		}
	};

	// ---------------------------------------------------------------------
	// Independent double-precision reference for angular momentum and rotational
	// energy. It uses only public state (world rotation, world angular velocity,
	// principal inverse-inertia diagonal, inertia rotation) and the textbook
	// inertia operator expressed as 6x6 skew matrices. Nothing from
	// MotionProperties' gyroscopic integrator or inverse-inertia operator is used.
	// ---------------------------------------------------------------------
	struct RefVec
	{
		double v[6] = {};
	};

	struct RefMat
	{
		double a[4][4] = {};
	};

	constexpr int cI[6] = { 0, 0, 0, 1, 1, 2 };
	constexpr int cJ[6] = { 1, 2, 3, 2, 3, 3 };

	RefMat ToMat(const Mat44 &inM)
	{
		RefMat r;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				r.a[i][j] = double(inM(i, j));
		return r;
	}

	RefMat Transpose(const RefMat &inM)
	{
		RefMat r;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				r.a[i][j] = inM.a[j][i];
		return r;
	}

	RefMat Mul(const RefMat &inA, const RefMat &inB)
	{
		RefMat r;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				for (int k = 0; k < 4; ++k)
					r.a[i][j] += inA.a[i][k] * inB.a[k][j];
		return r;
	}

	// Bivector -> antisymmetric 4x4 matrix in the engine's plane ordering.
	RefMat Skew(const RefVec &inB)
	{
		RefMat r;
		for (int k = 0; k < 6; ++k)
		{
			r.a[cI[k]][cJ[k]] = inB.v[k];
			r.a[cJ[k]][cI[k]] = -inB.v[k];
		}
		return r;
	}

	RefVec Components(const RefMat &inM)
	{
		RefVec r;
		for (int k = 0; k < 6; ++k)
			r.v[k] = 0.5 * (inM.a[cI[k]][cJ[k]] - inM.a[cJ[k]][cI[k]]);
		return r;
	}

	RefVec FromBivec(BivecArg inB)
	{
		RefVec r;
		for (int k = 0; k < 6; ++k)
			r.v[k] = double(inB[k]);
		return r;
	}

	double Dot(const RefVec &inA, const RefVec &inB)
	{
		double r = 0;
		for (int k = 0; k < 6; ++k)
			r += inA.v[k] * inB.v[k];
		return r;
	}

	double Norm(const RefVec &inA)
	{
		return std::sqrt(Dot(inA, inA));
	}

	// Kinetic energy: translational 0.5 m |v|^2 plus rotational 0.5 Omega . L.
	// All four position dimensions are spatial, and the translational term is
	// formed from the double components: Vec4::LengthSq() would evaluate the sum
	// of squares in float and feed a rounded value into the "independent"
	// reference, which is exactly the arithmetic the reference is meant to avoid.
	double KineticEnergy(double inMass, const Vec4 &inV, const RefVec &inOmega, const RefVec &inSpinL)
	{
		const double vx = double(inV.GetX()), vy = double(inV.GetY());
		const double vz = double(inV.GetZ()), vw = double(inV.GetW());
		const double v_sq = vx * vx + vy * vy + vz * vz + vw * vw;
		return 0.5 * inMass * v_sq + 0.5 * Dot(inOmega, inSpinL);
	}

	// Spin angular momentum L = Q diag(I) Q^T W(Omega) Q ... Q^T, computed by
	// conjugating the bivector into the principal frame and back.
	RefVec ApplyInertia(const RefMat &inQ, const RefVec &inInertiaDiag, const RefVec &inOmega)
	{
		RefMat qt = Transpose(inQ);
		RefVec local = Components(Mul(qt, Mul(Skew(inOmega), inQ)));
		for (int k = 0; k < 6; ++k)
			local.v[k] *= inInertiaDiag.v[k];
		return Components(Mul(inQ, Mul(Skew(local), qt)));
	}

	double RotationOrthogonalityError(const RefMat &inQ)
	{
		double worst = 0;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
			{
				double s = 0;
				for (int k = 0; k < 4; ++k)
					s += inQ.a[k][i] * inQ.a[k][j];
				worst = std::max(worst, std::abs(s - (i == j? 1.0 : 0.0)));
			}
		return worst;
	}

	// ---------------------------------------------------------------------
	// One sampled instantaneous property set, plus the run-wide extrema.
	// ---------------------------------------------------------------------
	struct BodySampleState
	{
		double mEnergy = 0;			// 0.5 m |v|^2 + 0.5 Omega . L_spin
		double mSpinL[6] = {};		// I(Omega)
		double mTotalL[6] = {};		// I(Omega) + m (v ^ x)
	};

	struct RunChecks
	{
		bool mFinite = true;
		bool mRotationsNormalized = true;
		bool mAllActive = true;
		bool mNoVelocityClamp = true;
		bool mCoreNonFiniteVelocity = true;	// MotionProperties' sticky diagnostic
		bool mReferenceOrthogonal = true;

		double mEnergyRelMax = 0;				// max over samples and bodies: |E - E0| / E0
		double mTotalEnergyRelMax = 0;			// |total E - total E0| / sum |E0|
		double mSpinLRelMax = 0;				// max over samples and bodies: |L - L0| / |L0|
		double mTotalLRelMax = 0;				// |total L - total L0| / sum |L0|
		double mTotalPRelMax = 0;				// |total p - total p0| / sum m|v0|
		double mReferenceEScale = 0;			// sum |E0|
		double mReferenceLScale = 0;			// sum |L0|
		double mReferencePScale = 0;			// sum m|v0|

		double mMinPairGap = 0;					// separated only: min center distance - (r_i + r_j)
		bool mPairGapChecked = false;

		uint64 mNonFiniteObservations = 0;
		uint64 mNonNormalizedObservations = 0;

		bool mEnergyGated = false;
		bool mAngularMomentumGated = false;
		double mEnergyGate = 0;
		double mAngularMomentumGate = 0;
	};

	struct Sample
	{
		double mTotalP[4] = {};
		RefVec mTotalL;
		double mTotalEnergy = 0;
		double mEnergyAbsSum = 0;		// sum |E_i|: well-conditioned scale for relative energy drift
		double mPAbsSum = 0;			// sum m_i |v_i|: scale for relative momentum drift
		double mLAbsSum = 0;			// sum |L_i|: scale for relative angular momentum drift
	};

	// Defined with the other helpers below; declared here for the self-test.
	double Percentile(const std::vector<uint64> &inSorted, double inQ);

	// ---------------------------------------------------------------------
	// Harness self-test. Validates the independent reference math and the
	// percentile definition against hand-computable fixtures, so a PASS on the
	// workload means something about the harness as well as the engine.
	// ---------------------------------------------------------------------
	int SelfTest()
	{
		int failures = 0;
		auto check = [&failures](bool inOk, const char *inWhat)
		{
			std::printf("  %s %s\n", inOk? "ok  " : "FAIL", inWhat);
			if (!inOk) ++failures;
		};

		// 1. Skew/components round-trip.
		RefVec b;
		for (int k = 0; k < 6; ++k) b.v[k] = 0.3 + 0.7 * double(k);
		check(Norm(RefVec{{ Components(Skew(b)).v[0] - b.v[0], Components(Skew(b)).v[1] - b.v[1],
			Components(Skew(b)).v[2] - b.v[2], Components(Skew(b)).v[3] - b.v[3],
			Components(Skew(b)).v[4] - b.v[4], Components(Skew(b)).v[5] - b.v[5] }}) < 1e-14,
			"skew/components round-trip");

		// 2. Identity orientation: L_k = I_k * omega_k and E = 0.5 sum I omega^2.
		RefMat identity;
		for (int i = 0; i < 4; ++i) identity.a[i][i] = 1.0;
		RefVec inertia;
		for (int k = 0; k < 6; ++k) inertia.v[k] = 2.0 + 0.5 * double(k);
		RefVec omega;
		for (int k = 0; k < 6; ++k) omega.v[k] = 0.1 * double(k + 1) * (k % 2 == 0? 1.0 : -1.0);
		RefVec l = ApplyInertia(identity, inertia, omega);
		double worst = 0, e = 0, e_expected = 0;
		for (int k = 0; k < 6; ++k)
		{
			worst = std::max(worst, std::abs(l.v[k] - inertia.v[k] * omega.v[k]));
			e += 0.5 * omega.v[k] * l.v[k];
			e_expected += 0.5 * inertia.v[k] * omega.v[k] * omega.v[k];
		}
		check(worst < 1e-14, "identity-orientation angular momentum equals I * omega");
		check(std::abs(e - e_expected) < 1e-13, "identity-orientation energy equals 0.5 * sum(I * omega^2)");

		// 2b. Energy composition: 2 kg at (1,0,0,0) with no spin is 1 J, and all
		// four position dimensions are spatial: |(1,2,2,4)|^2 = 25, so m = 1
		// gives 12.5 J exactly in double arithmetic.
		RefVec zero6;
		check(std::abs(KineticEnergy(2.0, Vec4(1, 0, 0, 0), zero6, zero6) - 1.0) < 1e-15,
			"translational kinetic energy equals 0.5 * m * |v|^2");
		check(std::abs(KineticEnergy(1.0, Vec4(1, 2, 2, 4), zero6, zero6) - 12.5) < 1e-15,
			"translational energy includes the fourth spatial component");

		// 3. Rotation matrix from a known rotor is orthogonal with determinant +1.
		Mat44 m = Mat44::sRotation(Rotor::sExp(Bivec(0.21f, -0.37f, 0.11f, 0.53f, -0.19f, 0.07f)));
		RefMat q = ToMat(m);
		check(RotationOrthogonalityError(q) < 1e-6, "rotor -> matrix is orthogonal");
		double det = q.a[0][0] * (q.a[1][1] * (q.a[2][2] * q.a[3][3] - q.a[2][3] * q.a[3][2]) - q.a[1][2] * (q.a[2][1] * q.a[3][3] - q.a[2][3] * q.a[3][1]) + q.a[1][3] * (q.a[2][1] * q.a[3][2] - q.a[2][2] * q.a[3][1]))
			- q.a[0][1] * (q.a[1][0] * (q.a[2][2] * q.a[3][3] - q.a[2][3] * q.a[3][2]) - q.a[1][2] * (q.a[2][0] * q.a[3][3] - q.a[2][3] * q.a[3][0]) + q.a[1][3] * (q.a[2][0] * q.a[3][2] - q.a[2][2] * q.a[3][0]))
			+ q.a[0][2] * (q.a[1][0] * (q.a[2][1] * q.a[3][3] - q.a[2][3] * q.a[3][1]) - q.a[1][1] * (q.a[2][0] * q.a[3][3] - q.a[2][3] * q.a[3][0]) + q.a[1][3] * (q.a[2][0] * q.a[3][1] - q.a[2][1] * q.a[3][0]))
			- q.a[0][3] * (q.a[1][0] * (q.a[2][1] * q.a[3][2] - q.a[2][2] * q.a[3][1]) - q.a[1][1] * (q.a[2][0] * q.a[3][2] - q.a[2][2] * q.a[3][0]) + q.a[1][2] * (q.a[2][0] * q.a[3][1] - q.a[2][1] * q.a[3][0]));
		check(std::abs(det - 1.0) < 1e-6, "rotor -> matrix determinant is +1");

		// 4. Inertia operator is symmetric: A(x).y == A(y).x.
		RefVec x, y;
		for (int k = 0; k < 6; ++k) { x.v[k] = 0.17 * double(k + 1); y.v[k] = -0.11 * double(6 - k); }
		RefVec ax = ApplyInertia(q, inertia, x), ay = ApplyInertia(q, inertia, y);
		check(std::abs(Dot(ax, y) - Dot(ay, x)) < 1e-12, "rotated inertia operator is symmetric");
		check(Dot(ax, x) > 0, "rotated inertia operator is positive definite");

		// 5. Energy is rotation invariant when omega is carried along: E(Q) with
		// omega expressed in the rotated body frame equals E(identity). The matrix
		// comes from a float rotor, so the invariance holds to the orthogonality
		// error of Q (~1e-7 here), not to double-precision exactness.
		RefVec omega_rot = Components(Mul(q, Mul(Skew(omega), Transpose(q))));
		RefVec l_rot = ApplyInertia(q, inertia, omega_rot);
		double e_rot = 0;
		for (int k = 0; k < 6; ++k) e_rot += 0.5 * omega_rot.v[k] * l_rot.v[k];
		std::printf("       energy invariance: rotated=%.17g identity=%.17g rel=%.3e\n",
			e_rot, e_expected, std::abs(e_rot - e_expected) / std::abs(e_expected));
		check(std::abs(e_rot - e_expected) / std::abs(e_expected) < 1e-5, "energy invariant under a rigid reorientation");

		// 6. Percentile fixture: values 1..100, nearest-rank.
		std::vector<uint64> fixture;
		for (uint64 i = 1; i <= 100; ++i) fixture.push_back(i);
		check(Percentile(fixture, 0.50) == 50 && Percentile(fixture, 0.95) == 95 && Percentile(fixture, 0.99) == 99 && Percentile(fixture, 1.0) == 100,
			"nearest-rank percentiles on 1..100");

		std::printf("%s\n", failures == 0? "SELFTEST: PASS" : "SELFTEST: FAIL");
		return failures == 0? 0 : 1;
	}

	// ---------------------------------------------------------------------
	// Command line.
	// ---------------------------------------------------------------------
	enum class EScenario { Separated, Interacting };

	struct Args
	{
		EScenario mScenario = EScenario::Separated;
		uint mBodies = 128;
		uint mSteps = cDefaultSteps;
		int mThreads = 0;			// 0 = JobSystemSingleThreaded (declared single-threaded run)
		double mRate = 120.0;
		uint64 mSeed = cSeed;
		std::string mJson;
		std::string mCsv;
		bool mSelfTest = false;
	};

	void Usage()
	{
		std::printf(
			"SpaceSimBenchmark --scenario separated|interacting --bodies N [--steps 1200] [--threads 0]\n"
			"                  [--rate 120] [--seed 0x...] [--json out.json] [--csv steps.csv]\n"
			"                  [--selftest] [--help]\n"
			"\n"
			"  separated    dynamic anisotropic tesseracts on a 4D grid, zero gravity, six-plane\n"
			"               spin, gyroscopic term on, common translational velocity, no contacts.\n"
			"  interacting  two lattice clusters closing on each other in zero gravity; contacts.\n"
			"  bodies       integer in [2, 4096].\n"
			"  steps        measured steps, integer in [1, 10000000] (1200 = 10.0 s at 120 Hz).\n"
			"  threads      integer in [0, 64]; 0 = single-threaded job system, N>0 = N workers.\n"
			"  rate         finite positive step rate in Hz. The float timestep and the simulated\n"
			"               duration must both be positive and finite; otherwise the run is refused.\n"
			"  seed         integer, base 0 (0x... hexadecimal accepted).\n"
			"\n"
			"Exit codes: 0 all checks passed and every requested artifact was written;\n"
			"            1 failed check, non-finite report value, or artifact I/O failure;\n"
			"            2 invalid command line.\n");
	}

	// Strict numeric tokens: the whole string must be consumed, with no sign on
	// unsigned options, no whitespace, no trailing characters and no overflow.
	// The previous strtoul/strtod calls without an end pointer accepted "128x",
	// "-1" and "inf" (the last one silently produced a zero timestep).
	constexpr unsigned long long cMinBodies = 2;
	constexpr unsigned long long cMaxBodies = 4096;
	constexpr unsigned long long cMaxSteps = 10000000;
	constexpr unsigned long long cMaxThreads = 64;

	bool ParseUnsigned(const char *inText, unsigned long long inMax, unsigned long long &outValue)
	{
		if (inText == nullptr || inText[0] == '\0' || inText[0] == '-' || inText[0] == '+'
			|| std::isspace(static_cast<unsigned char>(inText[0])))
			return false;
		errno = 0;
		char *end = nullptr;
		unsigned long long v = std::strtoull(inText, &end, 10);
		if (errno != 0 || end == inText || *end != '\0' || v > inMax)
			return false;
		outValue = v;
		return true;
	}

	bool ParseSeed(const char *inText, uint64 &outValue)
	{
		if (inText == nullptr || inText[0] == '\0' || inText[0] == '-' || inText[0] == '+'
			|| std::isspace(static_cast<unsigned char>(inText[0])))
			return false;
		errno = 0;
		char *end = nullptr;
		unsigned long long v = std::strtoull(inText, &end, 0);
		if (errno != 0 || end == inText || *end != '\0')
			return false;
		outValue = uint64(v);
		return true;
	}

	bool ParseDoubleToken(const char *inText, double &outValue)
	{
		if (inText == nullptr || inText[0] == '\0' || std::isspace(static_cast<unsigned char>(inText[0])))
			return false;
		errno = 0;
		char *end = nullptr;
		double v = std::strtod(inText, &end);
		if (errno != 0 || end == inText || *end != '\0')
			return false;
		outValue = v;
		return true;
	}

	// Inspect without opening or creating artifacts. weakly_canonical alone leaves
	// dangling symlinks unresolved: a prior artifact write could create their target
	// and turn two apparently distinct paths into aliases. Reject such links, even
	// in parents/chains or when only one output is requested. This is CLI preflight,
	// not protection against concurrent adversarial filesystem changes.
	bool ResolveArtifactDestination(const std::string &inName, const char *inOption,
		std::filesystem::path &outPath, bool &outExists)
	{
		outExists = false;
		if (inName.empty())
			return true;
		std::error_code error;
		auto fail = [&](const char *inReason)
		{
			std::fprintf(stderr, "cannot resolve %s destination '%s': %s (%s)\n",
				inOption, inName.c_str(), inReason, error.message().c_str());
			return false;
		};
		outPath = std::filesystem::absolute(inName, error);
		if (error)
			return fail("absolute path failed");

		// Do not normalize '..' before inspection: a preceding symlink must be
		// resolved by the filesystem, not lexically discarded.
		std::filesystem::path prefix = outPath.root_path();
		for (const auto &component : outPath.relative_path())
		{
			prefix /= component;
			auto status = std::filesystem::symlink_status(prefix, error);
			if (error && !(error == std::errc::no_such_file_or_directory
				&& status.type() == std::filesystem::file_type::not_found))
				return fail("filesystem inspection failed");
			if (std::filesystem::is_symlink(status))
			{
				// canonical follows the entire chain and fails on missing targets,
				// loops, inaccessible ancestors and other resolution errors.
				(void)std::filesystem::canonical(prefix, error);
				if (error)
					return fail("symlink must resolve to an existing target");
			}
		}
		outPath = std::filesystem::weakly_canonical(outPath, error);
		if (error)
			return fail("canonical path failed");
		outExists = std::filesystem::exists(outPath, error);
		if (error)
			return fail("destination inspection failed");
		return true;
	}

	bool ParseArgs(int argc, char **argv, Args &outArgs)
	{
		for (int i = 1; i < argc; ++i)
		{
			std::string a = argv[i];
			auto next = [&](const char *inWhat) -> const char *
			{
				if (i + 1 >= argc)
				{
					std::fprintf(stderr, "missing value after %s\n", inWhat);
					std::exit(2);
				}
				return argv[++i];
			};
			if (a == "--scenario")
			{
				std::string s = next("--scenario");
				if (s == "separated") outArgs.mScenario = EScenario::Separated;
				else if (s == "interacting") outArgs.mScenario = EScenario::Interacting;
				else { std::fprintf(stderr, "unknown scenario '%s'\n", s.c_str()); return false; }
			}
			else if (a == "--bodies")
			{
				const char *tok = next("--bodies");
				unsigned long long v = 0;
				if (!ParseUnsigned(tok, cMaxBodies, v) || v < cMinBodies)
				{
					std::fprintf(stderr, "--bodies must be an integer in [%llu, %llu] (got '%s')\n", cMinBodies, cMaxBodies, tok);
					return false;
				}
				outArgs.mBodies = uint(v);
			}
			else if (a == "--steps")
			{
				const char *tok = next("--steps");
				unsigned long long v = 0;
				if (!ParseUnsigned(tok, cMaxSteps, v) || v < 1)
				{
					std::fprintf(stderr, "--steps must be an integer in [1, %llu] (got '%s')\n", cMaxSteps, tok);
					return false;
				}
				outArgs.mSteps = uint(v);
			}
			else if (a == "--threads")
			{
				const char *tok = next("--threads");
				unsigned long long v = 0;
				if (!ParseUnsigned(tok, cMaxThreads, v))
				{
					std::fprintf(stderr, "--threads must be an integer in [0, %llu] (got '%s')\n", cMaxThreads, tok);
					return false;
				}
				outArgs.mThreads = int(v);
			}
			else if (a == "--rate")
			{
				const char *tok = next("--rate");
				double v = 0;
				if (!ParseDoubleToken(tok, v) || !std::isfinite(v) || v <= 0.0)
				{
					std::fprintf(stderr, "--rate must be a finite positive number (got '%s')\n", tok);
					return false;
				}
				outArgs.mRate = v;
			}
			else if (a == "--seed")
			{
				const char *tok = next("--seed");
				uint64 v = 0;
				if (!ParseSeed(tok, v))
				{
					std::fprintf(stderr, "--seed must be an unsigned integer (base 0, 0x... accepted) (got '%s')\n", tok);
					return false;
				}
				outArgs.mSeed = v;
			}
			else if (a == "--json") outArgs.mJson = next("--json");
			else if (a == "--csv") outArgs.mCsv = next("--csv");
			else if (a == "--selftest") outArgs.mSelfTest = true;
			else if (a == "--help" || a == "-h") { Usage(); std::exit(0); }
			else { std::fprintf(stderr, "unknown argument '%s'\n", a.c_str()); return false; }
		}

		std::filesystem::path json_path, csv_path;
		bool json_exists = false, csv_exists = false;
		if (!ResolveArtifactDestination(outArgs.mJson, "--json", json_path, json_exists)
			|| !ResolveArtifactDestination(outArgs.mCsv, "--csv", csv_path, csv_exists))
			return false;
		if (!outArgs.mJson.empty() && !outArgs.mCsv.empty())
		{
			// Only existing files can be hard-link aliases. Do not ignore errors
			// from equivalent(): missing destinations were accounted for above.
			bool same_file = json_path == csv_path;
			if (!same_file && json_exists && csv_exists)
			{
				std::error_code error;
				same_file = std::filesystem::equivalent(json_path, csv_path, error);
				if (error)
				{
					std::fprintf(stderr, "cannot compare --json and --csv destinations: %s\n", error.message().c_str());
					return false;
				}
			}
			if (same_file)
			{
				std::fprintf(stderr, "--json and --csv must name different files (got '%s' and '%s')\n",
					outArgs.mJson.c_str(), outArgs.mCsv.c_str());
				return false;
			}
		}

		// A run is only defined if the fixed step is a representable positive
		// float and the reported simulated duration is positive and finite.
		// Without this, --rate inf would produce dt = 0, an infinite rate in the
		// manifest and a zero deadline, and still exit 0.
		const double nominal_dt = 1.0 / outArgs.mRate;
		const float dt = float(nominal_dt);
		if (!std::isfinite(nominal_dt) || !(dt > 0.0f) || !std::isfinite(dt))
		{
			std::fprintf(stderr, "--rate %.17g does not yield a positive finite float timestep\n", outArgs.mRate);
			return false;
		}
		const double simulated_seconds = double(dt) * double(outArgs.mSteps);
		if (!std::isfinite(simulated_seconds) || !(simulated_seconds > 0.0))
		{
			std::fprintf(stderr, "--rate %.17g with --steps %u does not yield a positive finite simulated duration\n",
				outArgs.mRate, outArgs.mSteps);
			return false;
		}
		return true;
	}

	// ---------------------------------------------------------------------
	// Helpers.
	// ---------------------------------------------------------------------
	uint64 NowNs()
	{
		using namespace std::chrono;
		return uint64(duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count());
	}

	double Percentile(const std::vector<uint64> &inSorted, double inQ)
	{
		// Nearest-rank percentile: index = ceil(q * n) - 1.
		if (inSorted.empty())
			return 0;
		size_t idx = size_t(std::ceil(inQ * double(inSorted.size())));
		if (idx == 0) idx = 1;
		if (idx > inSorted.size()) idx = inSorted.size();
		return double(inSorted[idx - 1]);
	}

	uint64 Fnv1a(uint64 inHash, const void *inData, size_t inSize)
	{
		const uint8 *p = static_cast<const uint8 *>(inData);
		for (size_t i = 0; i < inSize; ++i)
		{
			inHash ^= uint64(p[i]);
			inHash *= 0x100000001B3ull;
		}
		return inHash;
	}

	template <class T> uint64 HashPod(uint64 inHash, const T &inValue) { return Fnv1a(inHash, &inValue, sizeof(T)); }

	// ---------------------------------------------------------------------
	// Artifact writers. A benchmark that cannot write the file it was asked to
	// write did not succeed, so every failure (open, short write, flush, close)
	// is reported on stderr and propagated. /dev/full is the portable Linux way
	// to exercise the flush/close path: fopen succeeds and the error only shows
	// up when the buffer is drained.
	// ---------------------------------------------------------------------
	bool WriteWholeFile(const std::string &inPath, const std::string &inContent, const char *inWhat)
	{
		FILE *f = std::fopen(inPath.c_str(), "wb");
		if (f == nullptr)
		{
			std::fprintf(stderr, "%s: cannot open '%s': %s\n", inWhat, inPath.c_str(), std::strerror(errno));
			return false;
		}
		const size_t written = std::fwrite(inContent.data(), 1, inContent.size(), f);
		bool ok = written == inContent.size();
		if (!ok)
			std::fprintf(stderr, "%s: short write to '%s' (%zu of %zu bytes)\n", inWhat, inPath.c_str(), written, inContent.size());
		if (std::fflush(f) != 0)
		{
			std::fprintf(stderr, "%s: flush failed for '%s': %s\n", inWhat, inPath.c_str(), std::strerror(errno));
			ok = false;
		}
		if (std::fclose(f) != 0)
		{
			std::fprintf(stderr, "%s: close failed for '%s': %s\n", inWhat, inPath.c_str(), std::strerror(errno));
			ok = false;
		}
		return ok;
	}

	bool WriteCsvFile(const std::string &inPath, const std::vector<uint64> &inStepNs, const char *inWhat)
	{
		FILE *f = std::fopen(inPath.c_str(), "wb");
		if (f == nullptr)
		{
			std::fprintf(stderr, "%s: cannot open '%s': %s\n", inWhat, inPath.c_str(), std::strerror(errno));
			return false;
		}
		bool ok = std::fprintf(f, "step,ns\n") >= 0;
		for (size_t i = 0; ok && i < inStepNs.size(); ++i)
			ok = std::fprintf(f, "%zu,%llu\n", i + 1, (unsigned long long)inStepNs[i]) >= 0;
		if (!ok)
			std::fprintf(stderr, "%s: write failed for '%s'\n", inWhat, inPath.c_str());
		if (std::fflush(f) != 0)
		{
			std::fprintf(stderr, "%s: flush failed for '%s': %s\n", inWhat, inPath.c_str(), std::strerror(errno));
			ok = false;
		}
		if (std::fclose(f) != 0)
		{
			std::fprintf(stderr, "%s: close failed for '%s': %s\n", inWhat, inPath.c_str(), std::strerror(errno));
			ok = false;
		}
		return ok;
	}

	// ---------------------------------------------------------------------
	// Main body of the benchmark.
	// ---------------------------------------------------------------------
	int Run(const Args &inArgs)
	{
		// ParseArgs has already rejected rates whose float timestep or simulated
		// duration is not a positive finite number.
		const float dt = float(1.0 / inArgs.mRate);
		const double deadline_ns = 1.0e9 / inArgs.mRate;
		const double simulated_seconds = double(dt) * double(inArgs.mSteps);

		// ---- setup (timed separately) -------------------------------------
		uint64 setup_start = NowNs();

		TempAllocatorImpl temp_allocator(32 * 1024 * 1024);

		std::unique_ptr<JobSystem> job_system;
		if (inArgs.mThreads <= 0)
			job_system = std::make_unique<JobSystemSingleThreaded>(cMaxPhysicsJobs);
		else
			job_system = std::make_unique<JobSystemThreadPool>(cMaxPhysicsJobs, cMaxPhysicsBarriers, inArgs.mThreads);

		BroadPhaseLayerInterfaceTable bp_layer_interface(cNumObjectLayers, cNumBroadPhaseLayers);
		bp_layer_interface.MapObjectToBroadPhaseLayer(cLayerNonMoving, BroadPhaseLayer(0));
		bp_layer_interface.MapObjectToBroadPhaseLayer(cLayerMoving, BroadPhaseLayer(1));

		ObjectLayerPairFilterTable object_layer_pair_filter(cNumObjectLayers);
		object_layer_pair_filter.EnableCollision(cLayerMoving, cLayerMoving);

		ObjectVsBroadPhaseLayerFilterTable object_vs_bp_filter(bp_layer_interface, cNumBroadPhaseLayers, object_layer_pair_filter, cNumObjectLayers);

		const uint max_bodies = inArgs.mBodies + 16;
		const uint max_body_pairs = std::max<uint>(4096, inArgs.mBodies * 16);
		const uint max_contact_constraints = std::max<uint>(4096, inArgs.mBodies * 16);

		PhysicsSystem physics;
		physics.Init(max_bodies, 0, max_body_pairs, max_contact_constraints, bp_layer_interface, object_vs_bp_filter, object_layer_pair_filter);

		PhysicsSettings settings = physics.GetPhysicsSettings();
		physics.SetGravity(Vec4::sZero());
		settings.mAllowSleeping = false;			// throughput workload: no step may be skipped by sleeping
		physics.SetPhysicsSettings(settings);

		CountingContactListener contact_listener;
		physics.SetContactListener(&contact_listener);

		BodyInterface &bi = physics.GetBodyInterface();

		// One shared shape. Convex radius 0 keeps the tesseract exact (no sphere rounding).
		BoxShapeSettings shape_settings(Vec4(cHalfX, cHalfY, cHalfZ, cHalfW), 0.0f);
		shape_settings.SetEmbedded();
		ShapeSettings::ShapeResult shape_result = shape_settings.Create();
		if (shape_result.HasError())
		{
			std::fprintf(stderr, "shape creation failed: %s\n", shape_result.GetError().c_str());
			return 2;
		}
		ShapeRefC shape = shape_result.Get();

		// Body placement and initial state.
		struct BodyInit
		{
			Vec4 mPosition;
			Vec4 mVelocity;
			Bivec mAngularVelocity;
		};
		std::vector<BodyInit> inits;
		inits.reserve(inArgs.mBodies);

		Rng rng(inArgs.mSeed);

		auto make_spin = [&rng]() {
			// Six nonzero planes. Magnitude is bounded well below the body's angular
			// velocity clamp; the gyroscopic term is active because the three
			// principal moments differ.
			return Bivec(rng.Uniform(0.9f, 2.1f), rng.Uniform(-1.9f, -0.6f), rng.Uniform(0.5f, 1.8f),
						 rng.Uniform(-1.4f, 0.4f), rng.Uniform(0.7f, 2.0f), rng.Uniform(-1.7f, -0.2f));
		};

		if (inArgs.mScenario == EScenario::Separated)
		{
			// 4D grid with ceil(N^(1/4)) cells per axis; unused cells are skipped.
			const uint side = uint(std::ceil(std::pow(double(inArgs.mBodies), 0.25)));
			const float half = float(side - 1) * 0.5f * cGridSpacing;
			for (uint i0 = 0; i0 < side && inits.size() < inArgs.mBodies; ++i0)
				for (uint i1 = 0; i1 < side && inits.size() < inArgs.mBodies; ++i1)
					for (uint i2 = 0; i2 < side && inits.size() < inArgs.mBodies; ++i2)
						for (uint i3 = 0; i3 < side && inits.size() < inArgs.mBodies; ++i3)
						{
							BodyInit b;
							b.mPosition = Vec4(float(i0) * cGridSpacing - half, float(i1) * cGridSpacing - half,
											   float(i2) * cGridSpacing - half, float(i3) * cGridSpacing - half);
							b.mVelocity = Vec4(cCommonVelocityX, cCommonVelocityY, cCommonVelocityZ, cCommonVelocityW);
							b.mAngularVelocity = make_spin();
							inits.push_back(b);
						}
		}
		else
		{
			// Two lattice clusters, one at -X and one at +X, closing at cClusterClosingSpeed.
			const uint per_cluster = inArgs.mBodies / 2;
			const uint side = uint(std::ceil(std::pow(double(per_cluster), 0.25)));
			const float half = float(side - 1) * 0.5f * cClusterSpacing;
			for (uint cluster = 0; cluster < 2 && inits.size() < inArgs.mBodies; ++cluster)
			{
				const float sign = cluster == 0? -1.0f : 1.0f;
				for (uint i0 = 0; i0 < side && inits.size() < (cluster + 1) * per_cluster; ++i0)
					for (uint i1 = 0; i1 < side && inits.size() < (cluster + 1) * per_cluster; ++i1)
						for (uint i2 = 0; i2 < side && inits.size() < (cluster + 1) * per_cluster; ++i2)
							for (uint i3 = 0; i3 < side && inits.size() < (cluster + 1) * per_cluster; ++i3)
							{
								BodyInit b;
								b.mPosition = Vec4(sign * cClusterCenterOffset + float(i0) * cClusterSpacing - half,
												   float(i1) * cClusterSpacing - half,
												   float(i2) * cClusterSpacing - half,
												   float(i3) * cClusterSpacing - half);
								b.mVelocity = Vec4(-sign * cClusterClosingSpeed, 0.0f, 0.0f, 0.0f);
								b.mAngularVelocity = make_spin();
								inits.push_back(b);
							}
			}
			// Odd body counts: add one more on the -X side.
			while (inits.size() < inArgs.mBodies)
			{
				BodyInit b;
				b.mPosition = Vec4(-cClusterCenterOffset, 0.0f, 0.0f, 0.0f);
				b.mVelocity = Vec4(cClusterClosingSpeed, 0.0f, 0.0f, 0.0f);
				b.mAngularVelocity = make_spin();
				inits.push_back(b);
			}
		}

		std::vector<BodyID> ids;
		ids.reserve(inits.size());
		for (const BodyInit &init : inits)
		{
			BodyCreationSettings body_settings(shape, RVec4(init.mPosition.GetX(), init.mPosition.GetY(), init.mPosition.GetZ(), init.mPosition.GetW()),
				Rotor::sExp(init.mAngularVelocity * 0.25f), EMotionType::Dynamic, cLayerMoving);
			body_settings.mLinearVelocity = init.mVelocity;
			body_settings.mAngularVelocity = init.mAngularVelocity;
			body_settings.mApplyGyroscopicForce = true;		// gyro on, every body
			body_settings.mAllowSleeping = false;				// no step may be skipped by sleeping
			body_settings.mGravityFactor = 0.0f;
			body_settings.mLinearDamping = 0.0f;
			body_settings.mAngularDamping = 0.0f;
			body_settings.mMaxLinearVelocity = 1.0e5f;			// no clamping in this workload
			body_settings.mMaxAngularVelocity = 1.0e4f;
			body_settings.mFriction = inArgs.mScenario == EScenario::Interacting? 0.0f : 0.2f;
			body_settings.mRestitution = inArgs.mScenario == EScenario::Interacting? 1.0f : 0.0f;
			BodyID id = bi.CreateAndAddBody(body_settings, EActivation::Activate);
			if (id.IsInvalid())
			{
				std::fprintf(stderr, "body creation failed at %u\n", uint(ids.size()));
				return 2;
			}
			ids.push_back(id);
		}

		physics.OptimizeBroadPhase();

		const uint64 setup_ns = NowNs() - setup_start;

		// ---- initial numerical reference state -----------------------------
		const double bounding_radius = std::sqrt(double(cHalfX * cHalfX + cHalfY * cHalfY + cHalfZ * cHalfZ + cHalfW * cHalfW));

		struct InitialBody
		{
			double mEnergy;
			double mSpinL[6];
		};
		std::vector<InitialBody> initial(inArgs.mBodies);
		RunChecks checks;
		checks.mEnergyGated = inArgs.mScenario == EScenario::Separated;
		checks.mAngularMomentumGated = inArgs.mScenario == EScenario::Separated;
		checks.mEnergyGate = cEnergyRelMax;
		checks.mAngularMomentumGate = cAngularMomentumRelMax;

		auto accumulate_body = [&](const Body &inBody, BodySampleState &out)
		{
			const MotionProperties *mp = inBody.GetMotionProperties();
			RVec4 rp = inBody.GetCenterOfMassPosition();
			Vec4 p((float)rp.GetX(), (float)rp.GetY(), (float)rp.GetZ(), (float)rp.GetW());
			Vec4 v = inBody.GetLinearVelocity();
			Bivec w = inBody.GetAngularVelocity();
			float inv_mass = mp->GetInverseMass();
			double mass = 1.0 / double(inv_mass);

			Rotor principal = inBody.GetRotation() * mp->GetInertiaRotation();
			RefMat q = ToMat(Mat44::sRotation(principal));
			checks.mReferenceOrthogonal = checks.mReferenceOrthogonal && RotationOrthogonalityError(q) < 1.0e-5;

			RefVec omega = FromBivec(w);
			RefVec inertia_diag;
			for (int k = 0; k < 6; ++k)
			{
				float inv_i = mp->GetInverseInertiaDiagonal()[k];
				inertia_diag.v[k] = inv_i > 0.0f? 1.0 / double(inv_i) : 0.0;
			}

			RefVec spin_l = ApplyInertia(q, inertia_diag, omega);

			RefVec orbital;
			for (int k = 0; k < 6; ++k)
				orbital.v[k] = mass * (double(v[cI[k]]) * double(p[cJ[k]]) - double(v[cJ[k]]) * double(p[cI[k]]));

			out.mEnergy = KineticEnergy(mass, v, omega, spin_l);
			for (int k = 0; k < 6; ++k)
			{
				out.mSpinL[k] = spin_l.v[k];
				out.mTotalL[k] = spin_l.v[k] + orbital.v[k];
			}

			// Finiteness and rotation-normalization checks.
			bool finite = std::isfinite(p.GetX()) && std::isfinite(p.GetY()) && std::isfinite(p.GetZ()) && std::isfinite(p.GetW())
				&& std::isfinite(v.GetX()) && std::isfinite(v.GetY()) && std::isfinite(v.GetZ()) && std::isfinite(v.GetW())
				&& std::isfinite(w[0]) && std::isfinite(w[1]) && std::isfinite(w[2])
				&& std::isfinite(w[3]) && std::isfinite(w[4]) && std::isfinite(w[5])
				&& std::isfinite(out.mEnergy);
			if (!finite) { checks.mFinite = false; ++checks.mNonFiniteObservations; }
			if (!inBody.GetRotation().IsNormalized(cRotationNormalizedTol)) { checks.mRotationsNormalized = false; ++checks.mNonNormalizedObservations; }
			if (!inBody.IsActive()) checks.mAllActive = false;
			if (mp->WasLinearVelocityClamped() || mp->WasAngularVelocityClamped()) checks.mNoVelocityClamp = false;
			// The clamp flags above are sticky, but the core's own non-finite
			// velocity diagnostic was previously not consulted at all: a step that
			// produced nan/inf and was repaired before the next sample would have
			// gone unnoticed. It is sticky too, and PhysicsBridge treats it as a
			// numerical error, so the benchmark reports it as a check.
			if (mp->HadNonFiniteVelocity()) checks.mCoreNonFiniteVelocity = false;
		};

		// One pass over every body: totals plus (optionally) per-body drift against
		// the initial reference. Used for the initial state and every sample.
		auto sample_bodies = [&](bool inRecordInitial) -> Sample
		{
			Sample out;
			for (uint i = 0; i < inArgs.mBodies; ++i)
			{
				const Body *body = physics.GetBodyLockInterfaceNoLock().TryGetBody(ids[i]);
				if (body == nullptr)
				{
					checks.mFinite = false;
					++checks.mNonFiniteObservations;
					continue;
				}

				BodySampleState s;
				accumulate_body(*body, s);

				out.mTotalEnergy += s.mEnergy;
				out.mEnergyAbsSum += std::abs(s.mEnergy);
				{
					double l_sq = 0;
					for (int k = 0; k < 6; ++k)
						l_sq += s.mTotalL[k] * s.mTotalL[k];
					out.mLAbsSum += std::sqrt(l_sq);
				}
				for (int k = 0; k < 6; ++k)
					out.mTotalL.v[k] += s.mTotalL[k];

				const MotionProperties *mp = body->GetMotionProperties();
				Vec4 v = body->GetLinearVelocity();
				double m = 1.0 / double(mp->GetInverseMass());
				// Double components, not the float Vec4::LengthSq(), for the same
				// reason as KineticEnergy().
				const double vx = double(v.GetX()), vy = double(v.GetY());
				const double vz = double(v.GetZ()), vw = double(v.GetW());
				out.mPAbsSum += m * std::sqrt(vx * vx + vy * vy + vz * vz + vw * vw);
				out.mTotalP[0] += m * double(v.GetX());
				out.mTotalP[1] += m * double(v.GetY());
				out.mTotalP[2] += m * double(v.GetZ());
				out.mTotalP[3] += m * double(v.GetW());

				if (inRecordInitial)
				{
					initial[i].mEnergy = s.mEnergy;
					for (int k = 0; k < 6; ++k)
						initial[i].mSpinL[k] = s.mSpinL[k];
				}
				else
				{
					if (checks.mEnergyGated)
						checks.mEnergyRelMax = std::max(checks.mEnergyRelMax,
							std::abs(s.mEnergy - initial[i].mEnergy) / std::max(std::abs(initial[i].mEnergy), 1.0e-30));
					if (checks.mAngularMomentumGated)
					{
						double n0 = 0, dn = 0;
						for (int k = 0; k < 6; ++k)
						{
							double d = s.mSpinL[k] - initial[i].mSpinL[k];
							dn += d * d;
							n0 += initial[i].mSpinL[k] * initial[i].mSpinL[k];
						}
						if (n0 > 0)
							checks.mSpinLRelMax = std::max(checks.mSpinLRelMax, std::sqrt(dn) / std::sqrt(n0));
					}
				}
			}
			return out;
		};

		// Initial reference, taken from the created bodies before the first Update.
		Sample initial_sample = sample_bodies(true);
		checks.mReferenceEScale = initial_sample.mEnergyAbsSum;
		checks.mReferenceLScale = initial_sample.mLAbsSum;
		checks.mReferencePScale = initial_sample.mPAbsSum;

		// Separation check (separated workload only): exact at start and end using a
		// conservative bounding-sphere radius for each body.
		auto check_pair_gap = [&]()
		{
			double min_gap = std::numeric_limits<double>::max();
			for (uint i = 0; i < inArgs.mBodies; ++i)
			{
				const Body *a = physics.GetBodyLockInterfaceNoLock().TryGetBody(ids[i]);
				if (a == nullptr) continue;
				RVec4 pa = a->GetCenterOfMassPosition();
				for (uint j = i + 1; j < inArgs.mBodies; ++j)
				{
					const Body *b = physics.GetBodyLockInterfaceNoLock().TryGetBody(ids[j]);
					if (b == nullptr) continue;
					RVec4 pb = b->GetCenterOfMassPosition();
					double dx = double(pa.GetX() - pb.GetX()), dy = double(pa.GetY() - pb.GetY());
					double dz = double(pa.GetZ() - pb.GetZ()), dw = double(pa.GetW() - pb.GetW());
					double d = std::sqrt(dx * dx + dy * dy + dz * dz + dw * dw);
					min_gap = std::min(min_gap, d - 2.0 * bounding_radius);
				}
			}
			checks.mMinPairGap = checks.mPairGapChecked? std::min(checks.mMinPairGap, min_gap) : min_gap;
			checks.mPairGapChecked = true;
		};

		if (inArgs.mScenario == EScenario::Separated)
			check_pair_gap();

		// ---- timing overhead of the instrument itself ----------------------
		double instrument_overhead_ns = 0;
		{
			uint64 acc = 0;
			constexpr uint cIters = 200000;
			for (uint i = 0; i < cIters; ++i)
			{
				uint64 a = NowNs();
				uint64 b = NowNs();
				acc += b - a;
			}
			instrument_overhead_ns = double(acc) / double(cIters);
			std::printf("instrument: steady_clock pair cost %.1f ns (mean of %u)\n", instrument_overhead_ns, cIters);
		}

		// ---- measured run --------------------------------------------------
		std::vector<uint64> step_ns;
		step_ns.reserve(inArgs.mSteps);
		uint64 sampling_ns = 0;
		uint64 update_errors = 0;
		uint64 deadline_misses = 0;
		uint64 deadline_misses_stretch = 0;
		uint64 sample_count = 0;

		const uint64 run_start = NowNs();
		for (uint step = 0; step < inArgs.mSteps; ++step)
		{
			// The measured interval covers the whole Update call, including any
			// first-call allocation or lazy broadphase work. Nothing is excluded.
			uint64 t1 = NowNs();
			EPhysicsUpdateError err = physics.Update(dt, 1, &temp_allocator, job_system.get());
			uint64 t2 = NowNs();
			uint64 elapsed = t2 - t1;
			step_ns.push_back(elapsed);
			if (err != EPhysicsUpdateError::None)
				++update_errors;
			if (double(elapsed) > deadline_ns)
				++deadline_misses;
			if (double(elapsed) > 1.0e9 / 150.0)
				++deadline_misses_stretch;

			// Sampling happens between timed calls so it never inflates a step time.
			if ((step + 1) % cSampleInterval == 0 || step + 1 == inArgs.mSteps)
			{
				uint64 s1 = NowNs();
				Sample s = sample_bodies(false);
				RefVec dl;
				for (int k = 0; k < 6; ++k)
					dl.v[k] = s.mTotalL.v[k] - initial_sample.mTotalL.v[k];
				double dpn = 0;
				for (int k = 0; k < 4; ++k)
				{
					double d = s.mTotalP[k] - initial_sample.mTotalP[k];
					dpn += d * d;
				}
				if (checks.mReferencePScale > 0)
					checks.mTotalPRelMax = std::max(checks.mTotalPRelMax, std::sqrt(dpn) / checks.mReferencePScale);
				if (checks.mReferenceLScale > 0)
					checks.mTotalLRelMax = std::max(checks.mTotalLRelMax, Norm(dl) / checks.mReferenceLScale);
				if (checks.mReferenceEScale > 0)
					checks.mTotalEnergyRelMax = std::max(checks.mTotalEnergyRelMax,
						std::abs(s.mTotalEnergy - initial_sample.mTotalEnergy) / checks.mReferenceEScale);
				sampling_ns += NowNs() - s1;
				++sample_count;
			}
		}
		const uint64 run_ns = NowNs() - run_start;

		if (inArgs.mScenario == EScenario::Separated)
			check_pair_gap();

		// ---- state hash (determinism / reproducibility) --------------------
		uint64 hash = 0xCBF29CE484222325ull;
		for (const BodyID &id : ids)
		{
			const Body *body = physics.GetBodyLockInterfaceNoLock().TryGetBody(id);
			if (body == nullptr) continue;
			RVec4 p = body->GetCenterOfMassPosition();
			float pv[4] = { (float)p.GetX(), (float)p.GetY(), (float)p.GetZ(), (float)p.GetW() };
			Rotor r = body->GetRotation();
			float rv[8] = { r.GetScalar(), r.GetE12(), r.GetE13(), r.GetE14(), r.GetE23(), r.GetE24(), r.GetE34(), r.GetPseudoscalar() };
			Vec4 v = body->GetLinearVelocity();
			float vv[4] = { v.GetX(), v.GetY(), v.GetZ(), v.GetW() };
			Bivec w = body->GetAngularVelocity();
			float wv[6] = { w[0], w[1], w[2], w[3], w[4], w[5] };
			hash = Fnv1a(hash, pv, sizeof(pv));
			hash = Fnv1a(hash, rv, sizeof(rv));
			hash = Fnv1a(hash, vv, sizeof(vv));
			hash = Fnv1a(hash, wv, sizeof(wv));
		}

		// ---- statistics ----------------------------------------------------
		std::vector<uint64> sorted = step_ns;
		std::sort(sorted.begin(), sorted.end());
		uint64 total_step_ns = 0;
		for (uint64 v : step_ns) total_step_ns += v;

		const double p50 = Percentile(sorted, 0.50);
		const double p95 = Percentile(sorted, 0.95);
		const double p99 = Percentile(sorted, 0.99);
		const double pmax = double(sorted.back());
		const double pmin = double(sorted.front());
		const double mean = double(total_step_ns) / double(step_ns.size());
		const double steps_per_second = double(inArgs.mSteps) / (double(total_step_ns) * 1.0e-9);

		// ---- gate evaluation ----------------------------------------------
		struct Gate
		{
			std::string mName;
			double mValue;
			double mLimit;
			bool mGated;
			bool mPass;
			const char *mNote;
		};
		std::vector<Gate> gates;
		gates.push_back({ "finite_state", checks.mFinite? 0.0 : 1.0, 0.0, true, checks.mFinite, "1 = violation" });
		gates.push_back({ "reference_rotation_orthogonal", checks.mReferenceOrthogonal? 0.0 : 1.0, 0.0, true, checks.mReferenceOrthogonal, "harness reference: |Q^T Q - I| < 1e-5" });
		gates.push_back({ "rotation_normalized", double(checks.mNonNormalizedObservations), 0.0, true, checks.mRotationsNormalized, "non-normalized observations" });
		gates.push_back({ "all_bodies_active", checks.mAllActive? 0.0 : 1.0, 0.0, true, checks.mAllActive, "sleeping disabled" });
		gates.push_back({ "velocity_not_clamped", checks.mNoVelocityClamp? 0.0 : 1.0, 0.0, true, checks.mNoVelocityClamp, "1 = clamp seen" });
		gates.push_back({ "core_nonfinite_velocity", checks.mCoreNonFiniteVelocity? 0.0 : 1.0, 0.0, true, checks.mCoreNonFiniteVelocity, "MotionProperties sticky diagnostic" });
		gates.push_back({ "energy_rel_drift_per_body", checks.mEnergyRelMax, checks.mEnergyGate, checks.mEnergyGated, !checks.mEnergyGated || checks.mEnergyRelMax <= checks.mEnergyGate, "gyro free rotation" });
		gates.push_back({ "total_energy_rel_drift", checks.mTotalEnergyRelMax, checks.mEnergyGate, checks.mEnergyGated, !checks.mEnergyGated || checks.mTotalEnergyRelMax <= checks.mEnergyGate, "sum |E0| scale" });
		gates.push_back({ "spin_angular_momentum_rel_drift", checks.mSpinLRelMax, checks.mAngularMomentumGate, checks.mAngularMomentumGated, !checks.mAngularMomentumGated || checks.mSpinLRelMax <= checks.mAngularMomentumGate, "per body, separated" });
		gates.push_back({ "total_angular_momentum_rel_drift", checks.mTotalLRelMax, cTotalAngularMomentumRelMax, true, checks.mTotalLRelMax <= cTotalAngularMomentumRelMax, "L about origin" });
		gates.push_back({ "total_linear_momentum_rel_drift", checks.mTotalPRelMax, cLinearMomentumRelMax, true, checks.mTotalPRelMax <= cLinearMomentumRelMax, "no external force" });
		if (inArgs.mScenario == EScenario::Separated)
			gates.push_back({ "pair_separation", checks.mMinPairGap, 0.0, true, checks.mMinPairGap > 0.0, "min gap (m), must stay positive" });
		else
			gates.push_back({ "contacts_observed", double(contact_listener.mAdded.load()), 1.0, true, contact_listener.mAdded.load() > 0, "interacting workload must actually interact" });

		bool gate_failed = false;
		for (const Gate &g : gates)
			if (g.mGated && !g.mPass)
				gate_failed = true;
		if (update_errors > 0)
			gate_failed = true;

		// ---- reported-value finiteness -------------------------------------
		// A NaN/Inf in any reported number makes the measurement invalid, and it
		// must never reach the JSON file (nan/inf are not JSON tokens). Every
		// double that is serialized is checked here; a violation fails the run
		// and suppresses the JSON artifact instead of emitting invalid JSON.
		bool report_finite = true;
		auto finite_value = [&](const char *inName, double inValue)
		{
			if (!std::isfinite(inValue))
			{
				report_finite = false;
				std::fprintf(stderr, "non-finite reported value: %s\n", inName);
			}
		};
		finite_value("rate_hz", inArgs.mRate);
		finite_value("timestep_s", double(dt));
		finite_value("simulated_seconds", simulated_seconds);
		finite_value("deadline_ns", deadline_ns);
		finite_value("instrument_overhead_ns", instrument_overhead_ns);
		finite_value("step_ns.min", pmin);
		finite_value("step_ns.p50", p50);
		finite_value("step_ns.p95", p95);
		finite_value("step_ns.p99", p99);
		finite_value("step_ns.max", pmax);
		finite_value("step_ns.mean", mean);
		finite_value("throughput_steps_per_second", steps_per_second);
		finite_value("energy_rel_drift", checks.mEnergyRelMax);
		finite_value("total_energy_rel_drift", checks.mTotalEnergyRelMax);
		finite_value("spin_angular_momentum_rel_drift", checks.mSpinLRelMax);
		finite_value("total_angular_momentum_rel_drift", checks.mTotalLRelMax);
		finite_value("total_linear_momentum_rel_drift", checks.mTotalPRelMax);
		finite_value("min_pair_gap_m", checks.mMinPairGap);
		finite_value("reference_energy_scale", checks.mReferenceEScale);
		finite_value("reference_angular_momentum_scale", checks.mReferenceLScale);
		finite_value("reference_linear_momentum_scale", checks.mReferencePScale);
		for (const Gate &g : gates)
		{
			finite_value(g.mName.c_str(), g.mValue);
			finite_value(g.mName.c_str(), g.mLimit);
		}
		if (!report_finite)
			gate_failed = true;

		// ---- report --------------------------------------------------------
		std::printf("scenario=%s bodies=%u steps=%u rate=%.1f Hz dt=%.9f s threads=%d\n",
			inArgs.mScenario == EScenario::Separated? "separated" : "interacting", inArgs.mBodies, inArgs.mSteps, inArgs.mRate, double(dt), inArgs.mThreads);
		std::printf("setup: %.3f ms\n", double(setup_ns) * 1.0e-6);
		std::printf("run:   %.3f ms total, %.3f ms sum(steps), sampling %.3f ms\n",
			double(run_ns) * 1.0e-6, double(total_step_ns) * 1.0e-6, double(sampling_ns) * 1.0e-6);
		std::printf("step ns: min=%.0f p50=%.0f p95=%.0f p99=%.0f max=%.0f mean=%.0f\n", pmin, p50, p95, p99, pmax, mean);
		std::printf("throughput: %.1f steps/s; deadline misses %.3f ms: %llu, %.3f ms: %llu\n",
			steps_per_second, deadline_ns * 1.0e-6, (unsigned long long)deadline_misses, 1000.0 / 150.0, (unsigned long long)deadline_misses_stretch);
		std::printf("checks: energy_rel=%.3e total_energy_rel=%.3e spin_L_rel=%.3e total_L_rel=%.3e total_p_rel=%.3e contacts(added/persisted/removed)=%llu/%llu/%llu\n",
			checks.mEnergyRelMax, checks.mTotalEnergyRelMax, checks.mSpinLRelMax, checks.mTotalLRelMax, checks.mTotalPRelMax,
			(unsigned long long)contact_listener.mAdded.load(), (unsigned long long)contact_listener.mPersisted.load(), (unsigned long long)contact_listener.mRemoved.load());
		std::printf("reference scales: sum|E0|=%.6e sum|L0|=%.6e sum m|v0|=%.6e\n",
			checks.mReferenceEScale, checks.mReferenceLScale, checks.mReferencePScale);
		if (checks.mPairGapChecked)
			std::printf("min pair gap: %.6f m\n", checks.mMinPairGap);
		std::printf("state hash: %016llx\n", (unsigned long long)hash);

		for (const Gate &g : gates)
			std::printf("  %-4s %-36s value=%.6e limit=%.6e %s\n", (!g.mGated? "info" : (g.mPass? "ok" : "FAIL")), g.mName.c_str(), g.mValue, g.mLimit, g.mNote);

		// ---- artifacts ------------------------------------------------------
		// Writing is part of the measurement contract: a run that cannot produce
		// the requested file exits non-zero and never prints RESULT: PASS.
		bool artifact_ok = true;
		if (!inArgs.mCsv.empty() && !WriteCsvFile(inArgs.mCsv, step_ns, "CSV"))
			artifact_ok = false;

		if (!inArgs.mJson.empty())
		{
			if (!report_finite)
			{
				std::fprintf(stderr, "JSON: refusing to write '%s' with non-finite values\n", inArgs.mJson.c_str());
				artifact_ok = false;
			}
			else
			{
				std::string j = "{\n";
				char buf[8192];
				bool format_ok = true;
				auto append = [&](const char *inFormat, auto &&... inArgsFmt)
				{
					int n = std::snprintf(buf, sizeof(buf), inFormat, inArgsFmt...);
					if (n < 0 || size_t(n) >= sizeof(buf))
						format_ok = false;
					else
						j.append(buf, size_t(n));
				};

				j += "  \"manifest\": {\n";
				append(
					"    \"scenario\": \"%s\", \"bodies\": %u, \"steps\": %u, \"timestep_s\": %.17g, \"rate_hz\": %.17g,\n"
					"    \"simulated_seconds\": %.17g, \"worker_threads\": %d, \"job_system\": \"%s\",\n"
					"    \"seed\": \"0x%016llx\", \"gyroscopic_term\": true, \"sleeping\": false, \"gravity\": [0,0,0,0],\n"
					"    \"shape_half_extents\": [%.4f, %.4f, %.4f, %.4f], \"linear_damping\": 0.0, \"angular_damping\": 0.0,\n"
					"    \"solver_velocity_steps\": %u, \"solver_position_steps\": %u,\n"
					"    \"sample_interval_steps\": %d, \"percentile_method\": \"nearest-rank\",\n"
					"    \"compiled_config\": \"%s\",\n",
					inArgs.mScenario == EScenario::Separated? "separated" : "interacting", inArgs.mBodies, inArgs.mSteps,
					double(dt), inArgs.mRate, simulated_seconds, inArgs.mThreads,
					inArgs.mThreads <= 0? "JobSystemSingleThreaded" : "JobSystemThreadPool",
					(unsigned long long)inArgs.mSeed, cHalfX, cHalfY, cHalfZ, cHalfW,
					settings.mNumVelocitySteps, settings.mNumPositionSteps, cSampleInterval,
#ifdef JPH_DOUBLE_PRECISION
					"JPH_DOUBLE_PRECISION"
#else
					"float"
#endif
				);
				// Provenance of the binary that produced the numbers. The library
				// and the harness TU are compiled with different flags, so both are
				// recorded; the source fingerprints cover the compiled harness
				// inputs and the Jolt sources linked into the binary.
				append(
					"    \"compiler\": \"%s\", \"build_type\": \"%s\",\n"
					"    \"library_lto\": \"%s\", \"benchmark_lto\": \"%s\", \"cross_platform_deterministic\": \"%s\", \"fast_math\": \"%s\",\n"
					"    \"harness_source_sha256\": \"%s\", \"core_source_sha256\": \"%s\", \"build_sha256\": \"%s\",\n"
					"    \"core_source_files\": %d, \"git_head_at_configure\": \"%s\", \"git_dirty_entries_at_configure\": %d,\n"
					"    \"compile_definitions\": \"%s\",\n"
					"    \"benchmark_compile_flags\": \"%s\",\n"
					"    \"library_compile_flags\": \"%s\",\n"
					"    \"cmake_cache\": \"%s\"\n",
					SPACESIMBENCH_COMPILER, SPACESIMBENCH_BUILD_TYPE,
					SPACESIMBENCH_LIB_LTO, SPACESIMBENCH_BENCH_LTO, SPACESIMBENCH_CROSS_PLATFORM_DETERMINISTIC, SPACESIMBENCH_FAST_MATH,
					SPACESIMBENCH_SOURCE_SHA256, SPACESIMBENCH_CORE_SHA256, SPACESIMBENCH_BUILD_SHA256,
					SPACESIMBENCH_CORE_SOURCE_FILES, SPACESIMBENCH_GIT_HEAD, SPACESIMBENCH_GIT_DIRTY_ENTRIES,
					SPACESIMBENCH_DEFINITIONS, SPACESIMBENCH_BENCH_COMPILE_FLAGS, SPACESIMBENCH_LIB_COMPILE_FLAGS,
					SPACESIMBENCH_CACHE
				);
				j += "  },\n  \"timing\": {\n";
				append(
					"    \"setup_ns\": %llu, \"run_wall_ns\": %llu, \"sum_step_ns\": %llu, \"sampling_ns\": %llu,\n"
					"    \"instrument_overhead_ns\": %.2f,\n"
					"    \"step_ns\": { \"min\": %.0f, \"p50\": %.0f, \"p95\": %.0f, \"p99\": %.0f, \"max\": %.0f, \"mean\": %.0f },\n"
					"    \"throughput_steps_per_second\": %.3f,\n"
					"    \"deadline_ns_8333\": %.3f, \"deadline_misses\": %llu, \"deadline_misses_150hz\": %llu\n",
					(unsigned long long)setup_ns, (unsigned long long)run_ns, (unsigned long long)total_step_ns, (unsigned long long)sampling_ns,
					instrument_overhead_ns, pmin, p50, p95, p99, pmax, mean, steps_per_second, deadline_ns,
					(unsigned long long)deadline_misses, (unsigned long long)deadline_misses_stretch);
				j += "  },\n  \"checks\": {\n";
				append(
					"    \"finite\": %s, \"rotations_normalized\": %s, \"all_bodies_active\": %s, \"velocity_not_clamped\": %s, \"core_nonfinite_velocity\": %s,\n"
					"    \"non_finite_observations\": %llu, \"non_normalized_observations\": %llu,\n"
					"    \"energy_rel_drift\": %.6e, \"total_energy_rel_drift\": %.6e, \"energy_gated\": %s,\n"
					"    \"spin_angular_momentum_rel_drift\": %.6e, \"spin_angular_momentum_gated\": %s,\n"
					"    \"total_angular_momentum_rel_drift\": %.6e, \"total_linear_momentum_rel_drift\": %.6e,\n"
					"    \"min_pair_gap_m\": %.6f, \"pair_gap_checked\": %s,\n"
					"    \"reference_energy_scale\": %.6e, \"reference_angular_momentum_scale\": %.6e, \"reference_linear_momentum_scale\": %.6e,\n"
					"    \"contacts_added\": %llu, \"contacts_persisted\": %llu, \"contacts_removed\": %llu,\n"
					"    \"update_errors\": %llu, \"samples\": %llu\n",
					checks.mFinite? "true" : "false", checks.mRotationsNormalized? "true" : "false",
					checks.mAllActive? "true" : "false", checks.mNoVelocityClamp? "true" : "false",
					checks.mCoreNonFiniteVelocity? "true" : "false",
					(unsigned long long)checks.mNonFiniteObservations, (unsigned long long)checks.mNonNormalizedObservations,
					checks.mEnergyRelMax, checks.mTotalEnergyRelMax, checks.mEnergyGated? "true" : "false",
					checks.mSpinLRelMax, checks.mAngularMomentumGated? "true" : "false",
					checks.mTotalLRelMax, checks.mTotalPRelMax,
					checks.mMinPairGap, checks.mPairGapChecked? "true" : "false",
					checks.mReferenceEScale, checks.mReferenceLScale, checks.mReferencePScale,
					(unsigned long long)contact_listener.mAdded.load(), (unsigned long long)contact_listener.mPersisted.load(),
					(unsigned long long)contact_listener.mRemoved.load(), (unsigned long long)update_errors, (unsigned long long)sample_count);
				j += "  },\n  \"gates\": [\n";
				for (size_t i = 0; i < gates.size(); ++i)
				{
					const Gate &g = gates[i];
					append(
						"    { \"name\": \"%s\", \"value\": %.6e, \"limit\": %.6e, \"gated\": %s, \"pass\": %s, \"note\": \"%s\" }%s\n",
						g.mName.c_str(), g.mValue, g.mLimit, g.mGated? "true" : "false", g.mPass? "true" : "false", g.mNote,
						i + 1 == gates.size()? "" : ",");
				}
				append("  ],\n  \"state_hash\": \"%016llx\",\n  \"result\": \"%s\"\n}\n",
					(unsigned long long)hash, (gate_failed || !artifact_ok)? "FAIL" : "PASS");

				if (!format_ok)
				{
					std::fprintf(stderr, "JSON: formatting failed for '%s'\n", inArgs.mJson.c_str());
					artifact_ok = false;
				}
				else if (!WriteWholeFile(inArgs.mJson, j, "JSON"))
					artifact_ok = false;
			}
		}

		const bool failed = gate_failed || !artifact_ok;
		std::printf("%s\n", failed? "RESULT: FAIL" : "RESULT: PASS");
		return failed? 1 : 0;
	}
}

int main(int argc, char **argv)
{
	Args args;
	if (!ParseArgs(argc, argv, args))
	{
		Usage();
		return 2;
	}

	RegisterDefaultAllocator();
	Factory::sInstance = new Factory();
	RegisterTypes();

	int result = args.mSelfTest? SelfTest() : Run(args);

	UnregisterTypes();
	delete Factory::sInstance;
	Factory::sInstance = nullptr;
	return result;
}
