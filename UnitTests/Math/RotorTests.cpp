// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "Math/Cl4Reference.h"
#include <Jolt/Math/Rotor.h>
#include <Jolt/Math/Mat44.h>

JPH_SUPPRESS_WARNINGS_STD_BEGIN
#include <random>
JPH_SUPPRESS_WARNINGS_STD_END

// Test-first suite for the Cl(4,0) rotor (CLAUDE.md Phase A step 2).
// Everything is validated against the brute force reference algebra in Cl4Reference.h.
// Tests marked "EXPECTED TO FAIL" pin down bugs from the CLAUDE.md P0 list; they must pass
// once the corresponding fix lands and must NOT be weakened to make the suite green.
TEST_SUITE("RotorTests")
{
	// Random unit vector
	static Vec4 sRandomUnit(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		for (;;)
		{
			Vec4 v(dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom));
			if (v.Length() > 0.1f)
				return Vec4(v.Normalized());
		}
	}

	// Random orthonormal basis of R^4 via Gram-Schmidt
	static void sRandomOrthoBasis(UnitTestRandom &ioRandom, Vec4 outBasis[4])
	{
		for (;;)
		{
			bool valid = true;
			for (int i = 0; i < 4 && valid; ++i)
			{
				Vec4 v = sRandomUnit(ioRandom);
				for (int j = 0; j < i; ++j)
					v = Vec4(v - outBasis[j] * v.Dot(outBasis[j]));
				if (v.Length() < 0.1f)
					valid = false;
				else
					outBasis[i] = Vec4(v.Normalized());
			}
			if (valid)
				return;
		}
	}

	// Random valid rotor (an actual element of Spin(4)): product of a double rotation and a simple rotation
	static Rotor sRandomRotor(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> angle(-JPH_PI, JPH_PI);
		Vec4 b[4];
		sRandomOrthoBasis(ioRandom, b);
		Rotor r = Rotor::sDoubleRotation(b[0], b[1], angle(ioRandom), b[2], b[3], angle(ioRandom));
		Vec4 b2[4];
		sRandomOrthoBasis(ioRandom, b2);
		return r * Rotor::sRotation(b2[0], b2[1], angle(ioRandom));
	}

	// Check a Rotor against a reference multivector (which must be even grade)
	static void sCheckEqualsRef(RotorArg inR, const Cl4Ref::MV &inRef, double inTolerance = 1.0e-5)
	{
		CHECK(Cl4Ref::FromRotor(inR).MaxAbsDiff(inRef) < inTolerance);
	}

	TEST_CASE("TestRotorIdentity")
	{
		Rotor identity = Rotor::sIdentity();
		CHECK(identity.GetScalar() == 1.0f);
		CHECK(identity.GetPseudoscalar() == 0.0f);
		CHECK(identity * Vec4(1, 2, 3, 4) == Vec4(1, 2, 3, 4));
		CHECK(Cl4Ref::RotorNormScalar(identity) == 1.0);
		CHECK(Cl4Ref::RotorNormDefect(identity) == 0.0);
	}

	// Validates the hand derived 8x8 geometric product table (see project_4d_rotor_formulas memory /
	// Rotor.inl) against the reference algebra, one basis element pair at a time
	TEST_CASE("TestRotorBasisProductTable")
	{
		for (int i = 0; i < 8; ++i)
			for (int j = 0; j < 8; ++j)
			{
				float ac[8] = { }, bc[8] = { };
				ac[i] = 1.0f;
				bc[j] = 1.0f;
				Rotor a(ac[0], ac[1], ac[2], ac[3], ac[4], ac[5], ac[6], ac[7]);
				Rotor b(bc[0], bc[1], bc[2], bc[3], bc[4], bc[5], bc[6], bc[7]);
				Cl4Ref::MV ref = Cl4Ref::Mul(Cl4Ref::FromRotor(a), Cl4Ref::FromRotor(b));
				sCheckEqualsRef(a * b, ref, 1.0e-12);
			}
	}

	// The product formula is bilinear, so testing on random dense even elements catches any
	// wrong-sign combination the basis test might hide through cancellation
	TEST_CASE("TestRotorGeometricProductRandom")
	{
		UnitTestRandom random(12345);
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		for (int n = 0; n < 100; ++n)
		{
			Rotor a(dist(random), dist(random), dist(random), dist(random), dist(random), dist(random), dist(random), dist(random));
			Rotor b(dist(random), dist(random), dist(random), dist(random), dist(random), dist(random), dist(random), dist(random));
			Cl4Ref::MV ref = Cl4Ref::Mul(Cl4Ref::FromRotor(a), Cl4Ref::FromRotor(b));
			sCheckEqualsRef(a * b, ref);
		}
	}

	TEST_CASE("TestRotorReversed")
	{
		UnitTestRandom random(23456);
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		for (int n = 0; n < 10; ++n)
		{
			Rotor r(dist(random), dist(random), dist(random), dist(random), dist(random), dist(random), dist(random), dist(random));
			sCheckEqualsRef(r.Reversed(), Cl4Ref::Reverse(Cl4Ref::FromRotor(r)), 1.0e-12);
		}
	}

	TEST_CASE("TestRotorSRotation")
	{
		UnitTestRandom random(34567);
		uniform_real_distribution<float> angle_dist(-JPH_PI, JPH_PI);
		for (int n = 0; n < 20; ++n)
		{
			Vec4 b[4];
			sRandomOrthoBasis(random, b);
			float angle = angle_dist(random);
			Rotor r = Rotor::sRotation(b[0], b[1], angle);

			// A simple rotation has no pseudoscalar part and is a valid unit rotor
			CHECK_APPROX_EQUAL(r.GetPseudoscalar(), 0.0f);
			CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormScalar(r)), 1.0f, 1.0e-5f);
			CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormDefect(r)), 0.0f, 1.0e-5f);

			// Vectors in the rotation plane stay in the plane and rotate by the given angle
			Vec4 rotated = r * b[0];
			CHECK_APPROX_EQUAL(rotated.Length(), 1.0f, 1.0e-5f);
			CHECK_APPROX_EQUAL(rotated.Dot(b[0]), cos(angle), 1.0e-5f);
			CHECK_APPROX_EQUAL(rotated.Dot(b[2]), 0.0f, 1.0e-5f);
			CHECK_APPROX_EQUAL(rotated.Dot(b[3]), 0.0f, 1.0e-5f);

			// Vectors orthogonal to the rotation plane are unchanged
			CHECK_APPROX_EQUAL(r * b[2], b[2], 1.0e-5f);
			CHECK_APPROX_EQUAL(r * b[3], b[3], 1.0e-5f);

			// Rotating in the same plane composes angles; reverse is the inverse rotation
			Rotor r2 = Rotor::sRotation(b[0], b[1], 0.5f * angle);
			CHECK((r2 * r2).IsClose(r, 1.0e-8f));
			CHECK_APPROX_EQUAL(r.Reversed() * rotated, b[0], 1.0e-5f);
		}
	}

	TEST_CASE("TestRotorSandwichVsReference")
	{
		UnitTestRandom random(45678);
		for (int n = 0; n < 50; ++n)
		{
			Rotor r = sRandomRotor(random);
			Vec4 v = sRandomUnit(random);
			Cl4Ref::MV ref = Cl4Ref::Sandwich(Cl4Ref::FromRotor(r), Cl4Ref::FromVec4(v));

			// For a valid rotor the sandwich product has no grade 3 remainder ...
			CHECK(Cl4Ref::GradePart(ref, 3).MaxAbs() < 1.0e-5);

			// ... and the grade 1 part is what Rotor::operator* must return
			CHECK(Cl4Ref::Grade1Diff(ref, r * v) < 1.0e-4);
		}
	}

	TEST_CASE("TestRotorSandwichPreservesGeometry")
	{
		UnitTestRandom random(56789);
		for (int n = 0; n < 20; ++n)
		{
			Rotor r = sRandomRotor(random);
			Vec4 u = sRandomUnit(random), v = sRandomUnit(random);
			CHECK_APPROX_EQUAL((r * u).Length(), 1.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL((r * u).Dot(r * v), u.Dot(v), 1.0e-4f);
		}
	}

	TEST_CASE("TestRotorComposition")
	{
		UnitTestRandom random(67890);
		for (int n = 0; n < 20; ++n)
		{
			Rotor r1 = sRandomRotor(random), r2 = sRandomRotor(random);
			Vec4 v = sRandomUnit(random);
			CHECK_APPROX_EQUAL((r1 * r2) * v, r1 * (r2 * v), 1.0e-3f);
		}
	}

	TEST_CASE("TestRotorInversed")
	{
		UnitTestRandom random(78901);
		for (int n = 0; n < 20; ++n)
		{
			Rotor r = sRandomRotor(random);
			CHECK((r * r.Inversed()).IsClose(Rotor::sIdentity(), 1.0e-8f));
			Vec4 v = sRandomUnit(random);
			CHECK_APPROX_EQUAL(r.Inversed() * (r * v), v, 1.0e-4f);
		}
	}

	TEST_CASE("TestRotorToRotationMatrix")
	{
		UnitTestRandom random(89012);
		for (int n = 0; n < 20; ++n)
		{
			Rotor r = sRandomRotor(random);
			Mat44 m = r.ToRotationMatrix();

			// Matrix and sandwich product must agree, also for the sRotation(Rotor) constructor
			Vec4 v = sRandomUnit(random);
			CHECK_APPROX_EQUAL(Vec4(m * Lane4(v)), r * v, 1.0e-4f);
			CHECK(Mat44::sRotation(r).IsClose(m, 1.0e-8f));

			// The matrix must be orthonormal
			CHECK((m.Transposed() * m).IsClose(Mat44::sIdentity(), 1.0e-8f));
		}
	}

	// Regression test for CLAUDE.md P0 bug #8 (discovered by this suite): Mat44::GetRotor used to
	// recover the wrong sign for the e23/e24/e34/e1234 (anti-self-dual) components, so it returned
	// the opposite rotation in the e23/e24/e34 planes and sRotation(GetRotor(M)) != M for any
	// rotation not confined to an e1-containing coordinate 2-plane. GetRotor must be a true inverse
	// of sRotation for every rotor, up to the overall +-R double-cover sign.
	TEST_CASE("TestRotorMatrixRoundTrip")
	{
		UnitTestRandom random(90123);
		uniform_real_distribution<float> angle_dist(-JPH_PI, JPH_PI);
		for (int n = 0; n < 20; ++n)
		{
			Vec4 b[4];
			sRandomOrthoBasis(random, b);
			Rotor rotors[3] = {
				Rotor::sRotation(b[0], b[1], angle_dist(random)),											// general-plane simple rotation
				Rotor::sDoubleRotation(b[0], b[1], angle_dist(random), b[2], b[3], angle_dist(random)),	// double rotation
				sRandomRotor(random) };																	// generic rotor
			for (RotorArg r : rotors)
			{
				Rotor back = Mat44::sRotation(r).GetRotor();

				// The recovered rotor must reproduce the same rotation matrix ...
				CHECK(Mat44::sRotation(back).IsClose(Mat44::sRotation(r), 1.0e-5f));

				// ... which is equivalent to being +-r
				bool close = back.IsClose(r, 1.0e-5f) || back.IsClose(-r, 1.0e-5f);
				CHECK(close);
			}
		}
	}

	// IsNormalized must check both rotor constraints, not just the 8-vector length: a unit-length
	// rotor with a nonzero e1234 defect (R ~R = a + b e1234, b != 0) is NOT a rotation.
	TEST_CASE("TestRotorIsNormalizedDetectsDefect")
	{
		// s = p = 1/sqrt(2): 8-vector length is exactly 1, but R ~R = 1 + e1234 (defect 2sp = 1),
		// i.e. |qL|^2 = (s+p)^2 = 2 and |qR|^2 = (s-p)^2 = 0 -- not an element of Spin(4).
		const float k = 1.0f / sqrt(2.0f);
		Rotor defective(k, 0, 0, 0, 0, 0, 0, k);
		CHECK_APPROX_EQUAL(defective.Length(), 1.0f, 1.0e-5f);
		CHECK(abs(Cl4Ref::RotorNormDefect(defective)) > 0.1); // shows the defect is real
		CHECK(!defective.IsNormalized());

		// A genuine unit rotor passes
		CHECK(Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), 0.7f).IsNormalized());
	}

	// Normalized() must project back onto the rotor manifold (both constraints of R ~R = 1), not
	// just rescale the 8-vector. Composition drift makes this the difference between a stable and a
	// slowly exploding simulation.
	TEST_CASE("TestRotorNormalizedRemovesDefect")
	{
		UnitTestRandom random(11223);
		uniform_real_distribution<float> noise(-0.05f, 0.05f);
		for (int n = 0; n < 20; ++n)
		{
			Rotor r = sRandomRotor(random);
			Rotor perturbed = r + Rotor(noise(random), noise(random), noise(random), noise(random), noise(random), noise(random), noise(random), noise(random));
			Rotor normalized = perturbed.Normalized();
			CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormScalar(normalized)), 1.0f, 1.0e-3f);
			CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormDefect(normalized)), 0.0f, 1.0e-3f);
		}
	}

	// LERP must return a valid rotor (unit length AND defect free), not the raw linear blend.
	TEST_CASE("TestRotorLERPReturnsValidRotor")
	{
		UnitTestRandom random(22334);
		for (int n = 0; n < 10; ++n)
		{
			Rotor r1 = sRandomRotor(random), r2 = sRandomRotor(random);
			if (r1.Dot(r2) < 0.0f)
				r2 = -r2;
			Rotor mid = r1.LERP(r2, 0.5f);
			CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormScalar(mid)), 1.0f, 1.0e-3f);
			CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormDefect(mid)), 0.0f, 1.0e-3f);
		}
	}

	TEST_CASE("TestRotorSLERP")
	{
		UnitTestRandom random(33445);
		for (int n = 0; n < 10; ++n)
		{
			Rotor r1 = sRandomRotor(random), r2 = sRandomRotor(random);
			if (r1.Dot(r2) < 0.0f)
				r2 = -r2;

			// End points
			CHECK(r1.SLERP(r2, 0.0f).IsClose(r1, 1.0e-6f));
			CHECK(r1.SLERP(r2, 1.0f).IsClose(r2, 1.0e-6f));
		}

		// Same plane interpolation must give the half angle rotation
		Rotor a = Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), 0.0f);
		Rotor b = Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), 1.0f);
		CHECK(a.SLERP(b, 0.5f).IsClose(Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), 0.5f), 1.0e-8f));

		// SLERP between rotors that differ by a double rotation must stay on the rotor manifold
		// (the defect-aware Normalized keeps the midpoint on Spin(4)). NOTE: this only checks
		// validity, not geodesic optimality; a fully correct 4D slerp would interpolate each SU(2)
		// factor via r1 * sExp(t * (~r1 * r2).Log()).
		Rotor identity = Rotor::sIdentity();
		Rotor iso = Rotor::sDoubleRotation(Vec4::sAxisX(), Vec4::sAxisY(), 1.2f, Vec4::sAxisZ(), Vec4::sAxisW(), 0.9f);
		Rotor mid = identity.SLERP(iso, 0.5f);
		CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormScalar(mid)), 1.0f, 1.0e-3f);
		CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormDefect(mid)), 0.0f, 1.0e-3f);
	}

#ifdef JPH_ROTOR_HAS_EXP_LOG
	// Rotor::sExp / Rotor::Log (general 4D exponential/logarithm, CLAUDE.md P0 bug #1). Unlike the
	// old inlined map in Body::AddRotationStep (which hardcoded the pseudoscalar to 0 and only
	// handled simple bivectors), sExp is correct for double/isoclinic rotations.
	TEST_CASE("TestRotorExpLog")
	{
		UnitTestRandom random(44556);
		// Keep both SU(2) angles |vL|, |vR| < pi (|v| <= sqrt(3)*2*0.6 ~ 2.08) so Log stays on the
		// principal branch; exp-vs-reference below is valid for any bivector regardless of range.
		uniform_real_distribution<float> dist(-0.6f, 0.6f);
		for (int n = 0; n < 50; ++n)
		{
			// Generic bivector: generally a double rotation, exp has a nonzero pseudoscalar part
			Bivec b(dist(random), dist(random), dist(random), dist(random), dist(random), dist(random));
			Rotor r = Rotor::sExp(b);
			sCheckEqualsRef(r, Cl4Ref::Exp(Cl4Ref::FromBivec(b)), 1.0e-4);

			// exp must produce a valid rotor
			CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormScalar(r)), 1.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL(float(Cl4Ref::RotorNormDefect(r)), 0.0f, 1.0e-4f);

			// log is the inverse on the principal branch
			Bivec back = r.Log();
			CHECK(Cl4Ref::FromBivec(back).MaxAbsDiff(Cl4Ref::FromBivec(b)) < 1.0e-3);
		}

		// exp(log(r)) round-trips exactly for any valid rotor, regardless of branch
		for (int n = 0; n < 20; ++n)
		{
			Rotor r = sRandomRotor(random);
			CHECK(Rotor::sExp(r.Log()).IsClose(r, 1.0e-4f));
		}

		// Simple bivector: must reduce to the quaternion-like formula
		Bivec simple(0.7f, 0, 0, 0, 0, 0);
		CHECK(Rotor::sExp(simple).IsClose(Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), 1.4f), 1.0e-8f));
	}
#endif // JPH_ROTOR_HAS_EXP_LOG

	// Extract the 6 bivector components (grade-2 part) of a reference multivector
	static Bivec sBivecFromMV(const Cl4Ref::MV &inMV)
	{
		return Bivec(
			float(inMV.mC[Cl4Ref::cBivecMask[0]]), float(inMV.mC[Cl4Ref::cBivecMask[1]]), float(inMV.mC[Cl4Ref::cBivecMask[2]]),
			float(inMV.mC[Cl4Ref::cBivecMask[3]]), float(inMV.mC[Cl4Ref::cBivecMask[4]]), float(inMV.mC[Cl4Ref::cBivecMask[5]]));
	}

	// Validates Rotor::RotateBivec / InverseRotateBivec (the adjoint action on bivectors) against the
	// brute-force sandwich R B ~R. This is the operation the world-space inverse-inertia operator and
	// the gyroscopic integrator rely on (CLAUDE.md bug #3).
	TEST_CASE("TestRotorRotateBivec")
	{
		UnitTestRandom random(778899);
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		for (int n = 0; n < 50; ++n)
		{
			Rotor r = sRandomRotor(random);
			Bivec b(dist(random), dist(random), dist(random), dist(random), dist(random), dist(random));

			// Brute force: the grade-2 part of R * B * ~R
			Cl4Ref::MV full = Cl4Ref::Sandwich(Cl4Ref::FromRotor(r), Cl4Ref::FromBivec(b));
			Bivec expected = sBivecFromMV(Cl4Ref::GradePart(full, 2));

			Bivec got = r.RotateBivec(b);
			CHECK((got - expected).Length() < 1.0e-4f);

			// The sandwich must be purely grade 2 (no scalar / pseudoscalar leakage)
			CHECK(std::abs(full.mC[0]) < 1.0e-4);
			CHECK(std::abs(full.mC[0b1111]) < 1.0e-4);

			// Inverse round-trips
			Bivec back = r.InverseRotateBivec(got);
			CHECK((back - b).Length() < 1.0e-4f);

			// The adjoint is orthogonal: it preserves the bivector inner product
			Bivec b2(dist(random), dist(random), dist(random), dist(random), dist(random), dist(random));
			CHECK_APPROX_EQUAL(r.RotateBivec(b).Dot(r.RotateBivec(b2)), b.Dot(b2), 1.0e-3f);

			// Consistency with vector rotation: Ad_R(a ^ c) == (R a) ^ (R c)
			Vec4 a = sRandomUnit(random), c = sRandomUnit(random);
			Bivec wedge_then_rotate = r.RotateBivec(Bivec::sWedge(a, c));
			Bivec rotate_then_wedge = Bivec::sWedge(r * a, r * c);
			CHECK((wedge_then_rotate - rotate_then_wedge).Length() < 1.0e-4f);
		}
	}

	// The world-space inverse-inertia operator I^{-1} = Ad(R) D^{-1} Ad(R^{-1}) must be symmetric
	// positive (semi)definite and reduce to the diagonal when the rotation is identity.
	TEST_CASE("TestInverseInertiaAdjoint")
	{
		UnitTestRandom random(445566);
		uniform_real_distribution<float> dist(0.5f, 3.0f);
		uniform_real_distribution<float> vdist(-1.0f, 1.0f);

		// Model the operator the same way MotionProperties does: Ad(R) (D * Ad(R^{-1}) v)
		auto apply = [](RotorArg inR, BivecArg inDiag, BivecArg inV) {
			return inR.RotateBivec(Bivec(inDiag.GetVec8() * inR.InverseRotateBivec(inV).GetVec8()));
		};

		for (int n = 0; n < 30; ++n)
		{
			Rotor r = sRandomRotor(random);
			Bivec diag(dist(random), dist(random), dist(random), dist(random), dist(random), dist(random));

			Bivec u(vdist(random), vdist(random), vdist(random), vdist(random), vdist(random), vdist(random));
			Bivec v(vdist(random), vdist(random), vdist(random), vdist(random), vdist(random), vdist(random));

			// Symmetry: <u, I^{-1} v> == <I^{-1} u, v>
			CHECK_APPROX_EQUAL(u.Dot(apply(r, diag, v)), apply(r, diag, u).Dot(v), 1.0e-3f);

			// Positive definite (diag > 0): <v, I^{-1} v> > 0 for v != 0
			CHECK(v.Dot(apply(r, diag, v)) > 0.0f);
		}

		// Identity rotation reduces to the plain diagonal
		Bivec diag(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f);
		Bivec v(0.3f, -0.7f, 0.2f, 0.9f, -0.1f, 0.5f);
		Bivec expected(diag.GetVec8() * v.GetVec8());
		CHECK((apply(Rotor::sIdentity(), diag, v) - expected).Length() < 1.0e-5f);
	}
}
