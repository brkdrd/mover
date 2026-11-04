// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "Math/Cl4Reference.h"
#include <Jolt/Math/Bivec.h>
#include <Jolt/Math/Rotor.h>

JPH_SUPPRESS_WARNINGS_STD_BEGIN
#include <random>
#include <cstring>
JPH_SUPPRESS_WARNINGS_STD_END

// Test-first suite for the 4D bivector (CLAUDE.md Phase A step 2). Angular velocity, torque and
// angular momentum are bivectors in 4D, so these operations sit under the whole dynamics stack.
TEST_SUITE("BivecTests")
{
	static Vec4 sRandomVec(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		return Vec4(dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom));
	}

	static Bivec sRandomBivec(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		return Bivec(dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom));
	}

	TEST_CASE("TestBivecWedgeVsReference")
	{
		UnitTestRandom random(13579);
		for (int n = 0; n < 1000; ++n)
		{
			Vec4 a = sRandomVec(random), b = sRandomVec(random);

			// a ^ b is the antisymmetric (grade 2) part of the geometric product a b
			Cl4Ref::MV ref = (Cl4Ref::Mul(Cl4Ref::FromVec4(a), Cl4Ref::FromVec4(b)) - Cl4Ref::Mul(Cl4Ref::FromVec4(b), Cl4Ref::FromVec4(a))) * 0.5;
			CHECK(Cl4Ref::FromBivec(Bivec::sWedge(a, b)).MaxAbsDiff(ref) < 1.0e-5);
		}
	}

	TEST_CASE("TestBivecWedgeProperties")
	{
		UnitTestRandom random(24680);
		for (int n = 0; n < 20; ++n)
		{
			Vec4 a = sRandomVec(random), b = sRandomVec(random);

			// a ^ a = 0, antisymmetry
			CHECK(Bivec::sWedge(a, a).IsNearZero(1.0e-10f));
			CHECK((Bivec::sWedge(a, b) + Bivec::sWedge(b, a)).IsNearZero(1.0e-10f));

			// Bilinearity
			CHECK((Bivec::sWedge(2.0f * Vec4(a), b) - 2.0f * Bivec::sWedge(a, b)).IsNearZero(1.0e-8f));
		}
	}

	// Contract is used as "angular velocity bivector -> velocity of a point": it must equal the
	// derivative of the rotor sandwich rotation, which in the reference algebra is the grade 1
	// part of the commutator (B v - v B) / 2. This ties Bivec::Contract, Rotor::operator* and
	// Body::AddRotationStep to one consistent orientation convention.
	TEST_CASE("TestBivecContractIsRotationGenerator")
	{
		UnitTestRandom random(35791);
		for (int n = 0; n < 1000; ++n)
		{
			Bivec b = sRandomBivec(random);
			Vec4 v = sRandomVec(random);
			Cl4Ref::MV commutator = (Cl4Ref::Mul(Cl4Ref::FromBivec(b), Cl4Ref::FromVec4(v)) - Cl4Ref::Mul(Cl4Ref::FromVec4(v), Cl4Ref::FromBivec(b))) * 0.5;
			CHECK(Cl4Ref::GradePart(commutator, 3).MaxAbs() < 1.0e-6); // commutator of bivector and vector is pure grade 1
			CHECK(Cl4Ref::Grade1Diff(commutator, b.Contract(v)) < 1.0e-5);
		}
	}

	TEST_CASE("TestBivecContractAntisymmetry")
	{
		UnitTestRandom random(46802);
		for (int n = 0; n < 20; ++n)
		{
			Bivec b = sRandomBivec(random);
			Vec4 u = sRandomVec(random), v = sRandomVec(random);

			// The generator is an antisymmetric linear map: v . (B v) = 0 (a point's rotational
			// velocity is perpendicular to its position) and u . (B v) = -v . (B u)
			CHECK_APPROX_EQUAL(v.Dot(b.Contract(v)), 0.0f, 1.0e-5f);
			CHECK_APPROX_EQUAL(u.Dot(b.Contract(v)), -v.Dot(b.Contract(u)), 1.0e-5f);
		}
	}

	TEST_CASE("TestBivecArithmetic")
	{
		Bivec a(1, 2, 3, 4, 5, 6), b(6, 5, 4, 3, 2, 1);
		CHECK(a + b == Bivec(7, 7, 7, 7, 7, 7));
		CHECK(a - b == Bivec(-5, -3, -1, 1, 3, 5));
		CHECK(-a == Bivec(-1, -2, -3, -4, -5, -6));
		CHECK(a * 2.0f == Bivec(2, 4, 6, 8, 10, 12));
		CHECK(a * b == Bivec(6, 10, 12, 12, 10, 6)); // component-wise (diagonal inertia application)
		CHECK_APPROX_EQUAL(a.Dot(b), 56.0f);
		CHECK_APPROX_EQUAL(a.LengthSq(), 91.0f);
		CHECK_APPROX_EQUAL(Bivec(3, 0, 0, 4, 0, 0).Length(), 5.0f);
		CHECK_APPROX_EQUAL(a.Normalized().Length(), 1.0f, 1.0e-5f);
		CHECK(Bivec::sZero().IsNearZero());
	}

	// Regression test for CLAUDE.md P1 bug #10: Bivec::sAnd used to type-pun a float[8] through a
	// uint32* and mutate it in place, which is strict-aliasing UB and was miscompiled at -O2/-O3
	// (sAnd returned its first argument unmodified). The memcpy-based implementation masks per
	// component correctly.
	TEST_CASE("TestBivecSAnd")
	{
		float all_bits;
		uint32 ones = 0xffffffffu;
		memcpy(&all_bits, &ones, sizeof(float));

		Bivec v(1, -2, 3, -4, 5, -6);
		Bivec mask(all_bits, 0, all_bits, 0, 0, all_bits);
		CHECK(Bivec::sAnd(v, mask) == Bivec(1, 0, 3, 0, 0, -6));
		CHECK(Bivec::sAnd(v, Bivec::sZero()) == Bivec::sZero());
	}

	// Reference-level validation of the invariant (self-dual / anti-self-dual) decomposition that
	// the general rotor exponential fix will use (CLAUDE.md P0 bug #1): any bivector B splits via
	// the pseudoscalar I = e1234 into B+- = (B +- B I) / 2 with B+ B- = 0 and [B+, B-] = 0, so
	// exp(B) = exp(B+) exp(B-) where each factor exponentiates like a quaternion.
	TEST_CASE("TestBivecSelfDualSplitReference")
	{
		using namespace Cl4Ref;
		MV pseudo;
		pseudo.mC[0b1111] = 1.0;

		UnitTestRandom random(57913);
		for (int n = 0; n < 20; ++n)
		{
			MV b = FromBivec(sRandomBivec(random));
			MV b_dual = Mul(b, pseudo);
			CHECK(GradePart(b_dual, 2).MaxAbsDiff(b_dual) < 1.0e-12); // duality maps bivectors to bivectors
			MV bp = (b + b_dual) * 0.5;
			MV bm = (b - b_dual) * 0.5;

			// Split reconstructs B, parts annihilate and commute
			CHECK((bp + bm).MaxAbsDiff(b) < 1.0e-12);
			CHECK(Mul(bp, bm).MaxAbs() < 1.0e-12);
			CHECK(Mul(bm, bp).MaxAbs() < 1.0e-12);

			// exp(B) = exp(B+) exp(B-)
			CHECK(Exp(b).MaxAbsDiff(Mul(Exp(bp), Exp(bm))) < 1.0e-12);

			// exp(B) is a valid rotor: R ~R = 1 exactly (scalar 1, no defect)
			MV r = Exp(b);
			MV rr = Mul(r, Reverse(r));
			CHECK(std::abs(rr.mC[0] - 1.0) < 1.0e-12);
			CHECK(std::abs(rr.mC[0b1111]) < 1.0e-12);
		}

		// Demonstrates why the simple-bivector formula in Body::AddRotationStep is insufficient
		// (premise of CLAUDE.md P0 bug #1): for a double rotation the true exponential has a
		// nonzero pseudoscalar part, which "cos|B| + sin|B| B / |B|" can never produce.
		MV b;
		b.mC[0b0011] = 1.0; // e12
		b.mC[0b1100] = 0.7; // e34
		MV true_exp = Exp(b);
		CHECK(std::abs(true_exp.mC[0b1111]) > 0.5); // = sin(1) * sin(0.7) ~ 0.542

		double len = std::sqrt(1.0 * 1.0 + 0.7 * 0.7);
		MV simple_formula = b * (std::sin(len) / len);
		simple_formula.mC[0] = std::cos(len);
		CHECK(true_exp.MaxAbsDiff(simple_formula) > 0.1);
	}
}
