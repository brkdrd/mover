// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Core/StringTools.h>
#include <Jolt/Math/RMat44.h>
#include <Jolt/Physics/Collision/RayCast.h>

namespace
{
// Do not use DVec4 equality/length to check DVec4: those once ignored W too.
void CheckComponents(DVec4Arg inV, double inX, double inY, double inZ, double inW)
{
	CHECK(inV.GetX() == inX);
	CHECK(inV.GetY() == inY);
	CHECK(inV.GetZ() == inZ);
	CHECK(inV.GetW() == inW);
}
}

TEST_SUITE("DVec4Tests")
{
	TEST_CASE("IndependentCoordinates")
	{
		CheckComponents(DVec4::sZero(), 0, 0, 0, 0);
		CheckComponents(DVec4::sOne(), 1, 1, 1, 1);
		CheckComponents(DVec4::sReplicate(2), 2, 2, 2, 2);
		CheckComponents(DVec4::sAxisX(), 1, 0, 0, 0);
		CheckComponents(DVec4::sAxisY(), 0, 1, 0, 0);
		CheckComponents(DVec4::sAxisZ(), 0, 0, 1, 0);
		CheckComponents(DVec4::sAxisW(), 0, 0, 0, 1);
		CheckComponents(DVec4(1, 2, 3), 1, 2, 3, 0); // legacy XYZ embedding
		CheckComponents(DVec4(Vec4(1, 2, 3, 7)), 1, 2, 3, 7);
		CheckComponents(DVec4(Lane4(1, 2, 3, 7)), 1, 2, 3, 7);
		DVec4 v(1, 2, 3, 7);
		v.SetZ(9);
		CheckComponents(v, 1, 2, 9, 7);
		v.SetComponent(2, 11);
		CheckComponents(v, 1, 2, 11, 7);
		v.SetComponent(3, 13);
		CHECK(v[3] == 13);
		v.SetX(4);
		v.SetY(5);
		v.SetW(6);
		CheckComponents(v, 4, 5, 11, 6);
		Double4 stored;
		v.StoreDouble4(&stored);
		CHECK(stored == Double4(4, 5, 11, 6));
		CheckComponents(DVec4(stored), 4, 5, 11, 6);
		CHECK(Double4(4, 5, 11, 6) != Double4(4, 5, 11, 7));
		CHECK(stored[3] == 6);
		Vec4 f(v);
		CHECK(f.GetW() == 6);
		CHECK(ConvertToString(v) == "4, 5, 11, 6");
	}

	TEST_CASE("FourLaneArithmeticAndReductions")
	{
		const DVec4 a(2, -3, 4, -8), b(5, 6, 7, 16);
		CheckComponents(a + b, 7, 3, 11, 8);
		CheckComponents(a - b, -3, -9, -3, -24);
		CheckComponents(a * b, 10, -18, 28, -128);
		CheckComponents(a / b, 2.0 / 5, -0.5, 4.0 / 7, -0.5);
		CheckComponents(-a, -2, 3, -4, 8);
		CheckComponents(2 * a, 4, -6, 8, -16);
		CheckComponents(a / 2, 1, -1.5, 2, -4);
		DVec4 v = a;
		v *= b;
		v *= 2;
		v /= 2;
		v += b;
		v -= b;
		CheckComponents(v, 10, -18, 28, -128);
		CheckComponents(DVec4::sFusedMultiplyAdd(a, b, DVec4::sOne()), 11, -17, 29, -127);
		CheckComponents(a.Abs(), 2, 3, 4, 8);
		CheckComponents(a.GetSign(), 1, -1, 1, -1);
		CheckComponents(DVec4(1, 4, 9, 16).Sqrt(), 1, 2, 3, 4);
		CheckComponents(DVec4(1, 2, 4, 8).Reciprocal(), 1, 0.5, 0.25, 0.125);
		CHECK(a.Dot(b) == 10 - 18 + 28 - 128);
		CHECK(a.LengthSq() == 4 + 9 + 16 + 64);
		CHECK(a.Length() == sqrt(93.0));
		CheckComponents(DVec4(0, 0, 0, 2).Normalized(), 0, 0, 0, 1);
		CHECK(DVec4::sAxisW().IsNormalized());
		CHECK_FALSE(DVec4(0, 0, 0, 2).IsNormalized());
		CHECK_FALSE(DVec4::sAxisW().IsNearZero());
		CHECK_FALSE(DVec4::sZero().IsClose(DVec4::sAxisW()));
		CheckComponents(DVec4::sMin(a, b), 2, -3, 4, -8);
		CheckComponents(DVec4::sMax(a, b), 5, 6, 7, 16);
		CheckComponents(DVec4::sClamp(a, DVec4::sZero(), DVec4::sReplicate(3)), 2, 0, 3, 0);
	}

	TEST_CASE("FourLaneComparisonsAndMasks")
	{
		for (uint mask = 0; mask < 16; ++mask)
		{
			DVec4 control(BitCast<double>(uint64(mask & 1) << 63), BitCast<double>(uint64((mask >> 1) & 1) << 63),
				BitCast<double>(uint64((mask >> 2) & 1) << 63), BitCast<double>(uint64((mask >> 3) & 1) << 63));
			CHECK(control.GetTrues() == int(mask));
			CHECK(control.TestAnyTrue() == (mask != 0));
			CHECK(control.TestAllTrue() == (mask == 15));
			DVec4 selected = DVec4::sSelect(DVec4(1, 2, 3, 4), DVec4(5, 6, 7, 8), control);
			for (uint axis = 0; axis < 4; ++axis)
				CHECK(selected[axis] == double(axis + ((mask & (1U << axis))? 5 : 1)));
		}
		DVec4 a(1, 2, 3, 4), b(1, 4, 2, 5);
		CHECK_FALSE(a == DVec4(1, 2, 3, 5));
		CHECK(a != DVec4(1, 2, 3, 5));
		CHECK(DVec4::sEquals(a, b).GetTrues() == 1);
		CHECK(DVec4::sLess(a, b).GetTrues() == 10);
		CHECK(DVec4::sLessOrEqual(a, b).GetTrues() == 11);
		CHECK(DVec4::sGreater(a, b).GetTrues() == 4);
		CHECK(DVec4::sGreaterOrEqual(a, b).GetTrues() == 5);
		DVec4 nan = DVec4::sZero();
		nan.SetW(numeric_limits<double>::quiet_NaN());
		CHECK(nan.IsNaN());
		CheckComponents(DVec4::sAnd(DVec4::sOne(), DVec4::sEquals(a, a)), 1, 1, 1, 1);
		CheckComponents(DVec4::sXor(a, a), 0, 0, 0, 0);
		CheckComponents(DVec4::sOr(DVec4::sZero(), a), 1, 2, 3, 4);
	}

#ifdef JPH_DOUBLE_PRECISION
	TEST_CASE("DoubleTransformsKeepSubFloatResolutionInEveryAxis")
	{
		// An exact SO(4) permutation, mixing W into X and X into W.
		const Mat44 rotation(Lane4(0, 0, 0, -1), Lane4(0, 1, 0, 0), Lane4(0, 0, 1, 0), Lane4(1, 0, 0, 0));
		const DVec4 p(0x1.0000000000001p40, -0x1.0000000000003p39, 0x1.0000000000005p38, -0x1.0000000000007p41);
		const DVec4 rotated(p.GetW(), p.GetY(), p.GetZ(), -p.GetX());
		CheckComponents(RMat44::sRotationTranslation(Rotor::sIdentity(), Vec4(1, 2, 3, 11)).GetTranslation(), 1, 2, 3, 11);
		CheckComponents(RMat44::sInverseRotationTranslation(Rotor::sIdentity(), Vec4(1, 2, 3, 11)).GetTranslation(), -1, -2, -3, -11);
		RMat44 transform(rotation, DVec4(1, 2, 3, 4));
		DVec4 result = transform * p;
		CheckComponents(result, rotated.GetX() + 1, rotated.GetY() + 2, rotated.GetZ() + 3, rotated.GetW() + 4);
		DVec4 round_trip = transform.InversedRotationTranslation() * result;
		CheckComponents(round_trip, p.GetX(), p.GetY(), p.GetZ(), p.GetW());
		RMat44 translated(rotation, p);
		CheckComponents(translated.InversedRotationTranslation().GetTranslation(), p.GetW(), -p.GetY(), -p.GetZ(), -p.GetX());
		RMat44 composed = RMat44(rotation) * RMat44::sTranslation(p);
		CheckComponents(composed.GetTranslation(), rotated.GetX(), rotated.GetY(), rotated.GetZ(), rotated.GetW());
		CheckComponents(translated.PostScaled(Vec4(2, 4, 8, 16)).GetTranslation(), 2 * p.GetX(), 4 * p.GetY(), 8 * p.GetZ(), 16 * p.GetW());
		CheckComponents(translated.Inversed().GetTranslation(), p.GetW(), -p.GetY(), -p.GetZ(), -p.GetX());
	}
#endif

	TEST_CASE("RayDirectionDoesNotSubtractLargeTranslatedEndpoints")
	{
		RMat44 transform = RMat44::sTranslation(RVec4::sReplicate(0x1p40));
		RayCast ray(Vec4(1, 2, 3, 4), Vec4(0.25f, -0.5f, 0.75f, -1.25f));
		RayCast transformed = ray.Transformed(transform);
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(transformed.mDirection[axis] == ray.mDirection[axis]);
	}

	TEST_CASE("DirectedRoundingPreservesFourthCoordinate")
	{
		const DVec4 value(2.0, 0x1.0000000000001p1, -0x1.0000000000001p1, 0x1.0000000000001p2);
		Vec4 up = value.ToVec4RoundUp(), down = value.ToVec4RoundDown();
		const Vec4 expected_up(2, 0x1.000002p1f, -2, 0x1.000002p2f);
		const Vec4 expected_down(2, 2, -0x1.000002p1f, 4);
		for (uint axis = 0; axis < 4; ++axis)
		{
			CHECK(up[axis] == expected_up[axis]);
			CHECK(down[axis] == expected_down[axis]);
			// The old method names are compatibility aliases, not XYZ projections.
			CHECK(value.ToVec3RoundUp().mF32[axis] == expected_up[axis]);
			CHECK(value.ToVec3RoundDown().mF32[axis] == expected_down[axis]);
		}
	}

	TEST_CASE("FloatDisplacementsAndHashPreserveW")
	{
		const DVec4 value(1, 2, 3, 17);
		const Vec4 delta(4, 5, 6, 19);
		CheckComponents(value + delta, 5, 7, 9, 36);
		CheckComponents(value - delta, -3, -3, -3, -2);
		DVec4 updated = value;
		updated += delta;
		CheckComponents(updated, 5, 7, 9, 36);
		updated -= delta;
		CheckComponents(updated, 1, 2, 3, 17);
		CHECK(Hash<DVec4>()(value) == Hash<Double4>()(Double4(1, 2, 3, 17)));
		CHECK(Hash<DVec4>()(value) != Hash<DVec4>()(DVec4(1, 2, 3, 19)));
	}

	TEST_CASE("CrossProductUsesAllFourAxes")
	{
		CheckComponents(DVec4::sCross(DVec4::sAxisX(), DVec4::sAxisY(), DVec4::sAxisZ()), 0, 0, 0, -1);
		const DVec4 a(1, 2, 3, 5), b(-2, 1, 4, -3), c(7, -1, 2, 6);
		DVec4 cross = DVec4::sCross(a, b, c);
		CHECK(cross.Dot(a) == 0);
		CHECK(cross.Dot(b) == 0);
		CHECK(cross.Dot(c) == 0);
	}

	// -----------------------------------------------------------------------
	// Compatibility coverage restored from the pre-four-lane tests. Expected
	// values are stated in the four-lane contract (XYZ-only inputs embed W = 0)
	// instead of being relaxed, so a lost W or a copied fourth lane still fails.

	TEST_CASE("LegacyXYZEmbeddingRemainsZeroInW")
	{
		// Three-lane constructors, setters and operations all embed W = 0.
		DVec4 v(1, 2, 3);
		CheckComponents(v, 1, 2, 3, 0);
		v.SetComponent(0, 4);
		v.SetComponent(1, 5);
		v.SetComponent(2, 6);
		CheckComponents(v, 4, 5, 6, 0);
		v.SetX(7);
		v.SetY(8);
		v.SetZ(9);
		CheckComponents(v, 7, 8, 9, 0);
		v.SetComponent(3, 13);
		v.Set(10, 11, 12); // Legacy XYZ setter resets W to 0
		CheckComponents(v, 10, 11, 12, 0);

		CheckComponents(DVec4(1, 2, 3) + DVec4(4, 5, 6), 5, 7, 9, 0);
		CheckComponents(DVec4(1, 2, 3) - DVec4(6, 5, 4), -5, -3, -1, 0);
		CheckComponents(DVec4(1, 2, 3) * DVec4(4, 5, 6), 4, 10, 18, 0);
		// Division is four-lane like every other operator: embedding both
		// operands with W = 0 makes the quotient's W a 0/0 NaN, which is the
		// documented four-lane result rather than a silently projected zero.
		DVec4 divided = DVec4(1, 2, 3) / DVec4(2, 8, 24);
		CHECK(divided.GetX() == 0.5);
		CHECK(divided.GetY() == 0.25);
		CHECK(divided.GetZ() == 0.125);
		CHECK(isnan(divided.GetW()));
		CheckComponents(DVec4(1, 2, 3, 8) / DVec4(2, 8, 24, 4), 0.5, 0.25, 0.125, 2);

		CheckComponents(DVec4::sMin(DVec4(1, 5, 3), DVec4(4, 2, 6)), 1, 2, 3, 0);
		CheckComponents(DVec4::sMax(DVec4(1, 5, 3), DVec4(4, 2, 6)), 4, 5, 6, 0);
		CheckComponents(DVec4::sClamp(DVec4(-1, 3, 7), DVec4(1, 2, 3), DVec4(4, 5, 6)), 1, 3, 6, 0);
	}

	TEST_CASE("Vec3ProjectionAndVec3OperandArithmetic")
	{
		// Explicit Vec3 conversion stays a deliberate XYZ projection.
		CHECK(Vec3(DVec4(1, 3, 5, 7)) == Vec3(1, 3, 5));

		// A Vec3 operand is embedded with W = 0, so it does not disturb W.
		CheckComponents(DVec4(1, 2, 3, 4) + Vec3(4, 5, 6), 5, 7, 9, 4);
		CheckComponents(DVec4(1, 2, 3, 4) - Vec3(6, 5, 4), -5, -3, -1, 4);
		DVec4 v(1, 2, 3, 4);
		v += Vec3(1, 2, 3);
		CheckComponents(v, 2, 4, 6, 4);
		v -= Vec3(1, 2, 3);
		CheckComponents(v, 1, 2, 3, 4);
	}

	TEST_CASE("NaNInEveryLane")
	{
		DVec4 v = DVec4::sNaN();
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(isnan(v[axis]));
		CHECK(v.IsNaN());

		// Clearing only XYZ must not report a NaN-free vector: W counts too.
		for (uint axis = 0; axis < 3; ++axis)
			v.SetComponent(axis, 0);
		CHECK(v.IsNaN());
		v.SetComponent(3, 0);
		CHECK_FALSE(v.IsNaN());
	}

	TEST_CASE("Double3LoadAndStoreAreXYZOnly")
	{
		struct
		{
			Double3 mValue;
			double mSentinel;
		} source = { Double3(1, 2, 3), 4.0 };
		DVec4 loaded = DVec4::sLoadDouble3Unsafe(source.mValue);
		CheckComponents(loaded, 1, 2, 3, 0); // The readable fourth double is not W

		struct
		{
			Double3 mValue;
			double mSentinel;
		} destination = { Double3(9, 9, 9), 9.0 };
		DVec4(1, 2, 3, 42).StoreDouble3(&destination.mValue);
		CHECK(destination.mValue == Double3(1, 2, 3));
		CHECK(destination.mSentinel == 9.0); // StoreDouble3 must not write a fourth double
	}

	TEST_CASE("BinaryCrossIsXYZOnly")
	{
		CHECK(DVec4(1, 0, 0).Cross(DVec4(0, 1, 0)) == DVec4(0, 0, 1));
		CHECK(DVec4(0, 1, 0).Cross(DVec4(1, 0, 0)) == DVec4(0, 0, -1));
		CHECK(DVec4(0, 1, 0).Cross(DVec4(0, 0, 1)) == DVec4(1, 0, 0));
		CHECK(DVec4(0, 0, 1).Cross(DVec4(0, 1, 0)) == DVec4(-1, 0, 0));
		CHECK(DVec4(0, 0, 1).Cross(DVec4(1, 0, 0)) == DVec4(0, 1, 0));
		CHECK(DVec4(1, 0, 0).Cross(DVec4(0, 0, 1)) == DVec4(0, -1, 0));

		// The 3D cross product ignores W and returns an XYZ result (W = 0).
		CheckComponents(DVec4(1, 0, 0, 7).Cross(DVec4(0, 1, 0, 9)), 0, 0, 1, 0);
		CheckComponents(DVec4(0, 1, 0, 7).Cross(DVec4(1, 0, 0, 9)), 0, 0, -1, 0);
	}

	TEST_CASE("ToleranceBoundariesIncludeW")
	{
		CHECK(DVec4(1, 2, 3).IsClose(DVec4(1.001, 2.001, 3.001), 1.0e-4));
		CHECK_FALSE(DVec4(1, 2, 3).IsClose(DVec4(1.001, 2.001, 3.001), 1.0e-6));

		// W participates in the distance; a W-only difference is not "close".
		CHECK(DVec4(1, 2, 3, 0).IsClose(DVec4(1, 2, 3, 0.02), 1.0e-4) == false);
		CHECK(DVec4(0, 0, 0, 0.5).IsClose(DVec4::sZero(), 0.25)); // Threshold is inclusive
		CHECK_FALSE(DVec4(0, 0, 0, 0.5).IsClose(DVec4::sZero(), 0x1.fffffffffffffp-3));

		CHECK(DVec4(1.001, 0, 0, 0).IsNormalized(1.0e-2));
		CHECK_FALSE(DVec4(0, 1.001, 0, 0).IsNormalized(1.0e-4));
		CHECK_FALSE(DVec4(0, 0, 0, 1.001).IsNormalized(1.0e-4)); // W counts towards length

		CHECK(DVec4(-1.0e-7, 1.0e-7, 1.0e-8).IsNearZero(1.0e-12));
		CHECK_FALSE(DVec4(-1.0e-7, 1.0e-7, -1.0e-5).IsNearZero(1.0e-12));
		CHECK_FALSE(DVec4(0, 0, 0, -1.0e-5).IsNearZero(1.0e-12));
		CHECK(DVec4(0, 0, 0, 0.5).IsNearZero(0.25)); // Threshold is inclusive
		CHECK_FALSE(DVec4(0, 0, 0, 0.5).IsNearZero(0x1.fffffffffffffp-3));
	}

	TEST_CASE("BitOperationsOnNontrivialFourLanePatterns")
	{
		// Distinct patterns in every lane: a lane copied from its neighbour or
		// a dropped fourth lane changes the expected result.
		const uint64 patterns[] = { 0b0011, 0b00110, 0b001100, 0b0011000 };
		const uint64 patterns2[] = { 0b0101, 0b01010, 0b010100, 0b0101000 };
		DVec4 v1(BitCast<double>(patterns[0]), BitCast<double>(patterns[1]), BitCast<double>(patterns[2]), BitCast<double>(patterns[3]));
		DVec4 v2(BitCast<double>(patterns2[0]), BitCast<double>(patterns2[1]), BitCast<double>(patterns2[2]), BitCast<double>(patterns2[3]));

		for (uint axis = 0; axis < 4; ++axis)
		{
			CHECK(DVec4::sOr(v1, v2)[axis] == BitCast<double>(patterns[axis] | patterns2[axis]));
			CHECK(DVec4::sXor(v1, v2)[axis] == BitCast<double>(patterns[axis] ^ patterns2[axis]));
			CHECK(DVec4::sAnd(v1, v2)[axis] == BitCast<double>(patterns[axis] & patterns2[axis]));
		}
	}

	TEST_CASE("NegatedZeroAndSignPolicy")
	{
		DVec4 neg_zero = -DVec4::sZero();
		CHECK(neg_zero == DVec4::sZero());

#if defined(JPH_CROSS_PLATFORM_DETERMINISTIC)
		// Deterministic builds normalize a negated zero to +0.
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(BitCast<uint64>(neg_zero[axis]) == 0);
#elif defined(JPH_USE_AVX) || defined(JPH_USE_SSE)
		// Preexisting x86 policy: negation is implemented as 0 - x, which maps
		// -0 to +0. This is the observed backend behavior, not a portable
		// bit-exact promise.
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(BitCast<uint64>(neg_zero[axis]) == 0);
#endif

		const DVec4 zeros(0.0, -0.0, 0.0, -0.0);
#ifdef JPH_USE_AVX512
		// Preexisting AVX512 policy: vfixupimm classifies -0 as positive.
		CheckComponents(zeros.GetSign(), 1, 1, 1, 1);
#else
		CheckComponents(zeros.GetSign(), 1, -1, 1, -1);
#endif

		DVec4 abs_zeros = zeros.Abs();
#if defined(JPH_USE_SSE) && !defined(JPH_USE_AVX512)
		// Preexisting SSE/AVX policy: Abs is max(-x, x), and for equal zeros the
		// max returns the second operand, so the sign of Abs(-0) stays negative.
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(BitCast<uint64>(abs_zeros[axis]) == (axis & 1? 0x8000000000000000ULL : 0));
#else
		// AVX512 vrangepd, NEON vabsq, RVV and scalar abs return +0.
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(BitCast<uint64>(abs_zeros[axis]) == 0);
#endif
	}

	TEST_CASE("IntermediateResultsBeforeCancellation")
	{
		// Check the products before subtracting them from each other.
		const DVec4 a(2, -3, 4, -8), b(5, 6, 7, 16);
		const DVec4 ab = a * b;
		CheckComponents(ab, 10, -18, 28, -128);
		const DVec4 ba = b * a;
		CheckComponents(ba, 10, -18, 28, -128);
		CheckComponents(ab - ba, 0, 0, 0, 0);

		// Check every in-place update, not only the state after it cancels.
		DVec4 v(1, 2, 3, 4);
		v *= b;
		CheckComponents(v, 5, 12, 21, 64);
		v *= 2;
		CheckComponents(v, 10, 24, 42, 128);
		v /= 2;
		CheckComponents(v, 5, 12, 21, 64);
		v += b;
		CheckComponents(v, 10, 18, 28, 80);
		v -= b;
		CheckComponents(v, 5, 12, 21, 64);

		// Large translated operands keep their sub-float increments before the
		// translation is subtracted again.
		const DVec4 large(0x1p40, -0x1p41, 0x1p42, -0x1p43);
		const DVec4 sum = large + b;
		CheckComponents(sum, 0x1p40 + 5, -0x1p41 + 6, 0x1p42 + 7, -0x1p43 + 16);
		CheckComponents(sum - b, 0x1p40, -0x1p41, 0x1p42, -0x1p43);

		// All four products are nonzero even though the dot product cancels.
		const DVec4 p(1, 2, 3, 4), q(-4, 3, -2, 1);
		const DVec4 pq = p * q;
		CheckComponents(pq, -4, 6, -6, 4);
		CHECK(p.Dot(q) == 0);
		CHECK(pq.GetX() + pq.GetY() + pq.GetZ() + pq.GetW() == 0);
	}
}
