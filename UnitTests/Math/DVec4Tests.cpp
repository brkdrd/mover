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
}
