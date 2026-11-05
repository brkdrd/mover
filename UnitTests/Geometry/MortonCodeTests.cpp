// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#include <Jolt/Geometry/MortonCode.h>
#include <array>

namespace
{
// Independent bit-by-bit reference: W is least significant, X most significant.
uint32 ReferenceMorton(const uint32 (&inCoordinates)[4])
{
	uint32 code = 0;
	for (uint bit = 0; bit < 8; ++bit)
		for (uint axis = 0; axis < 4; ++axis)
			code |= ((inCoordinates[axis] >> bit) & 1U) << (4 * bit + 3 - axis);
	return code;
}
}

TEST_SUITE("MortonCodeTests")
{
	TEST_CASE("EveryQuantizedBitHasItsOwnLane")
	{
		for (uint32 value = 0; value < 256; ++value)
		{
			uint32 expected = 0;
			for (uint bit = 0; bit < 8; ++bit)
				expected |= ((value >> bit) & 1U) << (4 * bit);
			CHECK(MortonCode::sExpandBits(float(value) / 255.0f) == expected);
		}
		const AABox bounds(Vec4::sZero(), Vec4::sReplicate(255));
		CHECK(MortonCode::sGetMortonCode(Vec4(1, 0, 0, 0), bounds) != MortonCode::sGetMortonCode(Vec4(0, 0, 0, 4), bounds));
		CHECK(MortonCode::sGetMortonCode(Vec4(2, 0, 0, 0), bounds) != MortonCode::sGetMortonCode(Vec4(0, 0, 0, 8), bounds));
		CHECK(MortonCode::sGetMortonCode(Vec4::sReplicate(255), bounds) == 0xffffffffU);
	}

	TEST_CASE("ExhaustiveFourAxisGridIsInjective")
	{
		const AABox bounds(Vec4::sZero(), Vec4::sReplicate(255));
		std::array<uint32, 65536> codes;
		for (uint32 index = 0; index < codes.size(); ++index)
		{
			uint32 q[4] = { (index >> 12) & 15U, (index >> 8) & 15U, (index >> 4) & 15U, index & 15U };
			codes[index] = MortonCode::sGetMortonCode(Vec4(float(q[0]), float(q[1]), float(q[2]), float(q[3])), bounds);
			CHECK(codes[index] == ReferenceMorton(q));
		}
		std::sort(codes.begin(), codes.end());
		CHECK(std::adjacent_find(codes.begin(), codes.end()) == codes.end());
	}

	TEST_CASE("FullRangeMatchesIndependentReference")
	{
		const AABox bounds(Vec4(-20, 100, 7, -1000), Vec4(235, 355, 262, -745));
		UnitTestRandom random(0x4d);
		for (uint sample = 0; sample < 10000; ++sample)
		{
			uint32 q[4];
			for (uint axis = 0; axis < 4; ++axis)
				q[axis] = uint32(random() & 255U);
			Vec4 point = bounds.mMin + Vec4(float(q[0]), float(q[1]), float(q[2]), float(q[3]));
			CHECK(MortonCode::sGetMortonCode(point, bounds) == ReferenceMorton(q));
		}
	}

	TEST_CASE("FlatBoundsDoNotDivideByZero")
	{
		const AABox flat(Vec4(1, 2, 3, 4), Vec4(1, 257, 3, 259));
		CHECK(MortonCode::sGetMortonCode(flat.mMin, flat) == 0);
		const uint32 q[4] = { 0, 255, 0, 255 };
		CHECK(MortonCode::sGetMortonCode(flat.mMax, flat) == ReferenceMorton(q));
		const AABox point(flat.mMin, flat.mMin);
		CHECK(MortonCode::sGetMortonCode(point.mMin, point) == 0);
	}
}
