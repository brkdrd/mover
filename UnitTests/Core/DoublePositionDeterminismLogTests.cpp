// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#ifndef JPH_ENABLE_DETERMINISM_LOG
#define JPH_ENABLE_DETERMINISM_LOG
#endif
#include <Jolt/Physics/DeterminismLog.h>
#include <cstdio>

TEST_SUITE("DoublePositionDeterminismLogTests")
{
	TEST_CASE("FourthCoordinateIsLoggedBitExactly")
	{
		const char *path = "double-position-determinism-test.log";
		{
			DeterminismLog log(path);
			log << DVec4(1, -2, 3, -4);
		}
		std::ifstream file(path);
		std::string contents;
		std::getline(file, contents);
		CHECK(contents == "3ff0000000000000 c000000000000000 4008000000000000 c010000000000000");
		file.close();
		std::remove(path);
	}
}
