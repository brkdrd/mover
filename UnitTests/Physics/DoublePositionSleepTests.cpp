// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#include <Jolt/Physics/Body/MotionProperties.h>

#ifdef JPH_DOUBLE_PRECISION
TEST_SUITE("DoublePositionSleepTests")
{
	TEST_CASE("SleepOffsetPreservesAllFourDoubleCoordinates")
	{
		const RVec4 points[3] = {
			RVec4(0x1.0000000000001p40, -0x1.0000000000003p39, 7.25, -0x1.0000000000005p41),
			RVec4(0x1.0000000000001p40 + 1, -0x1.0000000000003p39, 7.25, -0x1.0000000000005p41),
			RVec4(0x1.0000000000001p40, -0x1.0000000000003p39 + 1, 7.25, -0x1.0000000000005p41)
		};
		MotionProperties motion;
		motion.ResetSleepTestSpheres(points);
		const DVec4 offset = motion.GetSleepTestOffset();
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(offset[axis] == points[0][axis]);
	}
}
#endif
