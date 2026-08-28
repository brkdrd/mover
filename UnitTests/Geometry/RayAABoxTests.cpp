// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/AABox.h>
#include <Jolt/Geometry/RayAABox.h>

TEST_SUITE("RayAABoxTests")
{
	/// Helper to create a Vec4 with 1.0 in a single component
	static Vec4 sAxisVec(int inAxis)
	{
		Vec4 v = Vec4::sZero();
		v.SetX(inAxis == 0 ? 1.0f : 0.0f);
		v.SetY(inAxis == 1 ? 1.0f : 0.0f);
		v.SetZ(inAxis == 2 ? 1.0f : 0.0f);
		v.SetW(inAxis == 3 ? 1.0f : 0.0f);
		return v;
	}

	TEST_CASE("TestRayAABox")
	{
		AABox box(Vec4::sReplicate(-1.0f), Vec4::sOne());

		// Test all 4 axes in 4D
		for (int axis = 0; axis < 4; ++axis)
		{
			Vec4 axis_vec = sAxisVec(axis);

			{
				// Ray starting in the center of the box, pointing high
				float fraction = RayAABox(Vec4::sZero(), RayInvDirection(axis_vec), box.mMin, box.mMax);
				CHECK_APPROX_EQUAL(-1.0f, fraction, 1.0e-6f);
			}

			{
				// Ray starting in the center of the box, pointing low
				float fraction = RayAABox(Vec4::sZero(), RayInvDirection(-axis_vec), box.mMin, box.mMax);
				CHECK_APPROX_EQUAL(-1.0f, fraction, 1.0e-6f);
			}

			{
				// Ray starting high, pointing to low
				Vec4 origin = 1.1f * axis_vec;
				float fraction = RayAABox(origin, RayInvDirection(-axis_vec), box.mMin, box.mMax);
				CHECK_APPROX_EQUAL(0.1f, fraction, 1.0e-6f);
			}

			{
				// Ray starting high, pointing to high
				Vec4 origin = 1.1f * axis_vec;
				float fraction = RayAABox(origin, RayInvDirection(axis_vec), box.mMin, box.mMax);
				CHECK(fraction == FLT_MAX);
			}

			{
				// Ray starting low, pointing to high
				Vec4 origin = -1.1f * axis_vec;
				float fraction = RayAABox(origin, RayInvDirection(axis_vec), box.mMin, box.mMax);
				CHECK_APPROX_EQUAL(0.1f, fraction, 1.0e-6f);
			}

			{
				// Ray starting low, pointing to low
				Vec4 origin = -1.1f * axis_vec;
				float fraction = RayAABox(origin, RayInvDirection(-axis_vec), box.mMin, box.mMax);
				CHECK(fraction == FLT_MAX);
			}
		}

		{
			// Test ray that hits top plane under an angle
			Vec4 expected_hit(0, 1, 0, 0);
			float expected_fraction = 0.123f;
			Vec4 direction(4, -4, 0, 0);
			Vec4 origin = expected_hit - expected_fraction * direction;
			float fraction = RayAABox(origin, RayInvDirection(direction), box.mMin, box.mMax);
			CHECK_APPROX_EQUAL(expected_fraction, fraction, 1.0e-6f);
		}

		{
			// Test diagonal ray in 4D
			Vec4 origin(3, 3, 3, 3);
			Vec4 direction(-1, -1, -1, -1);
			float fraction = RayAABox(origin, RayInvDirection(direction), box.mMin, box.mMax);
			// Should hit at distance 2 (from (3,3,3,3) toward origin, hits at (1,1,1,1))
			CHECK_APPROX_EQUAL(2.0f, fraction, 1.0e-5f);
		}

		{
			// Test ray along W axis hitting the box
			Vec4 origin(0, 0, 0, 5);
			Vec4 direction(0, 0, 0, -1);
			float fraction = RayAABox(origin, RayInvDirection(direction), box.mMin, box.mMax);
			CHECK_APPROX_EQUAL(4.0f, fraction, 1.0e-5f);
		}

		{
			// Test ray missing the box in W dimension
			Vec4 origin(0, 0, 0, 5);
			Vec4 direction(0, 0, 0, 1);
			float fraction = RayAABox(origin, RayInvDirection(direction), box.mMin, box.mMax);
			CHECK(fraction == FLT_MAX);
		}
	}
}
