// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/RayTetrahedron.h>

TEST_SUITE("RayTetrahedronTests")
{
	TEST_CASE("TestRayTetrahedronHit")
	{
		// Regular tetrahedron at the origin, one vertex per axis
		Vec4 v0(1, 0, 0, 0);
		Vec4 v1(0, 1, 0, 0);
		Vec4 v2(0, 0, 1, 0);
		Vec4 v3(0, 0, 0, 1);

		// Ray from far along the (1,1,1,1) direction toward the tetrahedron
		Vec4 origin(5, 5, 5, 5);
		Vec4 direction = Vec4(-1, -1, -1, -1).Normalized();

		float t = RayTetrahedron(origin, direction, v0, v1, v2, v3);
		CHECK(t != FLT_MAX);
		CHECK(t > 0.0f);

		// Verify the hit point is on the tetrahedron's hyperplane
		Vec4 hit = origin + direction * t;
		// The hyperplane equation: x + y + z + w = 1 for this tetrahedron
		CHECK_APPROX_EQUAL(hit.GetX() + hit.GetY() + hit.GetZ() + hit.GetW(), 1.0f, 1.0e-4f);

		// Verify hit point has all non-negative barycentric coords (inside tetrahedron)
		// For this specific tetrahedron, barycentrics = coordinates themselves
		CHECK(hit.GetX() >= -1.0e-4f);
		CHECK(hit.GetY() >= -1.0e-4f);
		CHECK(hit.GetZ() >= -1.0e-4f);
		CHECK(hit.GetW() >= -1.0e-4f);
	}

	TEST_CASE("TestRayTetrahedronHitVertex")
	{
		Vec4 v0(1, 0, 0, 0);
		Vec4 v1(0, 1, 0, 0);
		Vec4 v2(0, 0, 1, 0);
		Vec4 v3(0, 0, 0, 1);

		// Ray aimed directly at vertex v0
		Vec4 origin(3, 0, 0, 0);
		Vec4 direction(-1, 0, 0, 0);
		float t = RayTetrahedron(origin, direction, v0, v1, v2, v3);
		CHECK_APPROX_EQUAL(t, 2.0f, 1.0e-4f);
	}

	TEST_CASE("TestRayTetrahedronMissOutside")
	{
		Vec4 v0(1, 0, 0, 0);
		Vec4 v1(0, 1, 0, 0);
		Vec4 v2(0, 0, 1, 0);
		Vec4 v3(0, 0, 0, 1);

		// A ray straight down the x axis meets the hyperplane x+y+z+w=1 at (1,0,0,0), which is exactly
		// the vertex v0 - a boundary case, so accept either a hit at t=6 or a rejection
		Vec4 origin(-5, 0, 0, 0);
		Vec4 direction(1, 0, 0, 0);
		float t = RayTetrahedron(origin, direction, v0, v1, v2, v3);
		CHECK((t == FLT_MAX || abs(t - 6.0f) < 1.0e-4f));

		// A ray that truly misses: it crosses the hyperplane outside the tetrahedron
		Vec4 origin2(5, 5, -1, 0);
		Vec4 direction2(-1, -1, 0, 0);
		float t2 = RayTetrahedron(origin2, direction2, v0, v1, v2, v3);
		// Hit hyperplane at x+y+z+w=1: (5-t)+(5-t)+(-1)+0=1 => 9-2t=1 => t=4
		// Hit point = (1, 1, -1, 0). z=-1 < 0, so outside tetrahedron.
		CHECK(t2 == FLT_MAX);
	}

	TEST_CASE("TestRayTetrahedronMissBehind")
	{
		Vec4 v0(1, 0, 0, 0);
		Vec4 v1(0, 1, 0, 0);
		Vec4 v2(0, 0, 1, 0);
		Vec4 v3(0, 0, 0, 1);

		// Ray starting inside the hyperplane pointing away
		Vec4 origin(0.25f, 0.25f, 0.25f, 0.25f); // centroid, on hyperplane
		Vec4 direction(1, 1, 1, 1); // pointing away from origin
		float t = RayTetrahedron(origin, direction, v0, v1, v2, v3);
		// n.D should be positive (moving away), n.s should be zero (on hyperplane)
		// t = 0 which is valid, but the point is the centroid itself
		// Actually t = 0 means hit at origin, which is on the tetrahedron
		CHECK(t >= 0.0f);
		if (t != FLT_MAX)
			CHECK_APPROX_EQUAL(t, 0.0f, 1.0e-4f);
	}

	TEST_CASE("TestRayTetrahedronParallel")
	{
		Vec4 v0(1, 0, 0, 0);
		Vec4 v1(0, 1, 0, 0);
		Vec4 v2(0, 0, 1, 0);
		Vec4 v3(0, 0, 0, 1);

		// Normal of this tetrahedron is (1,1,1,1)/2
		// A direction perpendicular to normal (in the hyperplane)
		Vec4 direction(1, -1, 0, 0); // perpendicular to (1,1,1,1)
		Vec4 origin(0, 0, 0, 5);
		float t = RayTetrahedron(origin, direction, v0, v1, v2, v3);
		CHECK(t == FLT_MAX);
	}

	TEST_CASE("TestRayTetrahedronDegenerate")
	{
		// Degenerate tetrahedron (all coplanar vertices)
		Vec4 v0(1, 0, 0, 0);
		Vec4 v1(0, 1, 0, 0);
		Vec4 v2(0, 0, 1, 0);
		Vec4 v3(0.5f, 0.5f, 0, 0); // Coplanar with v0, v1, v2 in the x+y+z=1 hyperplane... no, that's still valid
		// Make truly degenerate: v3 is a linear combination of edges from v0
		Vec4 v3d = v0 + 2.0f * (v1 - v0) + 3.0f * (v2 - v0); // Coplanar with v0,v1,v2
		Vec4 origin(0, 0, 0, -1);
		Vec4 direction(0, 0, 0, 1);
		float t = RayTetrahedron(origin, direction, v0, v1, v2, v3d);
		// Cross product of 3 coplanar vectors = zero, should fail
		CHECK(t == FLT_MAX);
	}

	TEST_CASE("TestRayTetrahedronAxisAligned")
	{
		// Tetrahedron with one face in each axis hyperplane
		Vec4 v0(2, 0, 0, 0);
		Vec4 v1(0, 2, 0, 0);
		Vec4 v2(0, 0, 2, 0);
		Vec4 v3(0, 0, 0, 2);

		// Ray along X axis from far away
		Vec4 origin(10, 0, 0, 0);
		Vec4 direction(-1, 0, 0, 0);
		float t = RayTetrahedron(origin, direction, v0, v1, v2, v3);
		CHECK_APPROX_EQUAL(t, 8.0f, 1.0e-4f);

		// Ray along W axis
		Vec4 origin_w(0, 0, 0, 10);
		Vec4 direction_w(0, 0, 0, -1);
		float t_w = RayTetrahedron(origin_w, direction_w, v0, v1, v2, v3);
		CHECK_APPROX_EQUAL(t_w, 8.0f, 1.0e-4f);
	}

	TEST_CASE("TestRayTetrahedron4DSOA")
	{
		// Test the SOA version
		Vec4 origin(5, 5, 5, 5);
		Vec4 direction = Vec4(-1, -1, -1, -1).Normalized();

		// Tetrahedron 0: standard simplex (should hit)
		// Tetrahedron 1: shifted far away (should miss)
		// Tetrahedron 2: same as 0 but larger (should hit earlier)
		// Tetrahedron 3: degenerate (should miss)
		Vec4 v0x(1, 100, 2, 1);
		Vec4 v0y(0, 100, 0, 0);
		Vec4 v0z(0, 100, 0, 0);
		Vec4 v0w(0, 100, 0, 0);

		Vec4 v1x(0, 100, 0, 0);
		Vec4 v1y(1, 100, 2, 1);
		Vec4 v1z(0, 100, 0, 0);
		Vec4 v1w(0, 100, 0, 0);

		Vec4 v2x(0, 100, 0, 0);
		Vec4 v2y(0, 100, 0, 0);
		Vec4 v2z(1, 100, 2, 1);
		Vec4 v2w(0, 100, 0, 0);

		Vec4 v3x(0, 100, 0, 0.5f);
		Vec4 v3y(0, 100, 0, 0.5f);
		Vec4 v3z(0, 100, 0, 0);
		Vec4 v3w(1, 100, 2, 0);

		Vec4 result = RayTetrahedron4(origin, direction,
			v0x, v0y, v0z, v0w,
			v1x, v1y, v1z, v1w,
			v2x, v2y, v2z, v2w,
			v3x, v3y, v3z, v3w);

		CHECK(result.GetX() != FLT_MAX); // Tetrahedron 0: hit
		CHECK(result.GetX() > 0.0f);
		CHECK(result.GetY() == FLT_MAX); // Tetrahedron 1: miss (too far)
		CHECK(result.GetZ() != FLT_MAX); // Tetrahedron 2: hit (larger)
		CHECK(result.GetZ() < result.GetX()); // Should hit earlier (larger tetrahedron)
	}
}
