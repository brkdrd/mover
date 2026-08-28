// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/RaySphere.h>

TEST_SUITE("RaySphereTests")
{
	TEST_CASE("TestRaySphereHitXYZ")
	{
		// Unit sphere at origin
		Vec4 center = Vec4::sZero();
		float radius = 1.0f;

		// Ray along X axis from far away
		Vec4 origin(5, 0, 0, 0);
		Vec4 direction(-1, 0, 0, 0);
		float t = RaySphere(origin, direction, center, radius);
		CHECK_APPROX_EQUAL(t, 4.0f, 1.0e-5f);
	}

	TEST_CASE("TestRaySphereHitW")
	{
		// Sphere at the origin
		Vec4 center = Vec4::sZero();
		float radius = 1.0f;

		// Ray along W axis
		Vec4 origin(0, 0, 0, 5);
		Vec4 direction(0, 0, 0, -1);
		float t = RaySphere(origin, direction, center, radius);
		CHECK_APPROX_EQUAL(t, 4.0f, 1.0e-5f);
	}

	TEST_CASE("TestRaySphereHitDiagonal")
	{
		// Sphere at (1,1,1,1) with radius 2
		Vec4 center(1, 1, 1, 1);
		float radius = 2.0f;

		// Ray along the diagonal from far away
		Vec4 dir = Vec4(-1, -1, -1, -1).Normalized();
		Vec4 origin = center - 10.0f * dir; // far from center along +diagonal
		float t = RaySphere(origin, dir, center, radius);
		CHECK_APPROX_EQUAL(t, 8.0f, 1.0e-4f);
	}

	TEST_CASE("TestRaySphereMiss")
	{
		Vec4 center = Vec4::sZero();
		float radius = 1.0f;

		// Ray that passes near but misses the sphere
		Vec4 origin(0, 2, 0, 0);
		Vec4 direction(1, 0, 0, 0);
		float t = RaySphere(origin, direction, center, radius);
		CHECK(t == FLT_MAX);
	}

	TEST_CASE("TestRaySphereMiss4D")
	{
		// Sphere at origin with radius 1
		Vec4 center = Vec4::sZero();
		float radius = 1.0f;

		// Ray that would hit in 3D but misses in 4D because of W offset
		Vec4 origin(5, 0, 0, 1.5f); // W=1.5 puts it outside sphere
		Vec4 direction(-1, 0, 0, 0);
		float t = RaySphere(origin, direction, center, radius);
		// Perpendicular distance from ray to center = 1.5 > radius 1.0
		CHECK(t == FLT_MAX);
	}

	TEST_CASE("TestRaySphereInsideOrigin")
	{
		Vec4 center = Vec4::sZero();
		float radius = 2.0f;

		// Ray starting inside the sphere
		Vec4 origin(0.5f, 0, 0, 0);
		Vec4 direction(1, 0, 0, 0);
		float t = RaySphere(origin, direction, center, radius);
		// RaySphere's documented contract (same as 3D Jolt) is that a ray starting inside the sphere
		// returns fraction 0, not a negative fraction - the caller is already in contact.
		CHECK(t == 0.0f);
	}

	TEST_CASE("TestRaySphereBehind")
	{
		Vec4 center = Vec4::sZero();
		float radius = 1.0f;

		// Ray pointing away from sphere
		Vec4 origin(5, 0, 0, 0);
		Vec4 direction(1, 0, 0, 0); // pointing away
		float t = RaySphere(origin, direction, center, radius);
		CHECK(t == FLT_MAX);
	}

	TEST_CASE("TestRaySphereOffset")
	{
		// Sphere offset in W
		Vec4 center(0, 0, 0, 5);
		float radius = 1.0f;

		// Ray along W from origin
		Vec4 origin = Vec4::sZero();
		Vec4 direction(0, 0, 0, 1);
		float t = RaySphere(origin, direction, center, radius);
		CHECK_APPROX_EQUAL(t, 4.0f, 1.0e-5f);
	}

	TEST_CASE("TestRaySphereTangent")
	{
		Vec4 center = Vec4::sZero();
		float radius = 1.0f;

		// Ray tangent to sphere (grazing, distance = radius)
		Vec4 origin(0, 1, 0, 0);
		Vec4 direction(1, 0, 0, 0);
		float t = RaySphere(origin, direction, center, radius);
		// Should just barely hit (t = 0, tangent point at origin's X position)
		CHECK(t != FLT_MAX);
		CHECK_APPROX_EQUAL(t, 0.0f, 1.0e-4f);
	}
}
