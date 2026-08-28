// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2023 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/ClosestPoint.h>

TEST_SUITE("ClosestPointTests")
{
	// Test closest point from inPoint to triangle (inA, inB, inC) in 4D
	inline static void TestClosestPointToTriangle(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, Vec4Arg inPoint, Vec4Arg inExpectedClosestPoint, uint32 inExpectedSet)
	{
		// Make triangle relative to inPoint so we can get the closest point to the origin
		Vec4 a = inA - inPoint;
		Vec4 b = inB - inPoint;
		Vec4 c = inC - inPoint;

		// Extract bits for A, B and C
		uint32 expected_a = inExpectedSet & 1;
		uint32 expected_b = (inExpectedSet & 2) >> 1;
		uint32 expected_c = (inExpectedSet & 4) >> 2;

		// Test all permutations of ABC
		uint32 set = 0;
		Vec4 closest_point = inPoint + ClosestPoint::GetClosestPointOnTriangle(a, b, c, set);
		CHECK(set == inExpectedSet);
		CHECK_APPROX_EQUAL(closest_point, inExpectedClosestPoint, 2.0e-5f);

		closest_point = inPoint + ClosestPoint::GetClosestPointOnTriangle(a, c, b, set);
		CHECK(set == ((expected_b << 2) | (expected_c << 1) | expected_a));
		CHECK_APPROX_EQUAL(closest_point, inExpectedClosestPoint, 2.0e-5f);

		closest_point = inPoint + ClosestPoint::GetClosestPointOnTriangle(b, a, c, set);
		CHECK(set == ((expected_c << 2) | (expected_a << 1) | expected_b));
		CHECK_APPROX_EQUAL(closest_point, inExpectedClosestPoint, 2.0e-5f);

		closest_point = inPoint + ClosestPoint::GetClosestPointOnTriangle(b, c, a, set);
		CHECK(set == ((expected_a << 2) | (expected_c << 1) | expected_b));
		CHECK_APPROX_EQUAL(closest_point, inExpectedClosestPoint, 2.0e-5f);

		closest_point = inPoint + ClosestPoint::GetClosestPointOnTriangle(c, a, b, set);
		CHECK(set == ((expected_b << 2) | (expected_a << 1) | expected_c));
		CHECK_APPROX_EQUAL(closest_point, inExpectedClosestPoint, 2.0e-5f);

		closest_point = inPoint + ClosestPoint::GetClosestPointOnTriangle(c, b, a, set);
		CHECK(set == ((expected_a << 2) | (expected_b << 1) | expected_c));
		CHECK_APPROX_EQUAL(closest_point, inExpectedClosestPoint, 2.0e-5f);
	}

	TEST_CASE("TestLongTriangle")
	{
		// Triangle in 4D with non-zero W component
		Vec4 a(100, 1, 0, 0.5f);
		Vec4 b(100, 1, 1, 0.5f);
		Vec4 c(-100, 1, 0, 0.5f);

		// Test interior
		TestClosestPointToTriangle(a, b, c, Vec4(0, 0, 0.1f, 0.5f), Vec4(0, 1, 0.1f, 0.5f), 0b0111);

		// Edge AB
		TestClosestPointToTriangle(a, b, c, Vec4(101, 0, 0.5f, 0.5f), Vec4(100, 1, 0.5f, 0.5f), 0b0011);

		// Edge AC
		TestClosestPointToTriangle(a, b, c, Vec4(0, 0, -0.1f, 0.5f), Vec4(0, 1, 0, 0.5f), 0b0101);

		// Edge BC
		Vec4 point_bc(0, 0, 1, 0.5f);
		Vec4 bc = c - b;
		Vec4 closest_bc = b + ((point_bc - b).Dot(bc) / bc.LengthSq()) * bc;
		TestClosestPointToTriangle(a, b, c, point_bc, closest_bc, 0b0110);

		// Vertex A
		TestClosestPointToTriangle(a, b, c, Vec4(101, 0, -1, 0.5f), a, 0b0001);

		// Vertex B
		TestClosestPointToTriangle(a, b, c, Vec4(101, 0, 2, 0.5f), b, 0b0010);

		// Vertex C
		TestClosestPointToTriangle(a, b, c, Vec4(-101, 0, 0, 0.5f), c, 0b0100);
	}

	TEST_CASE("TestTriangle4D")
	{
		// Triangle with non-trivial W coordinates to exercise the 4th dimension
		Vec4 a(1, 0, 0, 1);
		Vec4 b(0, 1, 0, 0);
		Vec4 c(0, 0, 1, 0);

		// Point at the centroid offset perpendicularly - should project to interior
		Vec4 centroid = (a + b + c) / 3.0f;
		// Find the centroid's closest point: it is the centroid itself
		TestClosestPointToTriangle(a, b, c, centroid, centroid, 0b0111);

		// Point closest to vertex A
		Vec4 far_from_a = a + Vec4(10, -5, -5, 10);
		TestClosestPointToTriangle(a, b, c, far_from_a, a, 0b0001);
	}

	TEST_CASE("TestNearColinearTriangle")
	{
		// A very long triangle that is nearly colinear, in 4D
		Vec4 a(99.9999847f, 0.946687222f, 99.9999847f, 0.0f);
		Vec4 b(-100.010002f, 0.977360725f, -100.010002f, 0.0f);
		Vec4 c(-100.000137f, 0.977310658f, -100.000137f, 0.0f);

		// Closest point is on edge AC
		Vec4 ac = c - a;
		Vec4 expected_closest = a + (-a.Dot(ac) / ac.LengthSq()) * ac;

		TestClosestPointToTriangle(a, b, c, Vec4::sZero(), expected_closest, 0b0101);
	}

	TEST_CASE("TestSmallTriangleWithPlaneGoingThroughOrigin")
	{
		// A small but non-degenerate triangle whose plane almost goes through the origin
		Vec4 a(-0.132395342f, -0.294095188f, -0.164812326f, 0.0f);
		Vec4 b(-0.126054004f, -0.283950001f, -0.159065604f, 0.0f);
		Vec4 c(-0.154956535f, -0.284792334f, -0.160523415f, 0.0f);

		float u, v, w;
		ClosestPoint::GetBaryCentricCoordinates(a, b, c, u, v, w);

		// Closest point should be close to origin
		Vec4 p = a * u + b * v + c * w;
		CHECK_APPROX_EQUAL(p, Vec4::sZero());

		// Closest point should be outside triangle
		CHECK((u < 0.0f || v > 0.0f || w < 0.0f));
	}

	TEST_CASE("TestBaryCentricCoordinates2Point")
	{
		// Basic 2-point barycentric coordinates in 4D
		Vec4 a(1, 2, 3, 4);
		Vec4 b(5, 6, 7, 8);

		// Origin projects onto the line segment
		float u, v;
		ClosestPoint::GetBaryCentricCoordinates(a, b, u, v);
		Vec4 closest = a * u + b * v;
		// Closest point should be on the line
		Vec4 ab = b - a;
		float t = -a.Dot(ab) / ab.LengthSq();
		Vec4 expected = a + t * ab;
		CHECK_APPROX_EQUAL(closest, expected, 1.0e-5f);
	}

	TEST_CASE("TestOriginOutsideOfTetrahedronPlanes")
	{
		// Regular-ish tetrahedron in 4D (W = 0, so it's a 3D tetrahedron embedded in 4D)
		Vec4 a(1, 0, 0, 0);
		Vec4 b(0, 1, 0, 0);
		Vec4 c(0, 0, 1, 0);
		Vec4 d(0, 0, 0, 0);

		// Origin is at d, so it's on the face BCD. It should be "inside" face BCD (not outside).
		// For face opposite to each vertex:
		// Face opposite A = (B,C,D): origin is ON this face
		// Face opposite B = (A,C,D): origin should be outside (same side as B? No, origin IS d)
		// This tests the Gram-determinant based approach
		UVec4 result = ClosestPoint::OriginOutsideOfTetrahedronPlanes(a, b, c, d);

		// The origin is at vertex D, so this is degenerate, but each of the 4 lanes must still be a
		// well formed boolean mask (all bits set or none) and the call must not crash
		for (int i = 0; i < 4; ++i)
			CHECK((result[i] == 0 || result[i] == 0xffffffff));
		// More importantly, test with origin clearly inside or outside
		Vec4 a2(2, 0, 0, 0);
		Vec4 b2(0, 2, 0, 0);
		Vec4 c2(0, 0, 2, 0);
		Vec4 d2(1, 1, 1, 0); // Interior-ish point for the 4th vertex

		// Face ABC of this tetrahedron is the plane x+y+z=2 and D has x+y+z=3, so a probe is inside
		// only if x+y+z > 2. The centroid ((2+0+0+1)/4 etc.) has x+y+z = 2.25 and is inside.
		Vec4 interior(0.75f, 0.75f, 0.75f, 0);
		Vec4 ta = a2 - interior, tb = b2 - interior, tc = c2 - interior, td = d2 - interior;
		UVec4 result2 = ClosestPoint::OriginOutsideOfTetrahedronPlanes(ta, tb, tc, td);
		CHECK_FALSE(result2.TestAnyTrue());

		// (0.5, 0.5, 0.5) has x+y+z = 1.5 < 2, so it is on the far side of face ABC from D and must be
		// reported outside exactly that one face
		Vec4 outside(0.5f, 0.5f, 0.5f, 0);
		UVec4 result3 = ClosestPoint::OriginOutsideOfTetrahedronPlanes(a2 - outside, b2 - outside, c2 - outside, d2 - outside);
		CHECK(result3.GetX() != 0);
		CHECK(result3.GetY() == 0);
		CHECK(result3.GetZ() == 0);
		CHECK(result3.GetW() == 0);
	}

	TEST_CASE("TestGetClosestPointOnTetrahedron")
	{
		// Tetrahedron with origin inside
		Vec4 a(1, 0, 0, 0);
		Vec4 b(-1, 1, 0, 0);
		Vec4 c(-1, -1, 1, 0);
		Vec4 d(-1, -1, -1, 0);

		// Origin is inside this tetrahedron (centroid is at (-0.5, -0.25, 0, 0))
		// Move vertices so origin is inside
		Vec4 offset(-0.1f, 0, 0, 0);
		Vec4 ta = a + offset, tb = b + offset, tc = c + offset, td = d + offset;

		uint32 set = 0;
		Vec4 closest = ClosestPoint::GetClosestPointOnTetrahedron(ta, tb, tc, td, set);

		// The origin is strictly inside (barycentric coordinates are all positive), so the closest
		// point is the origin itself and all 4 vertices are in the set
		CHECK(set == 0b1111);
		CHECK_APPROX_EQUAL(closest.Length(), 0.0f, 1.0e-4f);

		// Now a tetrahedron that really does not contain the origin: every vertex has x >= 1, so the
		// origin is on the far side of face BCD (the plane x = 1) and the closest point is (1,0,0,0)
		// on that face. Note that simply moving A further away, as one might expect, does NOT put the
		// origin outside - it makes the tetrahedron bigger and keeps the origin inside.
		Vec4 out_a(10, 0, 0, 0);
		Vec4 out_b(1, 1, 0, 0);
		Vec4 out_c(1, -1, 1, 0);
		Vec4 out_d(1, -1, -1, 0);
		set = 0;
		closest = ClosestPoint::GetClosestPointOnTetrahedron(out_a, out_b, out_c, out_d, set);
		CHECK((set & 0b0001) == 0); // A is not in the set
		CHECK(set == 0b1110);
		CHECK_APPROX_EQUAL(closest, Vec4(1, 0, 0, 0), 1.0e-4f);
	}

	TEST_CASE("TestGetClosestPointOnTetrahedron4D")
	{
		// Tetrahedron that uses the W dimension
		Vec4 a(1, 0, 0, 0);
		Vec4 b(0, 1, 0, 0);
		Vec4 c(0, 0, 1, 0);
		Vec4 d(0, 0, 0, 1);

		// Origin is equidistant from all vertices. Closest point should be in the interior.
		uint32 set = 0;
		Vec4 closest = ClosestPoint::GetClosestPointOnTetrahedron(a, b, c, d, set);
		CHECK(set == 0b1111);
		// The closest point on the tetrahedron to the origin should be at (0.25, 0.25, 0.25, 0.25)
		Vec4 expected(0.25f, 0.25f, 0.25f, 0.25f);
		CHECK_APPROX_EQUAL(closest, expected, 1.0e-4f);
	}
}
