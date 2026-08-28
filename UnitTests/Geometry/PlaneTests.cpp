// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/Plane.h>

TEST_SUITE("PlaneTests")
{
	TEST_CASE("TestPlaneSignedDistance")
	{
		// Hyperplane with normal along Y axis, passing through Y=2
		Plane p = Plane::sFromPointAndNormal(Vec4(0, 2, 0, 0), Vec4(0, 1, 0, 0));
		CHECK_APPROX_EQUAL(p.SignedDistance(Vec4(5, 7, 0, 0)), 5.0f);
		CHECK_APPROX_EQUAL(p.SignedDistance(Vec4(5, -3, 0, 0)), -5.0f);
		// W component should not affect distance for Y-normal plane
		CHECK_APPROX_EQUAL(p.SignedDistance(Vec4(5, 7, 0, 99.0f)), 5.0f);
	}

	TEST_CASE("TestPlaneSignedDistance4D")
	{
		// Hyperplane with normal along W axis
		Plane p = Plane::sFromPointAndNormal(Vec4(0, 0, 0, 3), Vec4(0, 0, 0, 1));
		CHECK_APPROX_EQUAL(p.SignedDistance(Vec4(0, 0, 0, 5)), 2.0f);
		CHECK_APPROX_EQUAL(p.SignedDistance(Vec4(0, 0, 0, 1)), -2.0f);
		CHECK_APPROX_EQUAL(p.SignedDistance(Vec4(99, 99, 99, 3)), 0.0f);
	}

	TEST_CASE("TestPlaneSignedDistanceDiagonal")
	{
		// Hyperplane with diagonal normal in 4D
		Vec4 normal = Vec4(1, 1, 1, 1).Normalized();
		Plane p = Plane::sFromPointAndNormal(Vec4::sZero(), normal);
		// Distance from origin to plane through origin = 0
		CHECK_APPROX_EQUAL(p.SignedDistance(Vec4::sZero()), 0.0f);
		// Distance from (1,1,1,1) = dot((1,1,1,1), normal) = 4/2 = 2
		CHECK_APPROX_EQUAL(p.SignedDistance(Vec4(1, 1, 1, 1)), 2.0f);
	}

	TEST_CASE("TestPlaneGetTransformed")
	{
		// Rotation matrix (pure rotation, no translation in Mat44)
		Mat44 rotation = Mat44::sIdentity(); // Use identity for simplicity
		Vec4 translation(5.0f, -7.0f, 9.0f, 1.0f);

		Vec4 point(11.0f, 13.0f, 15.0f, 2.0f);
		Vec4 normal = Vec4(1, 1, 0, 0).Normalized();

		Plane p1 = Plane::sFromPointAndNormal(point, normal).GetTransformed(rotation, translation);
		// With identity rotation: transformed normal = normal
		// transformed constant = original_constant - translation.Dot(normal)
		Vec4 transformed_point = Vec4(rotation * Lane4(point)) + translation;
		Plane p2 = Plane::sFromPointAndNormal(transformed_point, Vec4(rotation * Lane4(normal)));

		CHECK_APPROX_EQUAL(p1.GetNormal(), p2.GetNormal());
		CHECK_APPROX_EQUAL(p1.GetConstant(), p2.GetConstant(), 1.0e-4f);
	}

	TEST_CASE("TestPlaneFromPointsCCW")
	{
		// 4 points defining a hyperplane: the XYZ hyperplane (W=0)
		Vec4 v1(1, 0, 0, 0);
		Vec4 v2(0, 1, 0, 0);
		Vec4 v3(0, 0, 1, 0);
		Vec4 v4(1, 1, 0, 0);

		Plane p = Plane::sFromPointsCCW(v1, v2, v3, v4);
		// Normal should be along W axis (perpendicular to XYZ space)
		Vec4 n = p.GetNormal();
		// Check that the normal is perpendicular to all edge vectors
		CHECK_APPROX_EQUAL(n.Dot(v2 - v1), 0.0f, 1.0e-5f);
		CHECK_APPROX_EQUAL(n.Dot(v3 - v1), 0.0f, 1.0e-5f);
		CHECK_APPROX_EQUAL(n.Dot(v4 - v1), 0.0f, 1.0e-5f);
		// Normal should be unit length
		CHECK_APPROX_EQUAL(n.Length(), 1.0f, 1.0e-5f);
		// All 4 points should have signed distance 0
		CHECK_APPROX_EQUAL(p.SignedDistance(v1), 0.0f, 1.0e-5f);
		CHECK_APPROX_EQUAL(p.SignedDistance(v2), 0.0f, 1.0e-5f);
		CHECK_APPROX_EQUAL(p.SignedDistance(v3), 0.0f, 1.0e-5f);
		CHECK_APPROX_EQUAL(p.SignedDistance(v4), 0.0f, 1.0e-5f);
	}

	TEST_CASE("TestPlaneIntersectPlanes")
	{
		// 4 axis-aligned hyperplanes intersecting at a point
		Plane p1 = Plane::sFromPointAndNormal(Vec4(3, 0, 0, 0), Vec4(1, 0, 0, 0));
		Plane p2 = Plane::sFromPointAndNormal(Vec4(0, 2, 0, 0), Vec4(0, 1, 0, 0));
		Plane p3 = Plane::sFromPointAndNormal(Vec4(0, 0, 4, 0), Vec4(0, 0, 1, 0));
		Plane p4 = Plane::sFromPointAndNormal(Vec4(0, 0, 0, 5), Vec4(0, 0, 0, 1));

		{
			Vec4 point;
			bool found = Plane::sIntersectPlanes(p1, p2, p3, p4, point);
			CHECK(found);
			CHECK_APPROX_EQUAL(point, Vec4(3, 2, 4, 5), 1.0e-5f);
		}

		{
			// Two parallel planes should fail
			Plane p5 = Plane::sFromPointAndNormal(Vec4(0, 3, 0, 0), Vec4(0, 1, 0, 0));
			Vec4 point;
			bool found = Plane::sIntersectPlanes(p1, p2, p5, p4, point);
			CHECK_FALSE(found);
		}
	}

	TEST_CASE("TestPlaneIntersectPlanesDiagonal")
	{
		// 4 non-axis-aligned hyperplanes
		// We know the answer should be (1, 2, 3, 4)
		// Note: the normals must be linearly independent or there is no unique intersection point.
		// The obvious cyclic choice (1,1,0,0), (0,1,1,0), (0,0,1,1), (1,0,0,1) is NOT usable: it
		// satisfies n1 - n2 + n3 - n4 = 0, so those 4 hyperplanes meet in a line, not a point.
		Vec4 target(1, 2, 3, 4);
		Vec4 n1 = Vec4(1, 1, 0, 0).Normalized();
		Vec4 n2 = Vec4(0, 1, 1, 0).Normalized();
		Vec4 n3 = Vec4(0, 0, 1, 1).Normalized();
		Vec4 n4 = Vec4(1, 0, 0, -1).Normalized();

		Plane p1 = Plane::sFromPointAndNormal(target, n1);
		Plane p2 = Plane::sFromPointAndNormal(target, n2);
		Plane p3 = Plane::sFromPointAndNormal(target, n3);
		Plane p4 = Plane::sFromPointAndNormal(target, n4);

		Vec4 point;
		bool found = Plane::sIntersectPlanes(p1, p2, p3, p4, point);
		CHECK(found);
		CHECK_APPROX_EQUAL(point, target, 1.0e-4f);
	}

	TEST_CASE("TestPlaneProjectPointOnPlane")
	{
		Vec4 normal = Vec4(0, 0, 0, 1);
		Plane p = Plane::sFromPointAndNormal(Vec4(0, 0, 0, 2), normal);

		Vec4 projected = p.ProjectPointOnPlane(Vec4(1, 2, 3, 5));
		CHECK_APPROX_EQUAL(projected, Vec4(1, 2, 3, 2), 1.0e-5f);
	}
}
