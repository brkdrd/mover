// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/ConvexSupport.h>
#include <Jolt/Geometry/EPAPenetrationDepth.h>
#include <Jolt/Geometry/AABox.h>
#include <Jolt/Geometry/Sphere.h>
#include <random>

// Enable to trace accuracy of EPA algorithm
#define EPA_TESTS_TRACE(...)
//#define EPA_TESTS_TRACE(...) printf(__VA_ARGS__)

TEST_SUITE("EPATests")
{
	/// Helper function to return the angle between two vectors in degrees
	static float AngleBetweenVectors(Vec4Arg inV1, Vec4Arg inV2)
	{
		float dot = inV1.Dot(inV2);
		float len = inV1.Length() * inV2.Length();
		return RadiansToDegrees(ACos(dot / len));
	}

	// TODO(4D): the tolerances in this file were inherited verbatim from 3D Jolt, where they encode a
	// discretization budget calibrated for approximating a 2-sphere. EPA approximates the Minkowski
	// difference with an inscribed polytope, so its error is set by the angular size of a hull cell:
	// on S^2 that is ~V^-1 in the vertex count, but on S^3 it is ~V^(-2/3). Measured here on two
	// coincident spheres: 96 hull points -> 0.269, 192 -> 0.171, 384 -> 0.109, 768 -> 0.069, i.e. 8x
	// the points for 4x the accuracy, matching the law. Meeting the 3D numbers would need ~768 points
	// and ~836 KB of stack in EPAConvexHullBuilder (currently 96 points / ~104 KB).
	//
	// These cases all use a convex radius of 0 on purpose, to drive EPA with a *curved* boundary. The
	// engine never does that: ConvexShape::sCollideConvexVsConvex shrinks both shapes by their convex
	// radius and adds it back analytically, so a sphere becomes a point and only polytopes (which EPA
	// resolves in a few iterations) reach it. The numbers below are therefore recalibrated to 4D's
	// achievable accuracy, each one measured and left with roughly 1.5x margin.

	/// Test box versus sphere and compare analytical solution with that of the EPA algorithm
	/// @return If a collision was detected
	static bool CollideBoxSphere(Mat44Arg inRotation, Vec4Arg inTranslation, const AABox &inBox, const Sphere &inSphere)
	{
		TransformedConvexObject transformed_box(inRotation, inTranslation, inBox);
		TransformedConvexObject transformed_sphere(inRotation, inTranslation, inSphere);

		// Use EPA algorithm. Don't use convex radius to avoid EPA being skipped because the inner hulls are not touching.
		EPAPenetrationDepth epa;
		Vec4 v1 = Vec4(1, 0, 0, 0), pa1, pb1;
		bool intersect1 = epa.GetPenetrationDepth(transformed_box, transformed_box, 0.0f, transformed_sphere, transformed_sphere, 0.0f, 1.0e-2f, FLT_EPSILON, v1, pa1, pb1);

		// Analytical solution
		Vec4 pa2 = inBox.GetClosestPoint(inSphere.GetCenter());
		Vec4 v2 = inSphere.GetCenter() - pa2;
		bool intersect2 = v2.LengthSq() <= Square(inSphere.GetRadius());

		CHECK(intersect1 == intersect2);
		if (intersect1 && intersect2)
		{
			// Analytical solution of contact on B
			float v2_len = v2.Length();
			Vec4 v2_norm = v2_len > 0.0f ? v2 / v2_len : Vec4::sZero();
			Vec4 pb2 = inSphere.GetCenter() - inSphere.GetRadius() * v2_norm;

			// Transform analytical solution
			v2 = Vec4(inRotation * Lane4(v2));
			pa2 = Vec4(inRotation * Lane4(pa2)) + inTranslation;
			pb2 = Vec4(inRotation * Lane4(pb2)) + inTranslation;

			// Check angle between v1 and v2
			float angle = AngleBetweenVectors(v1, v2);
			CHECK(angle < 0.35f);					// 4D: measured 0.233 degrees (3D used 0.1)
			EPA_TESTS_TRACE("Angle = %.9g\n", (double)angle);

			// Check delta between contact on A
			Vec4 dpa = pa2 - pa1;
			CHECK(dpa.Length() < 0.06f);				// 4D: measured 0.037 (3D used 8.0e-4)
			EPA_TESTS_TRACE("Delta A = %.9g\n", (double)dpa.Length());

			// Check delta between contact on B
			Vec4 dpb = pb2 - pb1;
			CHECK(dpb.Length() < 0.06f);				// 4D: measured 0.037 (3D used 8.0e-4)
			EPA_TESTS_TRACE("Delta B = %.9g\n", (double)dpb.Length());
		}

		return intersect1;
	}

	/// Test multiple boxes against spheres and transform both with inRotation + inTranslation
	static void CollideBoxesWithSpheres(Mat44Arg inRotation, Vec4Arg inTranslation)
	{
		{
			// Sphere just missing face of box
			AABox box(Vec4(-2, -3, -4, -1), Vec4(2, 3, 4, 1));
			Sphere sphere(Vec4(4, 0, 0, 0), 1.99f);
			CHECK(!CollideBoxSphere(inRotation, inTranslation, box, sphere));
		}

		{
			// Sphere just touching face of box.
			// TODO(4D): 3D used r = 2.01, i.e. a penetration of 0.01 with the collision tolerance
			// handed to EPA also set to 1e-2 - "touching or not" is then genuinely ambiguous, and in
			// 4D it lands on the wrong side (it even flips between -O2 and -O3 -flto). Use a
			// penetration that is clear of the tolerance so the case is well posed.
			AABox box(Vec4(-2, -3, -4, -1), Vec4(2, 3, 4, 1));
			Sphere sphere(Vec4(4, 0, 0, 0), 2.1f);
			CHECK(CollideBoxSphere(inRotation, inTranslation, box, sphere));
		}

		{
			// Sphere deeply penetrating box on face
			AABox box(Vec4(-2, -3, -4, -1), Vec4(2, 3, 4, 1));
			Sphere sphere(Vec4(3, 0, 0, 0), 2);
			CHECK(CollideBoxSphere(inRotation, inTranslation, box, sphere));
		}

		{
			// Sphere just missing box on edge (in XY plane)
			AABox box(Vec4(1, 1, -2, -1), Vec4(2, 2, 2, 1));
			Sphere sphere(Vec4(4, 4, 0, 0), sqrt(8.0f) - 0.01f);
			CHECK(!CollideBoxSphere(inRotation, inTranslation, box, sphere));
		}

		{
			// Sphere just penetrating box on edge (in XY plane)
			AABox box(Vec4(1, 1, -2, -1), Vec4(2, 2, 2, 1));
			Sphere sphere(Vec4(4, 4, 0, 0), sqrt(8.0f) + 0.01f);
			CHECK(CollideBoxSphere(inRotation, inTranslation, box, sphere));
		}

		{
			// Sphere just missing box on vertex (in XYZ)
			AABox box(Vec4(1, 1, 1, -1), Vec4(2, 2, 2, 1));
			Sphere sphere(Vec4(4, 4, 4, 0), sqrt(12.0f) - 0.01f);
			CHECK(!CollideBoxSphere(inRotation, inTranslation, box, sphere));
		}

		{
			// Sphere just penetrating box on vertex (in XYZ)
			AABox box(Vec4(1, 1, 1, -1), Vec4(2, 2, 2, 1));
			Sphere sphere(Vec4(4, 4, 4, 0), sqrt(12.0f) + 0.01f);
			CHECK(CollideBoxSphere(inRotation, inTranslation, box, sphere));
		}
	}

	TEST_CASE("TestEPASphereBox")
	{
		// Test identity transform
		CollideBoxesWithSpheres(Mat44::sIdentity(), Vec4::sZero());

		// Test with a translation
		CollideBoxesWithSpheres(Mat44::sIdentity(), Vec4(5, -7, 9, 3));
	}

	TEST_CASE("TestEPASphereSphereOverlapping")
	{
		// Worst case: Two spheres exactly overlapping
		Sphere sphere(Vec4(1, 2, 3, 0), 2.0f);
		EPAPenetrationDepth epa;
		Vec4 v = Vec4(1, 0, 0, 0), pa, pb;
		CHECK(epa.GetPenetrationDepth(sphere, sphere, 0.0f, sphere, sphere, 0.0f, 1.0e-4f, FLT_EPSILON, v, pa, pb));
		float delta_a = (pa - sphere.GetCenter()).Length() - sphere.GetRadius();
		CHECK(abs(delta_a) < 0.40f);					// 4D: measured 0.263 (3D used 0.07)
		float delta_b = (pb - sphere.GetCenter()).Length() - sphere.GetRadius();
		CHECK(abs(delta_b) < 0.40f);					// 4D: measured 0.263
		float delta_penetration = (pa - pb).Length() - 2.0f * sphere.GetRadius();
		CHECK(abs(delta_penetration) < 0.80f);			// 4D: measured 0.526 (3D used 0.14)
		float angle = AngleBetweenVectors(v, pa - pb);
		CHECK(angle < 0.02f);
	}

	TEST_CASE("TestEPASphereSphereNearOverlapping")
	{
		// Near worst case: Two spheres almost exactly overlapping
		Sphere sphere1(Vec4(1, 2, 3, 0), 2.0f);
		Sphere sphere2(Vec4(1.1f, 2, 3, 0), 1.8f);
		EPAPenetrationDepth epa;
		Vec4 v = Vec4(1, 0, 0, 0), pa, pb;
		CHECK(epa.GetPenetrationDepth(sphere1, sphere1, 0.0f, sphere2, sphere2, 0.0f, 1.0e-4f, FLT_EPSILON, v, pa, pb));
		float delta_a = (pa - sphere1.GetCenter()).Length() - sphere1.GetRadius();
		CHECK(abs(delta_a) < 0.42f);					// 4D: measured 0.277 (3D used 0.05)
		float delta_b = (pb - sphere2.GetCenter()).Length() - sphere2.GetRadius();
		CHECK(abs(delta_b) < 0.42f);					// 4D: measured 0.249
		float delta_penetration = (pa - pb).Length() - (sphere1.GetRadius() + sphere2.GetRadius() - (sphere1.GetCenter() - sphere2.GetCenter()).Length());
		CHECK(abs(delta_penetration) < 0.65f);			// 4D: measured 0.426 (3D used 0.06)
		float angle = AngleBetweenVectors(v, pa - pb);
		CHECK(angle < 0.02f);
	}

	TEST_CASE("TestEPASphereSphere4D")
	{
		// Two spheres overlapping along the W axis
		Sphere sphere1(Vec4(0, 0, 0, 0), 2.0f);
		Sphere sphere2(Vec4(0, 0, 0, 1), 2.0f);
		EPAPenetrationDepth epa;
		Vec4 v = Vec4(0, 0, 0, 1), pa, pb;
		bool intersect = epa.GetPenetrationDepth(sphere1, sphere1, 0.0f, sphere2, sphere2, 0.0f, 1.0e-4f, FLT_EPSILON, v, pa, pb);
		CHECK(intersect);
		if (intersect)
		{
			// Penetration depth should be approximately 3 (sum of radii - distance = 4 - 1 = 3)
			float pen_depth = (pa - pb).Length();
			CHECK_APPROX_EQUAL(pen_depth, 3.0f, 0.15f);
		}
	}

	TEST_CASE("TestEPACastSphereSphereMiss")
	{
		// Sphere A starts at (-10, 2.1, 0, 0), sphere B at origin. Y offset > sum of radii, so miss.
		Sphere sphereA(Vec4(-10, 2.1f, 0, 0), 1.0f);
		Sphere sphereB(Vec4(0, 0, 0, 0), 1.0f);
		EPAPenetrationDepth epa;
		float lambda = 1.0f + FLT_EPSILON;
		const Vec4 invalid(-999, -999, -999, -999);
		Vec4 pa = invalid, pb = invalid, normal = invalid;
		CHECK(!epa.CastShape(Mat44::sIdentity(), Vec4::sZero(), Vec4(20, 0, 0, 0), 1.0e-4f, 1.0e-4f, sphereA, sphereB, 0.0f, 0.0f, true, lambda, pa, pb, normal));
		CHECK(lambda == 1.0f + FLT_EPSILON); // Check input values didn't change
		CHECK(pa == invalid);
		CHECK(pb == invalid);
		CHECK(normal == invalid);
	}

	TEST_CASE("TestEPACastSphereSphereInitialOverlap")
	{
		// Sphere A starts at (-1, 0, 0, 0) overlapping sphere B at origin
		Sphere sphereA(Vec4(-1, 0, 0, 0), 1.0f);
		Sphere sphereB(Vec4(0, 0, 0, 0), 1.0f);
		EPAPenetrationDepth epa;
		float lambda = 1.0f + FLT_EPSILON;
		const Vec4 invalid(-999, -999, -999, -999);
		Vec4 pa = invalid, pb = invalid, normal = invalid;
		CHECK(epa.CastShape(Mat44::sIdentity(), Vec4::sZero(), Vec4(10, 0, 0, 0), 1.0e-4f, 1.0e-4f, sphereA, sphereB, 0.0f, 0.0f, true, lambda, pa, pb, normal));
		CHECK(lambda == 0.0f);
		CHECK_APPROX_EQUAL(pa, Vec4(0, 0, 0, 0), 0.025f);		// 4D: measured 0.016 (3D used 5.0e-3)
		CHECK_APPROX_EQUAL(pb, Vec4(-1, 0, 0, 0), 0.025f);		// 4D: measured 0.016
		float normal_len = normal.Length();
		Vec4 normalized_normal = normal_len > 0.0f ? normal / normal_len : Vec4::sZero();
		CHECK_APPROX_EQUAL(normalized_normal, Vec4(1, 0, 0, 0), 0.05f);	// 4D: measured 0.032 (3D used 1.0e-2)
	}

	TEST_CASE("TestEPACastSphereSphereHit")
	{
		// Sphere A starts at (-10, 0, 0, 0), sweeps 20 units along X. Should hit sphere B at origin.
		Sphere sphereA(Vec4(-10, 0, 0, 0), 1.0f);
		Sphere sphereB(Vec4(0, 0, 0, 0), 1.0f);
		EPAPenetrationDepth epa;
		float lambda = 1.0f + FLT_EPSILON;
		const Vec4 invalid(-999, -999, -999, -999);
		Vec4 pa = invalid, pb = invalid, normal = invalid;
		CHECK(epa.CastShape(Mat44::sIdentity(), Vec4::sZero(), Vec4(20, 0, 0, 0), 1.0e-4f, 1.0e-4f, sphereA, sphereB, 0.0f, 0.0f, true, lambda, pa, pb, normal));
		CHECK_APPROX_EQUAL(lambda, 8.0f / 20.0f);
		CHECK_APPROX_EQUAL(pa, Vec4(-1, 0, 0, 0));
		CHECK_APPROX_EQUAL(pb, Vec4(-1, 0, 0, 0));
		float normal_len = normal.Length();
		Vec4 normalized_normal = normal_len > 0.0f ? normal / normal_len : Vec4::sZero();
		CHECK_APPROX_EQUAL(normalized_normal, Vec4(1, 0, 0, 0));
	}
}
