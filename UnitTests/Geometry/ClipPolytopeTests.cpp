// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/ClipPoly.h>

JPH_SUPPRESS_WARNINGS_STD_BEGIN
#include <random>
JPH_SUPPRESS_WARNINGS_STD_END

// TODO(4D): tests for the 4D contact-patch clip (conceptual gap 1). In 3D the patch is a 2D polygon
// and can be clipped as an ordered loop; in 4D it is a 3D polytope given as an unordered vertex set.
// These check the clip against independently computed answers rather than against itself.
TEST_SUITE("ClipPolytopeTests")
{
	using Face = StaticArray<Vec4, 32>;
	using Out = StaticArray<Vec4, 64>;

	/// The 8 vertices of an axis aligned cube cell lying in the w = inW hyperplane
	static Face sCubeCell(Vec4Arg inCenter, float inHalf, float inW)
	{
		Face f;
		for (int mask = 0; mask < 8; ++mask)
			f.push_back(Vec4(inCenter.GetX() + ((mask & 1)? inHalf : -inHalf),
							 inCenter.GetY() + ((mask & 2)? inHalf : -inHalf),
							 inCenter.GetZ() + ((mask & 4)? inHalf : -inHalf),
							 inW));
		return f;
	}

	/// Is inPoint inside the axis aligned cube (ignoring w)?
	static bool sInCube(Vec4Arg inPoint, Vec4Arg inCenter, float inHalf, float inTol)
	{
		for (int i = 0; i < 3; ++i)
			if (abs(inPoint[i] - inCenter[i]) > inHalf + inTol)
				return false;
		return true;
	}

	/// Axis aligned bounds of a point set over x,y,z
	static void sBounds(const Out &inPoints, Vec4 &outMin, Vec4 &outMax)
	{
		outMin = Vec4(FLT_MAX, FLT_MAX, FLT_MAX, 0);
		outMax = Vec4(-FLT_MAX, -FLT_MAX, -FLT_MAX, 0);
		for (const Vec4 &p : inPoints)
			for (int i = 0; i < 3; ++i)
			{
				outMin.SetX(i == 0? min(outMin.GetX(), p.GetX()) : outMin.GetX());
				outMin.SetY(i == 1? min(outMin.GetY(), p.GetY()) : outMin.GetY());
				outMin.SetZ(i == 2? min(outMin.GetZ(), p.GetZ()) : outMin.GetZ());
				outMax.SetX(i == 0? max(outMax.GetX(), p.GetX()) : outMax.GetX());
				outMax.SetY(i == 1? max(outMax.GetY(), p.GetY()) : outMax.GetY());
				outMax.SetZ(i == 2? max(outMax.GetZ(), p.GetZ()) : outMax.GetZ());
			}
	}

	/// A cube cell has 6 facets, so the halfspace derivation must find exactly 6 planes
	TEST_CASE("TestClipHalfSpacesOfCubeCell")
	{
		Face cube = sCubeCell(Vec4::sZero(), 1.0f, 0.0f);
		ClipHalfSpaces planes;
		GetClipHalfSpaces(cube, Vec4(0, 0, 0, 1), planes);
		CHECK(planes.size() == 6);

		// Every normal must be perpendicular to the axis and axis aligned in x/y/z
		for (const ClipHalfSpace &h : planes)
		{
			CHECK_APPROX_EQUAL(h.mNormal.GetW(), 0.0f, 1.0e-5f);
			CHECK_APPROX_EQUAL(h.mNormal.Length(), 1.0f, 1.0e-5f);
			CHECK_APPROX_EQUAL(abs(h.mDistance), 1.0f, 1.0e-4f);
		}

		// The centre is inside, a point well outside is not
		CHECK(IsInsideClipHalfSpaces(Vec4(0, 0, 0, 0), planes, 1.0e-4f));
		CHECK(IsInsideClipHalfSpaces(Vec4(0.9f, -0.9f, 0.9f, 5.0f), planes, 1.0e-4f));	// w is ignored
		CHECK_FALSE(IsInsideClipHalfSpaces(Vec4(1.5f, 0, 0, 0), planes, 1.0e-4f));
		CHECK_FALSE(IsInsideClipHalfSpaces(Vec4(0, 0, -1.5f, 0), planes, 1.0e-4f));
	}

	/// A small cell fully inside a big one must clip to exactly the small one
	TEST_CASE("TestClipPolytopeFullyContained")
	{
		Face small_cell = sCubeCell(Vec4::sZero(), 0.5f, 0.0f);
		Face big = sCubeCell(Vec4::sZero(), 4.0f, -1.0f);
		Out clipped;
		ClipPolytopeVsPolytope(small_cell, big, Vec4(0, 0, 0, 1), clipped);

		// All 8 corners survive, and they are lifted onto the small cell's hyperplane (w = 0)
		CHECK(clipped.size() == 8);
		for (const Vec4 &p : clipped)
		{
			CHECK(sInCube(p, Vec4::sZero(), 0.5f, 1.0e-4f));
			CHECK_APPROX_EQUAL(p.GetW(), 0.0f, 1.0e-4f);
		}
	}

	/// Two cells overlapping on half their extent: the patch must be the overlap region, which has
	/// vertices that are NOT vertices of either input - they come from edge/facet crossings.
	TEST_CASE("TestClipPolytopePartialOverlap")
	{
		// A spans x in [-1,1], B spans x in [0,2]; overlap is x in [0,1]
		Face a = sCubeCell(Vec4::sZero(), 1.0f, 0.0f);
		Face b = sCubeCell(Vec4(1.0f, 0, 0, 0), 1.0f, -0.5f);
		Out clipped;
		ClipPolytopeVsPolytope(a, b, Vec4(0, 0, 0, 1), clipped);

		CHECK(clipped.size() >= 8);

		Vec4 bmin, bmax;
		sBounds(clipped, bmin, bmax);
		// The clipped patch must be exactly the overlap box [0,1] x [-1,1] x [-1,1]
		CHECK_APPROX_EQUAL(bmin.GetX(), 0.0f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmax.GetX(), 1.0f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmin.GetY(), -1.0f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmax.GetY(), 1.0f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmin.GetZ(), -1.0f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmax.GetZ(), 1.0f, 1.0e-3f);

		// Every point must be in both polytopes (projected along the axis) and on A's hyperplane
		for (const Vec4 &p : clipped)
		{
			CHECK(sInCube(p, Vec4::sZero(), 1.0f, 1.0e-3f));
			CHECK(sInCube(p, Vec4(1.0f, 0, 0, 0), 1.0f, 1.0e-3f));
			CHECK_APPROX_EQUAL(p.GetW(), 0.0f, 1.0e-4f);
		}
	}

	/// A "plus sign" overlap: neither polytope has a vertex inside the other, so the patch is made
	/// entirely of crossings. This is the case a vertices-only clip would get completely wrong.
	TEST_CASE("TestClipPolytopeCrossOverlapHasNoSharedVertices")
	{
		// A is long in x and thin in y; B is long in y and thin in x. Overlap is the small centre box.
		Face a, b;
		for (int mask = 0; mask < 8; ++mask)
		{
			a.push_back(Vec4((mask & 1)? 3.0f : -3.0f, (mask & 2)? 0.5f : -0.5f, (mask & 4)? 1.0f : -1.0f, 0.0f));
			b.push_back(Vec4((mask & 1)? 0.5f : -0.5f, (mask & 2)? 3.0f : -3.0f, (mask & 4)? 1.0f : -1.0f, 0.7f));
		}
		Out clipped;
		ClipPolytopeVsPolytope(a, b, Vec4(0, 0, 0, 1), clipped);

		CHECK(clipped.size() >= 8);
		Vec4 bmin, bmax;
		sBounds(clipped, bmin, bmax);
		CHECK_APPROX_EQUAL(bmin.GetX(), -0.5f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmax.GetX(), 0.5f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmin.GetY(), -0.5f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmax.GetY(), 0.5f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmin.GetZ(), -1.0f, 1.0e-3f);
		CHECK_APPROX_EQUAL(bmax.GetZ(), 1.0f, 1.0e-3f);
	}

	/// Disjoint polytopes must produce nothing
	TEST_CASE("TestClipPolytopeDisjoint")
	{
		Face a = sCubeCell(Vec4::sZero(), 1.0f, 0.0f);
		Face b = sCubeCell(Vec4(10.0f, 0, 0, 0), 1.0f, 0.0f);
		Out clipped;
		ClipPolytopeVsPolytope(a, b, Vec4(0, 0, 0, 1), clipped);
		CHECK(clipped.empty());
	}

	/// The result must not depend on which polytope is passed first (up to which hyperplane the
	/// points are lifted onto), and must be invariant to the order of the vertices in the array.
	TEST_CASE("TestClipPolytopeOrderInvariance")
	{
		Face a = sCubeCell(Vec4::sZero(), 1.0f, 0.0f);
		Face b = sCubeCell(Vec4(0.75f, 0.25f, 0, 0), 1.0f, -0.5f);

		Out forward;
		ClipPolytopeVsPolytope(a, b, Vec4(0, 0, 0, 1), forward);
		Vec4 fmin, fmax;
		sBounds(forward, fmin, fmax);

		// Shuffle a's vertices
		Face shuffled;
		for (int i = 7; i >= 0; --i)
			shuffled.push_back(a[i]);
		Out shuffled_out;
		ClipPolytopeVsPolytope(shuffled, b, Vec4(0, 0, 0, 1), shuffled_out);
		Vec4 smin, smax;
		sBounds(shuffled_out, smin, smax);
		CHECK_APPROX_EQUAL(fmin, smin, 1.0e-3f);
		CHECK_APPROX_EQUAL(fmax, smax, 1.0e-3f);

		// Swapping the roles must give the same patch in x/y/z
		Out reverse;
		ClipPolytopeVsPolytope(b, a, Vec4(0, 0, 0, 1), reverse);
		Vec4 rmin, rmax;
		sBounds(reverse, rmin, rmax);
		CHECK_APPROX_EQUAL(fmin.GetX(), rmin.GetX(), 1.0e-3f);
		CHECK_APPROX_EQUAL(fmax.GetX(), rmax.GetX(), 1.0e-3f);
		CHECK_APPROX_EQUAL(fmin.GetY(), rmin.GetY(), 1.0e-3f);
		CHECK_APPROX_EQUAL(fmax.GetY(), rmax.GetY(), 1.0e-3f);
	}

	/// A tetrahedral cell (what a ConvexHullShape returns) clipped by a big cube
	TEST_CASE("TestClipPolytopeTetraCell")
	{
		Face tetra;
		tetra.push_back(Vec4(0, 0, 0, 0));
		tetra.push_back(Vec4(1, 0, 0, 0));
		tetra.push_back(Vec4(0, 1, 0, 0));
		tetra.push_back(Vec4(0, 0, 1, 0));
		Face big = sCubeCell(Vec4::sZero(), 5.0f, -2.0f);

		ClipHalfSpaces planes;
		GetClipHalfSpaces(tetra, Vec4(0, 0, 0, 1), planes);
		CHECK(planes.size() == 4);		// A tetrahedron has 4 facets

		Out clipped;
		ClipPolytopeVsPolytope(tetra, big, Vec4(0, 0, 0, 1), clipped);

		// The tetra is fully contained, so all 4 of its corners must survive
		for (int i = 0; i < 4; ++i)
		{
			bool found = false;
			for (const Vec4 &p : clipped)
				if ((p - tetra[i]).LengthSq() < 1.0e-6f)
					found = true;
			CHECK(found);
		}

		// Every output point must lie in the intersection. Note the output may legitimately contain
		// points that are on the patch boundary without being vertices of it - a chord of the
		// clipper can cross a facet of the clipped polytope - which is why we assert membership
		// rather than an exact count. PruneContactPoints selects extremes, so these cost nothing.
		for (const Vec4 &p : clipped)
			CHECK(IsInsideClipHalfSpaces(p, planes, 1.0e-3f));
		CHECK(clipped.size() <= 16);
	}
}
