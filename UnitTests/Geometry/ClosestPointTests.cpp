// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2023 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/ClosestPoint.h>
#include <algorithm>
#include <random>

TEST_SUITE("ClosestPointTests")
{
	TEST_CASE("TestLowerDimensionalProjectionScaleInvariance")
	{
		// Exact dyadic geometry: a nonzero segment has no dimensionful minimum length;
		// the triangle projects to weights (1/2, 1/4, 1/4), independently of scale.
		for (float scale : { 0x1p-60f, 0x1p-30f, 0x1p-20f, 0x1p-10f, 1.0f, 0x1p10f, 0x1p60f })
			for (Vec4 translation : { Vec4::sZero(), Vec4(4, -8, 16, -32) * scale })
			{
				CAPTURE(scale);
				Vec4 offset = Vec4(0.25f, -0.125f, 0, 0) * scale;
				Vec4 direction = Vec4(1, 2, 3, 4) * scale;
				for (bool reverse : { false, true })
				{
					Vec4 a = translation + offset - direction, b = translation + offset + direction;
					if (reverse) std::swap(a, b);
					a -= translation;
					b -= translation;
					float u, v;
					CHECK(ClosestPoint::GetBaryCentricCoordinates(a, b, u, v));
					CHECK_APPROX_EQUAL(u, 0.5f, 1.0e-6f);
					CHECK_APPROX_EQUAL(v, 0.5f, 1.0e-6f);
					uint32 set;
					CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnLine(a, b, set) / scale, offset / scale, 2.0e-6f);
					CHECK(set == 3);
				}
				const Vec4 basis[] = { Vec4(-1, 0, 0.125f, 0.25f), Vec4(1, -1, 0.125f, 0.25f), Vec4(1, 1, 0.125f, 0.25f) };
				const float weights[] = { 0.5f, 0.25f, 0.25f };
				int permutation[] = { 0, 1, 2 };
				do
				{
					Vec4 a = (translation + basis[permutation[0]] * scale) - translation;
					Vec4 b = (translation + basis[permutation[1]] * scale) - translation;
					Vec4 c = (translation + basis[permutation[2]] * scale) - translation;
					float u, v, w;
					CHECK(ClosestPoint::GetBaryCentricCoordinates(a, b, c, u, v, w));
					CHECK_APPROX_EQUAL(u, weights[permutation[0]], 1.0e-6f);
					CHECK_APPROX_EQUAL(v, weights[permutation[1]], 1.0e-6f);
					CHECK_APPROX_EQUAL(w, weights[permutation[2]], 1.0e-6f);
					Vec4 expected = Vec4(0, 0, 0.125f, 0.25f) * scale;
					uint32 set;
					CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTriangle(a, b, c, set) / scale, expected / scale, 2.0e-6f);
					CHECK(set == 7);
					CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTriangle<true>(a, b, c, set) / scale, expected / scale, 2.0e-6f);
					CHECK(set == 7);
				} while (std::next_permutation(permutation, permutation + 3));
			}
	}

	TEST_CASE("TestLowerDimensionalDegenerateAndConditionedProjections")
	{
		for (float scale : { 0x1p-30f, 1.0f, 0x1p30f })
		{
			CAPTURE(scale);
			const Vec4 point = Vec4(1, -2, 3, -4) * scale;
			float u, v, w;
			CHECK_FALSE(ClosestPoint::GetBaryCentricCoordinates(point, point, u, v));
			CHECK(u + v == 1.0f);
			CHECK_FALSE(ClosestPoint::GetBaryCentricCoordinates(point, point, point, u, v, w));
			CHECK(u + v + w == 1.0f);
			uint32 set;
			CHECK(ClosestPoint::GetClosestPointOnTriangle(point, point, point, set) == point);
			CHECK(set == 4);
			CHECK(ClosestPoint::GetClosestPointOnTriangle<true>(point, point, point, set) == point);
			CHECK(set == 4);
			// Endpoints and an interior projection, checked under endpoint exchange.
			for (float query : { -2.0f, 0.0f, 2.0f })
				for (bool reverse : { false, true })
				{
					Vec4 a = Vec4(-1 - query, 0.125f, 0, 0) * scale;
					Vec4 b = Vec4(1 - query, 0.125f, 0, 0) * scale;
					if (reverse) std::swap(a, b);
					Vec4 expected(Clamp(query, -1.0f, 1.0f) - query, 0.125f, 0, 0);
					CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnLine(a, b, set) / scale, expected, 1.0e-6f);
					CHECK(set == (query == 0? 3u : ((query < 0) != reverse? 1u : 2u)));
				}
			for (float height : { 0.0f, 0x1p-12f, 0x1p-8f })
			{
				CAPTURE(height);
				const Vec4 basis[] = { Vec4(-1, -height * 0.25f, 0.125f, 0.25f),
					Vec4(1, -height * 0.25f, 0.125f, 0.25f), Vec4(0, height * 0.75f, 0.125f, 0.25f) };
				int permutation[] = { 0, 1, 2 };
				do
				{
					Vec4 a = basis[permutation[0]] * scale, b = basis[permutation[1]] * scale, c = basis[permutation[2]] * scale;
					bool conditioned = height == 0x1p-8f;
					CHECK(ClosestPoint::GetBaryCentricCoordinates(a, b, c, u, v, w) == conditioned);
					Vec4 expected(0, conditioned? 0 : -height * 0.25f, 0.125f, 0.25f);
					CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTriangle(a, b, c, set) / scale, expected, 1.0e-6f);
					uint32 base_set = (1u << int(std::find(permutation, permutation + 3, 0) - permutation))
						| (1u << int(std::find(permutation, permutation + 3, 1) - permutation));
					if (height > 0) CHECK(set == (conditioned? 7u : base_set));
					// Only invoke the optimized path when its stated C-in-feature precondition holds.
					if ((set & 4u) != 0)
						CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTriangle<true>(a, b, c, set) / scale, expected, 1.0e-6f);
				} while (std::next_permutation(permutation, permutation + 3));
			}
		}
	}

	TEST_CASE("TestTriangleSmallScaleKnownSupportProjection")
	{
		// Last pentachoron of envelope sample 0 / scale .01 / margin .001 / W / sphere first.
		// Exact-rational active-set solve chooses B,C,E: its area^2 is 5.0801358435e-11
		// but area^2/(|BC|^2 |BE|^2) is .849985, so this is NOT a sliver.
		const Vec4 fixture[] = {
			UVec4(3124139324u, 3124762040u, 3125614128u, 3113572352u).ReinterpretAsFloat(),
			UVec4(980414140u, 970477680u, 975393232u, 941508608u).ReinterpretAsFloat(),
			UVec4(3137814712u, 977921368u, 986919216u, 948463104u).ReinterpretAsFloat()
		};
		const float weights[] = { 0.4413719031941673f, 0.494423786144555f, 0.0642043106612777f };
		const Vec4 reference(0.000012611873664899887f, -0.00006997230677438714f,
			0.00006617920941013659f, -0.00010289104727044411f);
		for (float scale : { 0.125f, 1.0f, 8.0f, 128.0f })
		{
			CAPTURE(scale);
			int permutation[] = { 0, 1, 2 };
			do
			{
				Vec4 a = fixture[permutation[0]] * scale, b = fixture[permutation[1]] * scale, c = fixture[permutation[2]] * scale;
				float u, v, w;
				CHECK(ClosestPoint::GetBaryCentricCoordinates(a, b, c, u, v, w));
				CHECK_APPROX_EQUAL(u, weights[permutation[0]], 1.0e-6f);
				CHECK_APPROX_EQUAL(v, weights[permutation[1]], 1.0e-6f);
				CHECK_APPROX_EQUAL(w, weights[permutation[2]], 1.0e-6f);
				uint32 set;
				CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTriangle(a, b, c, set), reference * scale, 1.0e-9f * scale);
				CHECK(set == 7);
				CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTriangle<true>(a, b, c, set), reference * scale, 1.0e-9f * scale);
				CHECK(set == 7);
			} while (std::next_permutation(permutation, permutation + 3));
		}
	}

	TEST_CASE("TestTetrahedronSmallScaleProjection")
	{
		// The origin projects to the centroid with barycentrics (1/4, 1/4, 1/4, 1/4).
		// A fixed 1e-12 cutoff on a length^6 determinant incorrectly discarded this interior.
		const Vec4 basis[] = { Vec4(-1, -1, -1, 0.125f), Vec4(1, 0, 0, 0.125f),
			Vec4(0, 1, 0, 0.125f), Vec4(0, 0, 1, 0.125f) };
		for (float scale : { 0.001f, 0.01f, 0.1f, 1.0f, 100.0f })
		{
			CAPTURE(scale);
			int permutation[] = { 0, 1, 2, 3 };
			do
			{
				Vec4 a = basis[permutation[0]] * scale, b = basis[permutation[1]] * scale;
				Vec4 c = basis[permutation[2]] * scale, d = basis[permutation[3]] * scale;
				CHECK(ClosestPoint::OriginOutsideOfTetrahedronPlanes(a, b, c, d).GetTrues() == 0);
				float u, v, w, x;
				CHECK(ClosestPoint::GetBaryCentricCoordinates(a, b, c, d, u, v, w, x));
				for (float weight : { u, v, w, x })
					CHECK_APPROX_EQUAL(weight, 0.25f, 1.0e-5f);
				uint32 set;
				Vec4 expected(0, 0, 0, 0.125f * scale);
				CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTetrahedron(a, b, c, d, set), expected, 1.0e-5f * scale);
				CHECK(set == 15);
				CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTetrahedron<true>(a, b, c, d, set), expected, 1.0e-5f * scale);
				CHECK(set == 15);
			} while (std::next_permutation(permutation, permutation + 4));
		}
	}

	TEST_CASE("TestTetrahedronSmallScaleKnownSupportProjection")
	{
		// Exact support simplex at the no-progress exit of seeded GJK sample 0, scale 0.01,
		// initial axis W. Independent 65-digit solve (not float Gram cofactors)
		// gives barycentrics [0.00832342377, 0.37596215505, 0.49389370775, 0.12182071343].
		const Vec4 fixture[] = {
			UVec4(1010285740u, 3163555020u, 1015486452u, 996166808u).ReinterpretAsFloat(),
			UVec4(3147128350u, 989130896u, 961850688u, 3127518464u).ReinterpretAsFloat(),
			UVec4(996926212u, 3123856728u, 3127008288u, 3098048256u).ReinterpretAsFloat(),
			UVec4(3133775915u, 3140428114u, 983431248u, 989491424u).ReinterpretAsFloat()
		};
		const Vec4 reference(-0.00003012045801648818f, -0.00010516357624502538f,
			-0.00006613909852123448f, -0.00011422888476000912f);
		for (float scale : { 0.1f, 1.0f, 10.0f, 100.0f })
		{
			CAPTURE(scale);
			int permutation[] = { 0, 1, 2, 3 };
			do
			{
				Vec4 a = fixture[permutation[0]] * scale, b = fixture[permutation[1]] * scale;
				Vec4 c = fixture[permutation[2]] * scale, d = fixture[permutation[3]] * scale;
				CHECK(ClosestPoint::OriginOutsideOfTetrahedronPlanes(a, b, c, d).GetTrues() == 0);
				float u, v, w, x;
				CHECK(ClosestPoint::GetBaryCentricCoordinates(a, b, c, d, u, v, w, x));
				uint32 set;
				CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTetrahedron(a, b, c, d, set), reference * scale, 1.0e-7f * scale);
				CHECK(set == 15);
				CHECK_APPROX_EQUAL(ClosestPoint::GetClosestPointOnTetrahedron<true>(a, b, c, d, set), reference * scale, 1.0e-7f * scale);
				CHECK(set == 15);
			} while (std::next_permutation(permutation, permutation + 4));
		}
	}

	TEST_CASE("TestPentachoronFacesAgainstKnownBarycentrics")
	{
		UnitTestRandom random(98765);
		std::uniform_real_distribution<float> angle(-1.0f, 1.0f);
		std::uniform_int_distribution<int> weight(-8, 8);
		for (int sample = 0; sample < 1000; ++sample)
		{
			float bary[5];
			bary[0] = 1.0f;
			for (int i = 1; i < 5; ++i)
			{
				int w = weight(random);
				bary[i] = float(w == 0? 1 : w) * 0.125f;
				bary[0] -= bary[i];
			}
			if (bary[0] == 0.0f) continue; // Keep away from the rounding-sensitive boundary.
			Rotor rotation = Rotor::sExp(Bivec(angle(random), angle(random), angle(random), angle(random), angle(random), angle(random)));
			for (float scale : { 0.001f, 0.1f, 1.0f, 10.0f, 1000.0f })
			{
				Vec4 query(bary[1], bary[2], bary[3], bary[4]);
				Vec4 vertices[5] = { Vec4::sZero(), Vec4::sAxisX(), Vec4::sAxisY(), Vec4::sAxisZ(), Vec4::sAxisW() };
				for (Vec4 &vertex : vertices) vertex = rotation * ((vertex - query) * scale);
				uint32 expected = (bary[4] < 0.0f? 1u : 0u) | (bary[0] < 0.0f? 2u : 0u)
					| (bary[1] < 0.0f? 4u : 0u) | (bary[2] < 0.0f? 8u : 0u) | (bary[3] < 0.0f? 16u : 0u);
				CHECK(ClosestPoint::OriginOutsideOfPentachoronPlanes(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4]) == expected);
			}
		}
		CHECK(ClosestPoint::OriginOutsideOfPentachoronPlanes(Vec4::sZero(), Vec4::sAxisX(), Vec4::sAxisY(), Vec4::sAxisZ(), Vec4(0, 0, 1, 1.0e-5f)) == 0x1f);
	}

	TEST_CASE("TestPentachoronKnownInteriorCancellation")
	{
		// Iteration 6 of TestGJKSphereVsBoxKnownFalseNegative, preserved without decimal rounding.
		// Independent high-precision solve gives barycentrics
		// [0.29662915, 0.30951689, 0.28358813, 0.09363796, 0.01662787]: all strictly inside.
		const Vec4 fixture[] = {
			UVec4(1025978056u, 1018103232u, 1075501340u, 1018264664u).ReinterpretAsFloat(),
			UVec4(1021253184u, 3175558144u, 3215635887u, 1013765576u).ReinterpretAsFloat(),
			UVec4(3197819333u, 1016904320u, 3208892514u, 3190056164u).ReinterpretAsFloat(),
			UVec4(1057573245u, 1036647744u, 3208960858u, 1056745198u).ReinterpretAsFloat(),
			UVec4(1063728008u, 3199673500u, 3212183460u, 3207743627u).ReinterpretAsFloat()
		};
		for (float scale : { 0.01f, 1.0f, 100.0f })
		{
			CAPTURE(scale);
			int permutation[] = { 0, 1, 2, 3, 4 };
			int count = 0;
			do
			{
				CAPTURE(count);
				Vec4 a = fixture[permutation[0]] * scale, b = fixture[permutation[1]] * scale;
				Vec4 c = fixture[permutation[2]] * scale, d = fixture[permutation[3]] * scale, e = fixture[permutation[4]] * scale;
				CHECK(ClosestPoint::OriginOutsideOfPentachoronPlanes(a, b, c, d, e) == 0);
				uint32 set = 0;
				Vec4 closest = ClosestPoint::GetClosestPointOnPentachoron(a, b, c, d, e, set);
				CHECK(set == 0x1f);
				CHECK(closest == Vec4::sZero());
				closest = ClosestPoint::GetClosestPointOnPentachoron<true>(a, b, c, d, e, set);
				CHECK(set == 0x1f);
				CHECK(closest == Vec4::sZero());
				++count;
			} while (std::next_permutation(permutation, permutation + 5));
			CHECK(count == 120);
		}
	}

	TEST_CASE("TestPentachoronSliversAgainstKnownBarycentrics")
	{
		float boundary_weight = 0.03125f;
		SUBCASE("Interior and exterior") { boundary_weight = 0.03125f; }
		SUBCASE("Near face interior and exterior") { boundary_weight = 1.0f / 32768.0f; }
		UnitTestRandom random(0x4d51a7u);
		for (int sample = 0; sample < 48; ++sample)
		{
			CAPTURE(sample);
			// Nearly parallel edges: the Gram determinant/Hadamard ratio is about 3.6e-6,
			// just above the existing 1e-6 degeneracy threshold, not a rank-deficient simplex.
			const Vec4 basis[] = { Vec4::sZero(), Vec4(1, 0, 0, 0), Vec4(1, 0.125f, 0, 0),
				Vec4(1, 0, 0.125f, 0), Vec4(1, 0, 0, 0.125f) };
			float bary[5] = { 1, 0, 0, 0, 0 };
			for (int i = 1; i < 5; ++i)
			{
				bary[i] = float(1 + (random() & 7u)) / 64.0f;
				bary[0] -= bary[i];
			}
			// Exercise every outside-face bit as well as strict interior containment.
			int outside = sample % 6;
			int near_face = outside < 5? outside : 4;
			float near_weight = outside < 5? -boundary_weight : boundary_weight;
			bary[(near_face + 1) % 5] += bary[near_face] - near_weight;
			bary[near_face] = near_weight;
			Vec4 query = Vec4::sZero();
			for (int i = 0; i < 5; ++i)
				query += bary[i] * basis[i];
			float angles[6];
			for (float &a : angles)
				a = float(int(random() & 1023u) - 512) / 512.0f;
			Rotor rotation = Rotor::sExp(Bivec(angles[0], angles[1], angles[2], angles[3], angles[4], angles[5]));
			for (float scale : { 0.001f, 0.01f, 1.0f, 100.0f, 1000.0f })
			{
				CAPTURE(scale);
				Vec4 vertices[5];
				for (int i = 0; i < 5; ++i)
					vertices[i] = rotation * ((basis[i] - query) * scale);
				int permutation[] = { 0, 1, 2, 3, 4 };
				int count = 0;
				do
				{
					CAPTURE(count);
					Vec4 a = vertices[permutation[0]], b = vertices[permutation[1]], c = vertices[permutation[2]];
					Vec4 d = vertices[permutation[3]], e = vertices[permutation[4]];
					// The oracle is the prescribed barycentrics, not another Gram/cofactor solve.
					uint32 expected = (bary[permutation[4]] < 0? 1u : 0u) | (bary[permutation[0]] < 0? 2u : 0u)
						| (bary[permutation[1]] < 0? 4u : 0u) | (bary[permutation[2]] < 0? 8u : 0u)
						| (bary[permutation[3]] < 0? 16u : 0u);
					CHECK(ClosestPoint::OriginOutsideOfPentachoronPlanes(a, b, c, d, e) == expected);
					if (expected == 0)
					{
						uint32 set = 0;
						CHECK(ClosestPoint::GetClosestPointOnPentachoron<true>(a, b, c, d, e, set) == Vec4::sZero());
						CHECK(set == 0x1f);
					}
					++count;
				} while (std::next_permutation(permutation, permutation + 5));
				CHECK(count == 120);
			}
		}
	}

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
