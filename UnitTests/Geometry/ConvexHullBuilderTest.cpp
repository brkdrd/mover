// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Geometry/ConvexHullBuilder.h>

TEST_SUITE("ConvexHullBuilderTest")
{
	constexpr float cTolerance = 1.0e-3f;
	using Cell = ConvexHullBuilder::Cell;
	using Positions = ConvexHullBuilder::Positions;

	TEST_CASE("TestDegenerate4D")
	{
		const char *error = nullptr;

		{
			// Too few points (< 5)
			Positions positions { Vec4(1, 2, 3, 4) };

			ConvexHullBuilder builder(positions);
			CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::TooFewPoints);
			positions.push_back(Vec4(2, 2, 3, 4));
			CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::TooFewPoints);
			positions.push_back(Vec4(1, 3, 3, 4));
			CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::TooFewPoints);
			positions.push_back(Vec4(1, 2, 4, 4));
			CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::TooFewPoints);
		}

		{
			// 5 collinear points → degenerate
			Positions positions;
			for (float v = 0.0f; v < 1.01f; v += 0.25f)
				positions.push_back(Vec4(v, 0, 0, 0));

			ConvexHullBuilder builder(positions);
			CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::Degenerate);
		}

		{
			// 5 coplanar points (in a 2D plane) → degenerate
			Positions positions;
			positions.push_back(Vec4(0, 0, 0, 0));
			positions.push_back(Vec4(1, 0, 0, 0));
			positions.push_back(Vec4(0, 1, 0, 0));
			positions.push_back(Vec4(1, 1, 0, 0));
			positions.push_back(Vec4(0.5f, 0.5f, 0, 0));

			ConvexHullBuilder builder(positions);
			CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::Degenerate);
		}

		{
			// 5 co-hyperplanar points (in 3D subspace) → degenerate
			Positions positions;
			positions.push_back(Vec4(0, 0, 0, 0));
			positions.push_back(Vec4(1, 0, 0, 0));
			positions.push_back(Vec4(0, 1, 0, 0));
			positions.push_back(Vec4(0, 0, 1, 0));
			positions.push_back(Vec4(0.5f, 0.5f, 0.5f, 0));

			ConvexHullBuilder builder(positions);
			CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::Degenerate);
		}
	}

	TEST_CASE("TestPentachoron")
	{
		const char *error = nullptr;

		// Standard pentachoron (4-simplex): 5 vertices, 5 tetrahedral cells
		Positions positions;
		positions.push_back(Vec4(1, 0, 0, 0));
		positions.push_back(Vec4(0, 1, 0, 0));
		positions.push_back(Vec4(0, 0, 1, 0));
		positions.push_back(Vec4(0, 0, 0, 1));
		positions.push_back(Vec4(0, 0, 0, 0));

		ConvexHullBuilder builder(positions);
		CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::Success);
		CHECK(builder.GetNumVerticesUsed() == 5);
		CHECK(builder.GetCells().size() == 5);

		// Verify each cell is a valid tetrahedron
		for (const Cell *c : builder.GetCells())
		{
			CHECK_FALSE(c->mRemoved);
			// All 4 neighbours must be set
			for (int j = 0; j < 4; ++j)
				CHECK(c->mNeighbour[j] != nullptr);
		}

		// Calculate error
		Cell *error_cell;
		float max_error, coplanar_distance;
		int error_position_idx;
		builder.DetermineMaxError(error_cell, max_error, error_position_idx, coplanar_distance);
		CHECK(max_error < max(coplanar_distance, cTolerance));
	}

	TEST_CASE("TestHypercube")
	{
		const char *error = nullptr;

		// 4D hypercube (tesseract): 16 vertices
		Positions positions;
		for (int x = 0; x < 2; ++x)
			for (int y = 0; y < 2; ++y)
				for (int z = 0; z < 2; ++z)
					for (int w = 0; w < 2; ++w)
						positions.push_back(Vec4((float)x, (float)y, (float)z, (float)w));

		ConvexHullBuilder builder(positions);
		CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::Success);
		CHECK(builder.GetNumVerticesUsed() == 16);

		// A tesseract has 24 octahedral cells in its boundary... but since we use simplicial cells
		// (tetrahedra), each octahedral cell is triangulated. The exact count depends on the
		// triangulation. Just verify we got a reasonable number and low error.
		CHECK(builder.GetCells().size() >= 24);

		Cell *error_cell;
		float max_error, coplanar_distance;
		int error_position_idx;
		builder.DetermineMaxError(error_cell, max_error, error_position_idx, coplanar_distance);
		CHECK(max_error < max(coplanar_distance, cTolerance));
	}

	TEST_CASE("TestHypercubeWithInterior")
	{
		const char *error = nullptr;

		// 4D cube with interior points (3^4 = 81 points)
		Positions positions;
		for (int x = 0; x < 3; ++x)
			for (int y = 0; y < 3; ++y)
				for (int z = 0; z < 3; ++z)
					for (int w = 0; w < 3; ++w)
						positions.push_back(Vec4(0.1f * x, 0.2f * y, 0.3f * z, 0.4f * w));

		ConvexHullBuilder builder(positions);
		CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::Success);
		CHECK(builder.GetNumVerticesUsed() >= 16); // At least all 16 corner vertices must be used

		Cell *error_cell;
		float max_error, coplanar_distance;
		int error_position_idx;
		builder.DetermineMaxError(error_cell, max_error, error_position_idx, coplanar_distance);
		CHECK(max_error < max(coplanar_distance, cTolerance));
	}

	TEST_CASE("TestCenterOfMassAndVolume")
	{
		const char *error = nullptr;

		// Unit pentachoron: vertices at origin and along each axis
		Positions positions;
		positions.push_back(Vec4(0, 0, 0, 0));
		positions.push_back(Vec4(1, 0, 0, 0));
		positions.push_back(Vec4(0, 1, 0, 0));
		positions.push_back(Vec4(0, 0, 1, 0));
		positions.push_back(Vec4(0, 0, 0, 1));

		ConvexHullBuilder builder(positions);
		CHECK(builder.Initialize(INT_MAX, cTolerance, error) == ConvexHullBuilder::EResult::Success);

		Vec4 com;
		float volume;
		builder.GetCenterOfMassAndVolume(com, volume);

		// Center of mass of a regular simplex = average of vertices
		CHECK(abs(com.GetX() - 0.2f) < 0.01f);
		CHECK(abs(com.GetY() - 0.2f) < 0.01f);
		CHECK(abs(com.GetZ() - 0.2f) < 0.01f);
		CHECK(abs(com.GetW() - 0.2f) < 0.01f);

		// Volume of a unit 4-simplex = 1/24
		CHECK(abs(abs(volume) - 1.0f / 24.0f) < 0.001f);
	}

	TEST_CASE("TestRandom4DHull")
	{
		const char *error = nullptr;

		UnitTestRandom random(0x4d4d4d4d);

		uniform_real_distribution<float> scale_start(0.1f, 0.5f);
		uniform_real_distribution<float> scale_range(0.1f, 2.0f);
		uniform_real_distribution<float> zero_one(0.0f, 1.0f);
		uniform_real_distribution<float> zero_two(0.0f, 2.0f);

		for (int iteration = 0; iteration < 50; ++iteration)
		{
			float start = scale_start(random);
			uniform_real_distribution<float> vertex_scale(start, start + scale_range(random));

			// Non-uniform scaling to avoid near-spherical shapes
			uniform_real_distribution<float> shape_scale(0.1f, 1.0f);
			float sx = shape_scale(random), sy = shape_scale(random);
			float sz = shape_scale(random), sw = shape_scale(random);

			Positions positions;
			for (int i = 0; i < 60; ++i)
			{
				// Random 4D point on the 3-sphere, scaled
				float r = vertex_scale(random);
				float phi1 = 2.0f * JPH_PI * zero_one(random);
				float phi2 = 2.0f * JPH_PI * zero_one(random);
				float theta = JPH_PI * zero_one(random);
				float st = Sin(theta), ct = Cos(theta);
				float sp1 = Sin(phi1), cp1 = Cos(phi1);
				float sp2 = Sin(phi2), cp2 = Cos(phi2);
				Vec4 p1 = Vec4(r * st * cp1 * cp2 * sx,
							   r * st * cp1 * sp2 * sy,
							   r * st * sp1 * sz,
							   r * ct * sw);
				positions.push_back(p1);

				// Point close to p1
				Vec4 offset(cTolerance * zero_two(random) * (zero_one(random) - 0.5f),
							cTolerance * zero_two(random) * (zero_one(random) - 0.5f),
							cTolerance * zero_two(random) * (zero_one(random) - 0.5f),
							cTolerance * zero_two(random) * (zero_one(random) - 0.5f));
				positions.push_back(p1 + offset);
			}

			ConvexHullBuilder builder(positions);
			ConvexHullBuilder::EResult result = builder.Initialize(INT_MAX, cTolerance, error);
			CHECK(result == ConvexHullBuilder::EResult::Success);

			if (result == ConvexHullBuilder::EResult::Success)
			{
				Cell *error_cell;
				float max_error, coplanar_distance;
				int error_position_idx;
				builder.DetermineMaxError(error_cell, max_error, error_position_idx, coplanar_distance);
				CHECK(max_error < max(coplanar_distance, 1.2f * cTolerance));
			}
		}
	}

	TEST_CASE("TestMaxVertices")
	{
		const char *error = nullptr;

		// Generate many 4D points, limit vertices
		Positions positions;
		UnitTestRandom random(42);
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		for (int i = 0; i < 200; ++i)
			positions.push_back(Vec4(dist(random), dist(random), dist(random), dist(random)));

		ConvexHullBuilder builder(positions);
		ConvexHullBuilder::EResult result = builder.Initialize(10, cTolerance, error);
		CHECK((result == ConvexHullBuilder::EResult::MaxVerticesReached || result == ConvexHullBuilder::EResult::Success));
	}
}
