// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Math/Mat44.h>

JPH_SUPPRESS_WARNINGS_STD_BEGIN
#include <random>
JPH_SUPPRESS_WARNINGS_STD_END

// TODO(4D): Mat44::GetDeterminant is new in the 4D fork (upstream only ever needed the 3x3
// determinant of the rotation part). It gates the invertibility of the 4x4 effective mass in
// PointConstraintPart, so it gets a reference test like the rest of the new math.
TEST_SUITE("Mat44DeterminantTests")
{
	/// Brute force determinant from the Leibniz formula: sum over all 24 permutations of
	/// sign(sigma) * prod_i M(i, sigma(i)). Deliberately as close to the definition as possible.
	static double sReferenceDeterminant(Mat44Arg inM)
	{
		int perm[4] = { 0, 1, 2, 3 };
		double det = 0.0;

		// Iterate over all 24 permutations in lexicographic order
		do
		{
			// Sign of the permutation = (-1)^(number of inversions)
			int inversions = 0;
			for (int i = 0; i < 4; ++i)
				for (int j = i + 1; j < 4; ++j)
					if (perm[i] > perm[j])
						++inversions;

			double term = (inversions & 1) != 0? -1.0 : 1.0;
			for (int i = 0; i < 4; ++i)
				term *= double(inM(uint(i), uint(perm[i])));
			det += term;
		} while (std::next_permutation(perm, perm + 4));

		return det;
	}

	static Mat44 sRandomMat(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> dist(-2.0f, 2.0f);
		Mat44 m;
		for (uint c = 0; c < 4; ++c)
			m.SetColumn4(c, Lane4(dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom)));
		return m;
	}

	TEST_CASE("TestMat44DeterminantVsReference")
	{
		UnitTestRandom random(0xDE7);
		for (int n = 0; n < 100; ++n)
		{
			Mat44 m = sRandomMat(random);
			CHECK(abs(double(m.GetDeterminant()) - sReferenceDeterminant(m)) < 1.0e-3);
		}
	}

	TEST_CASE("TestMat44DeterminantIdentityAndScale")
	{
		CHECK_APPROX_EQUAL(Mat44::sIdentity().GetDeterminant(), 1.0f, 1.0e-6f);
		CHECK_APPROX_EQUAL(Mat44::sZero().GetDeterminant(), 0.0f, 1.0e-6f);

		// A diagonal matrix has the product of its diagonal as determinant
		Mat44 d;
		d.SetColumn4(0, Lane4(2, 0, 0, 0));
		d.SetColumn4(1, Lane4(0, 3, 0, 0));
		d.SetColumn4(2, Lane4(0, 0, 5, 0));
		d.SetColumn4(3, Lane4(0, 0, 0, 7));
		CHECK_APPROX_EQUAL(d.GetDeterminant(), 210.0f, 1.0e-4f);
	}

	TEST_CASE("TestMat44DeterminantOfRotationIsOne")
	{
		// Rotations preserve 4-volume and orientation
		Mat44 r = Mat44::sRotationXY(0.3f) * Mat44::sRotationZW(-1.1f) * Mat44::sRotationYW(2.4f);
		CHECK_APPROX_EQUAL(r.GetDeterminant(), 1.0f, 1.0e-5f);
	}

	TEST_CASE("TestMat44DeterminantIsMultiplicative")
	{
		UnitTestRandom random(0xDE8);
		for (int n = 0; n < 25; ++n)
		{
			Mat44 a = sRandomMat(random), b = sRandomMat(random);
			float lhs = (a * b).GetDeterminant();
			float rhs = a.GetDeterminant() * b.GetDeterminant();
			CHECK(abs(lhs - rhs) < 1.0e-3f * max(1.0f, abs(rhs)));
		}
	}

	TEST_CASE("TestMat44DeterminantDetectsSingular")
	{
		UnitTestRandom random(0xDE9);
		for (int n = 0; n < 25; ++n)
		{
			// Make column 3 a linear combination of the others: the matrix is now rank deficient
			Mat44 m = sRandomMat(random);
			m.SetColumn4(3, m.GetColumn4(0) * 2.0f - m.GetColumn4(1));
			CHECK(abs(m.GetDeterminant()) < 1.0e-3f);
		}
	}

	TEST_CASE("TestMat44NonZeroDeterminantImpliesInvertible")
	{
		UnitTestRandom random(0xDEA);
		for (int n = 0; n < 25; ++n)
		{
			Mat44 m = sRandomMat(random);
			if (abs(m.GetDeterminant()) > 0.1f)
				CHECK((m * m.Inversed()).IsClose(Mat44::sIdentity(), 1.0e-4f));
		}
	}
}
