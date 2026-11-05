// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Core/FPException.h>
#include <cmath>

JPH_NAMESPACE_BEGIN

/// Function to determine the eigen vectors and values of a N x N real symmetric matrix
/// by Jacobi transformations. This method is most suitable for N < 10.
///
/// Taken and adapted from Numerical Recipes paragraph 11.1
///
/// An eigen vector is a vector v for which \f$A \: v = \lambda \: v\f$
///
/// Where:
/// A: A square matrix.
/// \f$\lambda\f$: a non-zero constant value.
///
/// @see https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors
///
/// Matrix is a matrix type, which has dimensions N x N.
/// @param inMatrix is the matrix of which to return the eigenvalues and vectors
/// @param outEigVec will contain a matrix whose columns contain the normalized eigenvectors (must be identity before call)
/// @param outEigVal will contain the eigenvalues
/// @param inCheckAccuracy Convert numerical sanity-check / convergence failures
/// to false instead of assertions, in all builds. Unsupported intermediate float
/// ranges are rejected before arithmetic that could overflow; output is
/// unspecified on failure. No additional exception masking or relaxed tolerances.
/// The default preserves the historical assertion-based diagnostic behavior.
template <class Vector, class Matrix>
bool EigenValueSymmetric(const Matrix &inMatrix, Matrix &outEigVec, Vector &outEigVal, bool inCheckAccuracy = false)
{
	// This algorithm can generate infinite values, see comment below
	FPExceptionDisableInvalid disable_invalid;
	JPH_UNUSED(disable_invalid);

	// Maximum number of sweeps to make
	const int cMaxSweeps = 50;

	// Get problem dimension
	const uint n = inMatrix.GetRows();

	// Make sure the dimensions are right
	JPH_ASSERT(inMatrix.GetRows() == n);
	JPH_ASSERT(inMatrix.GetCols() == n);
	JPH_ASSERT(outEigVec.GetRows() == n);
	JPH_ASSERT(outEigVec.GetCols() == n);
	JPH_ASSERT(outEigVal.GetRows() == n);
	JPH_ASSERT(outEigVec.IsIdentity());

	// Checked callers must survive overflow traps, including in the historical
	// float squared-norm diagnostics. Reject unsupported ranges BEFORE float
	// arithmetic; masking overflow would let an infinity tolerance pass.
	// The n^2 margin also bounds the Jacobi sums/products and diagnostic
	// residuals. This is a range restriction, NOT an accuracy/rank tolerance.
	if (inCheckAccuracy)
	{
		double max_entry = std::sqrt(double(FLT_MAX)) / (16.0 * n * n);
		for (uint r = 0; r < n; ++r)
			for (uint c = 0; c < n; ++c)
				if (!std::isfinite(inMatrix(r, c)) || std::abs(double(inMatrix(r, c))) > max_entry)
					return false;
	}

	// Get the matrix in a so we can mess with it
	Matrix a = inMatrix;

	Vector b, z;

	for (uint ip = 0; ip < n; ++ip)
	{
		// Initialize b to diagonal of a
		b[ip] = a(ip, ip);

		// Initialize output to diagonal of a
		outEigVal[ip] = a(ip, ip);

		// Reset z
		z[ip] = 0.0f;
	}

	for (int sweep = 0; sweep < cMaxSweeps; ++sweep)
	{
		// Get the sum of the off-diagonal elements of a
		float sm = 0.0f;
		for (uint ip = 0; ip < n - 1; ++ip)
			for (uint iq = ip + 1; iq < n; ++iq)
				sm += abs(a(ip, iq));
		float avg_sm = sm / Square(n);

		// Normal return, convergence to machine underflow
		if (avg_sm < FLT_MIN) // Original code: sm == 0.0f, when the average is denormal, we also consider it machine underflow
		{
			// Checked callers use exactly the historical diagnostics and
			// tolerances, including in release builds, but get recoverable failure.
			bool check_accuracy = inCheckAccuracy;
			JPH_IF_ENABLE_ASSERTS(check_accuracy = true;)
			if (check_accuracy)
				for (uint c = 0; c < n; ++c)
				{
					bool normalized = outEigVec.GetColumn(c).IsNormalized();
					if (inCheckAccuracy && !normalized)
						return false;
					JPH_ASSERT(normalized);

					Vector mat_eigvec = inMatrix * outEigVec.GetColumn(c);
					Vector eigval_eigvec = outEigVal[c] * outEigVec.GetColumn(c);
					float mat_length_sq = mat_eigvec.LengthSq(), eigval_length_sq = eigval_eigvec.LengthSq();
					if (inCheckAccuracy && (!std::isfinite(mat_length_sq) || !std::isfinite(eigval_length_sq)))
						return false;
					bool accurate = mat_eigvec.IsClose(eigval_eigvec, max(mat_length_sq, eigval_length_sq) * 1.0e-6f);
					if (inCheckAccuracy && !accurate)
						return false;
					JPH_ASSERT(accurate);
				}

			// Success
			return true;
		}

		// On the first three sweeps use a fraction of the sum of the off diagonal elements as threshold
		// Note that we pick a minimum threshold of FLT_MIN because dividing by a denormalized number is likely to result in infinity.
		float thresh = sweep < 4? 0.2f * avg_sm : FLT_MIN; // Original code: 0.0f instead of FLT_MIN

		for (uint ip = 0; ip < n - 1; ++ip)
			for (uint iq = ip + 1; iq < n; ++iq)
			{
				float &a_pq = a(ip, iq);
				float &eigval_p = outEigVal[ip];
				float &eigval_q = outEigVal[iq];

				float abs_a_pq = abs(a_pq);
				float g = 100.0f * abs_a_pq;

				// After four sweeps, skip the rotation if the off-diagonal element is small
				if (sweep > 4
					&& abs(eigval_p) + g == abs(eigval_p)
					&& abs(eigval_q) + g == abs(eigval_q))
				{
					a_pq = 0.0f;
				}
				else if (abs_a_pq > thresh)
				{
					float h = eigval_q - eigval_p;
					float abs_h = abs(h);

					float t;
					if (abs_h + g == abs_h)
					{
						t = a_pq / h;
					}
					else
					{
						// A finite input can still produce a huge theta. Checked
						// success must imply that the following ordinary float
						// division AND theta square are safe on the decoder's rerun.
						if (inCheckAccuracy && std::abs(0.5 * double(h) / a_pq) > 0.5 * std::sqrt(double(FLT_MAX)))
							return false;
						float theta = 0.5f * h / a_pq; // Warning: Can become infinite if a(ip, iq) is very small which may trigger an invalid float exception
						t = 1.0f / (abs(theta) + sqrt(1.0f + theta * theta)); // If theta becomes inf, t will be 0 so the infinite is not a problem for the algorithm
						if (theta < 0.0f) t = -t;
					}

					float c = 1.0f / sqrt(1.0f + t * t);
					float s = t * c;
					float tau = s / (1.0f + c);
					h = t * a_pq;

					a_pq = 0.0f;

					z[ip] -= h;
					z[iq] += h;

					eigval_p -= h;
					eigval_q += h;

					#define JPH_EVS_ROTATE(a, i, j, k, l)		\
						g = a(i, j),							\
						h = a(k, l),							\
						a(i, j) = g - s * (h + g * tau),		\
						a(k, l) = h + s * (g - h * tau)

					uint j;
					for (j = 0; j < ip; ++j)		JPH_EVS_ROTATE(a, j, ip, j, iq);
					for (j = ip + 1; j < iq; ++j)	JPH_EVS_ROTATE(a, ip, j, j, iq);
					for (j = iq + 1; j < n; ++j)	JPH_EVS_ROTATE(a, ip, j, iq, j);
					for (j = 0; j < n; ++j)			JPH_EVS_ROTATE(outEigVec, j, ip, j, iq);

					#undef JPH_EVS_ROTATE
				}
			}

		// Update eigenvalues with the sum of ta_pq and reinitialize z
		for (uint ip = 0; ip < n; ++ip)
		{
			b[ip] += z[ip];
			outEigVal[ip] = b[ip];
			z[ip] = 0.0f;
		}
	}

	// Failure. Checked callers must not enter the assertion handler.
	if (!inCheckAccuracy)
		JPH_ASSERT(false, "Too many iterations");
	return false;
}

JPH_NAMESPACE_END
