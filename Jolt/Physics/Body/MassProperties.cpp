// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Body/MassProperties.h>
#include <Jolt/Math/Matrix.h>
#include <Jolt/Math/Vector.h>
#include <Jolt/Math/EigenValueSymmetric.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>
#include <Jolt/Core/StreamIn.h>
#include <Jolt/Core/StreamOut.h>
#include <Jolt/Core/InsertionSort.h>
#include <cmath>

JPH_NAMESPACE_BEGIN

JPH_IMPLEMENT_SERIALIZABLE_NON_VIRTUAL(MassProperties)
{
	JPH_ADD_ATTRIBUTE(MassProperties, mMass)
	JPH_ADD_ATTRIBUTE(MassProperties, mInertia)
}

bool MassProperties::DecomposePrincipalMomentsOfInertia(Mat44 &outRotation, Vec4 &outDiagonal) const
{
	Mat44 rotation = Mat44::sZero();
	Vec4 diagonal = Vec4::sZero();

	// Certify exact rank one BEFORE the float eigensolver. Its relative
	// eigenpair residual test is undefined on the null space, and roundoff
	// there must not become a positive moment with an enormous reciprocal.
	// Products of two binary32 values are exact in binary64 (48 <= 53 bits):
	// these minor tests have NO tolerance and cannot erase a skinny positive
	// eigenvalue. This recognizes rank of the submitted matrix, not of an
	// unknown unrounded particle distribution.
	Lane4 pivot_column = mInertia.GetColumn4(0);
	double p = pivot_column[0];
	auto select_pivot = [&](int inAxis, Lane4Arg inColumn)
	{
		if (inColumn[inAxis] > p)
		{
			p = inColumn[inAxis];
			pivot_column = inColumn;
		}
	};
	select_pivot(1, mInertia.GetColumn4(1));
	select_pivot(2, mInertia.GetColumn4(2));
	select_pivot(3, mInertia.GetColumn4(3));
	bool rank_one = p > 0.0 && std::isfinite(p);
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			rank_one = rank_one && std::isfinite(mInertia(i, j))
				&& mInertia(i, j) == mInertia(j, i)
				&& double(mInertia(i, j)) * p == double(pivot_column[i]) * pivot_column[j];
	if (rank_one)
	{
		double trace = double(mInertia(0, 0)) + mInertia(1, 1) + mInertia(2, 2) + mInertia(3, 3);
		// Stable Householder completion of the occupied direction. Its first
		// column is +/- the normalized pivot column; the remaining columns
		// span its orthogonal complement before rounding. Choose the sign
		// to avoid cancellation.
		double w[4], norm_sq = 0.0;
		double scale = 1.0 / std::sqrt(p * trace);
		double sign = pivot_column[0] < 0.0f? -1.0 : 1.0;
		for (int i = 0; i < 4; ++i)
		{
			w[i] = sign * pivot_column[i] * scale + (i == 0? 1.0 : 0.0);
			norm_sq += w[i] * w[i];
		}
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				rotation(i, j) = float((i == j? 1.0 : 0.0) - 2.0 * w[i] * w[j] / norm_sq);
		diagonal[0] = float(trace);
	}
	else
	{
		// Using eigendecomposition on the 4x4 symmetric second moment tensor.
		Matrix<4, 4> inertia;
		inertia.CopyPart(mInertia, 0, 0, 4, 4, 0, 0);
		Matrix<4, 4> eigen_vec = Matrix<4, 4>::sIdentity();
		Vector<4> eigen_val;
		if (!EigenValueSymmetric(inertia, eigen_vec, eigen_val))
			return false;

		// Sort so that the largest eigenvalue goes first
		int indices[] = { 0, 1, 2, 3 };
		InsertionSort(indices, indices + 4, [&eigen_val](int inLeft, int inRight) { return eigen_val[inLeft] > eigen_val[inRight]; });
		for (int i = 0; i < 4; ++i)
		{
			const auto &col = eigen_vec.GetColumn(indices[i]);
			rotation.SetColumn4(i, Lane4(col[0], col[1], col[2], col[3]));
			diagonal[i] = eigen_val[indices[i]];
		}
	}

	// Ensure the rotation matrix is a proper rotation (determinant +1) rather than a reflection.
	// In 4D, flipping a single axis changes determinant sign. We flip the last column if det < 0.
	// sCross(c0, c1, c2).Dot(c3) has the OPPOSITE sign of the determinant in our
	// 4D cross convention (already -1 for the identity basis). A reflection cannot
	// be represented by GetRotor and silently corrupts a rotated inertia tensor.
	if (rotation.GetDeterminant() < 0.0f)
		rotation.SetColumn4(3, -rotation.GetColumn4(3));

	outRotation = rotation;
	outDiagonal = diagonal;

#ifdef JPH_ENABLE_ASSERTS
	// Validate: M ~= R * diag(D) * R^T
	Mat44 diag_mat = Mat44::sZero();
	for (int i = 0; i < 4; ++i)
		diag_mat(i, i) = diagonal[i];
	Mat44 reconstructed = outRotation * diag_mat * outRotation.Transposed();
	for (int i = 0; i < 4; ++i)
	{
		Vec4 col_orig = Vec4(mInertia.GetColumn4(i));
		Vec4 col_recon = Vec4(reconstructed.GetColumn4(i));
		JPH_ASSERT(col_recon.IsClose(col_orig, max(col_orig.LengthSq(), 1.0f) * 1.0e-6f));
	}
#endif

	return true;
}

Bivec MassProperties::sBivecInertiaFromSecondMomentDiagonal(Vec4Arg inDiagonal)
{
	// In the principal frame the bivector inertia is diagonal with I_ij = M_ii + M_jj
	float m1 = inDiagonal.GetX();
	float m2 = inDiagonal.GetY();
	float m3 = inDiagonal.GetZ();
	float m4 = inDiagonal.GetW();
	return Bivec(
		m1 + m2,   // I_12 (XY plane)
		m1 + m3,   // I_13 (XZ plane)
		m1 + m4,   // I_14 (XW plane)
		m2 + m3,   // I_23 (YZ plane)
		m2 + m4,   // I_24 (YW plane)
		m3 + m4);  // I_34 (ZW plane)
}

void MassProperties::SetMassAndInertiaOfSolidBox(Vec4Arg inBoxSize, float inDensity)
{
	// Mass = density * 4D hypervolume
	mMass = inBoxSize.GetX() * inBoxSize.GetY() * inBoxSize.GetZ() * inBoxSize.GetW() * inDensity;

	// For a solid tesseract with edges S_i and uniform density the second moment tensor is diagonal:
	//   M_ii = m * S_i^2 / 12
	Vec4 size_sq = inBoxSize * inBoxSize;
	Vec4 diag = size_sq * (mMass / 12.0f);

	mInertia = Mat44::sZero();
	mInertia(0, 0) = diag.GetX();
	mInertia(1, 1) = diag.GetY();
	mInertia(2, 2) = diag.GetZ();
	mInertia(3, 3) = diag.GetW();
}

void MassProperties::ScaleToMass(float inMass)
{
	if (mMass > 0.0f)
	{
		// Inertia scales linearly with mass
		float mass_scale = inMass / mMass;
		mMass = inMass;
		mInertia *= mass_scale;
	}
	else
	{
		// No prior inertia to scale — just set the mass
		mMass = inMass;
	}
}

Vec4 MassProperties::sGetEquivalentSolidBoxSize(float inMass, Vec4Arg inInertiaDiagonal)
{
	// Invert M_ii = m * S_i^2 / 12  →  S_i = sqrt(12 * M_i / m)
	Vec4 s_sq = inInertiaDiagonal * (12.0f / inMass);
	return Vec4(sqrt(s_sq.GetX()), sqrt(s_sq.GetY()), sqrt(s_sq.GetZ()), sqrt(s_sq.GetW()));
}

void MassProperties::Rotate(Mat44Arg inRotation)
{
	// Similarity transform: M' = R M R^T
	mInertia = inRotation * mInertia * inRotation.Transposed();
}

void MassProperties::Translate(Vec4Arg inTranslation)
{
	// Parallel axis theorem for the second moment tensor (assuming original COM at origin):
	//   M'_ij = M_ij + m * t_i * t_j
	// (cross terms involving the first moment vanish because COM is at origin)
	float t0 = inTranslation.GetX();
	float t1 = inTranslation.GetY();
	float t2 = inTranslation.GetZ();
	float t3 = inTranslation.GetW();

	// Add m * t t^T (rank-1 symmetric update)
	for (int i = 0; i < 4; ++i)
	{
		float ti = inTranslation[i];
		mInertia(i, 0) += mMass * ti * t0;
		mInertia(i, 1) += mMass * ti * t1;
		mInertia(i, 2) += mMass * ti * t2;
		mInertia(i, 3) += mMass * ti * t3;
	}
}

void MassProperties::Scale(Vec4Arg inScale)
{
	// Mass scales with 4D hypervolume. Negative components flip axes but mass stays positive.
	float mass_scale = abs(inScale.GetX() * inScale.GetY() * inScale.GetZ() * inScale.GetW());

	// Second moment tensor transforms under anisotropic scaling as:
	//   M'_ij = (s_i * s_j) * M_ij
	// plus an extra factor of mass_scale because the density integral scales with hypervolume.
	for (int i = 0; i < 4; ++i)
	{
		float si = inScale[i];
		for (int j = 0; j < 4; ++j)
		{
			float sj = inScale[j];
			mInertia(i, j) *= mass_scale * si * sj;
		}
	}

	mMass *= mass_scale;
}

void MassProperties::SaveBinaryState(StreamOut &inStream) const
{
	inStream.Write(mMass);
	inStream.Write(mInertia);
}

void MassProperties::RestoreBinaryState(StreamIn &inStream)
{
	inStream.Read(mMass);
	inStream.Read(mInertia);
}

JPH_NAMESPACE_END
