// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Geometry/OrientedBox.h>
#include <Jolt/Geometry/AABox.h>

JPH_NAMESPACE_BEGIN

/// Helper: test a single separating axis. Returns true if the boxes are separated along this axis.
static JPH_INLINE bool sTestSeparatingAxis(Vec4Arg inAxis, Vec4Arg inD, const Vec4 *inAxesA, Vec4Arg inHalfA, const Vec4 *inAxesB, Vec4Arg inHalfB, float inEpsilon)
{
	float len_sq = inAxis.LengthSq();
	if (len_sq < 1.0e-12f)
		return false; // Degenerate axis, cannot separate

	float ra = 0, rb = 0;
	for (int i = 0; i < 4; i++)
	{
		ra += abs(inAxis.Dot(inAxesA[i])) * inHalfA[i];
		rb += abs(inAxis.Dot(inAxesB[i])) * inHalfB[i];
	}
	return abs(inAxis.Dot(inD)) > ra + rb + inEpsilon;
}

bool OrientedBox::Overlaps(const AABox &inBox, float inEpsilon) const
{
	// 4D SAT: OBB vs AABox requires 56 separating axes
	// 4 face normals of A (standard basis), 4 face normals of B (OBB axes),
	// 24 edge-face axes (A edge x B face pair), 24 face-edge axes (B edge x A face pair)

	Vec4 a_center = inBox.GetCenter();
	Vec4 a_half = inBox.GetExtent();

	// A axes = standard basis (identity columns)
	Vec4 a_axes[4] = { Vec4(1, 0, 0, 0), Vec4(0, 1, 0, 0), Vec4(0, 0, 1, 0), Vec4(0, 0, 0, 1) };

	// B axes = OBB rotation columns
	Vec4 b_axes[4] = { Vec4(mOrientation.GetColumn4(0)), Vec4(mOrientation.GetColumn4(1)), Vec4(mOrientation.GetColumn4(2)), Vec4(mOrientation.GetColumn4(3)) };

	// Translation vector between centers
	Vec4 d = mTranslation - a_center;

	// Test 4 face normals of A (standard basis axes)
	for (int i = 0; i < 4; i++)
		if (sTestSeparatingAxis(a_axes[i], d, a_axes, a_half, b_axes, mHalfExtents, inEpsilon))
			return false;

	// Test 4 face normals of B (OBB axes)
	for (int i = 0; i < 4; i++)
		if (sTestSeparatingAxis(b_axes[i], d, a_axes, a_half, b_axes, mHalfExtents, inEpsilon))
			return false;

	// Test 24 edge-face axes: sCross(A_edge_i, B_edge_j, B_edge_k) for j < k
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int k = j + 1; k < 4; k++)
			{
				Vec4 axis = Vec4::sCross(a_axes[i], b_axes[j], b_axes[k]);
				if (sTestSeparatingAxis(axis, d, a_axes, a_half, b_axes, mHalfExtents, inEpsilon))
					return false;
			}

	// Test 24 face-edge axes: sCross(B_edge_i, A_edge_j, A_edge_k) for j < k
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int k = j + 1; k < 4; k++)
			{
				Vec4 axis = Vec4::sCross(b_axes[i], a_axes[j], a_axes[k]);
				if (sTestSeparatingAxis(axis, d, a_axes, a_half, b_axes, mHalfExtents, inEpsilon))
					return false;
			}

	// No separating axis found
	return true;
}

bool OrientedBox::Overlaps(const OrientedBox &inBox, float inEpsilon) const
{
	// 4D SAT: OBB vs OBB requires 56 separating axes

	// A axes
	Vec4 a_axes[4] = { Vec4(mOrientation.GetColumn4(0)), Vec4(mOrientation.GetColumn4(1)), Vec4(mOrientation.GetColumn4(2)), Vec4(mOrientation.GetColumn4(3)) };

	// B axes
	Vec4 b_axes[4] = { Vec4(inBox.mOrientation.GetColumn4(0)), Vec4(inBox.mOrientation.GetColumn4(1)), Vec4(inBox.mOrientation.GetColumn4(2)), Vec4(inBox.mOrientation.GetColumn4(3)) };

	// Translation vector between centers
	Vec4 d = inBox.mTranslation - mTranslation;

	// Test 4 face normals of A
	for (int i = 0; i < 4; i++)
		if (sTestSeparatingAxis(a_axes[i], d, a_axes, mHalfExtents, b_axes, inBox.mHalfExtents, inEpsilon))
			return false;

	// Test 4 face normals of B
	for (int i = 0; i < 4; i++)
		if (sTestSeparatingAxis(b_axes[i], d, a_axes, mHalfExtents, b_axes, inBox.mHalfExtents, inEpsilon))
			return false;

	// Test 24 edge-face axes: sCross(A_edge_i, B_edge_j, B_edge_k) for j < k
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int k = j + 1; k < 4; k++)
			{
				Vec4 axis = Vec4::sCross(a_axes[i], b_axes[j], b_axes[k]);
				if (sTestSeparatingAxis(axis, d, a_axes, mHalfExtents, b_axes, inBox.mHalfExtents, inEpsilon))
					return false;
			}

	// Test 24 face-edge axes: sCross(B_edge_i, A_edge_j, A_edge_k) for j < k
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int k = j + 1; k < 4; k++)
			{
				Vec4 axis = Vec4::sCross(b_axes[i], a_axes[j], a_axes[k]);
				if (sTestSeparatingAxis(axis, d, a_axes, mHalfExtents, b_axes, inBox.mHalfExtents, inEpsilon))
					return false;
			}

	// No separating axis found
	return true;
}

JPH_NAMESPACE_END
