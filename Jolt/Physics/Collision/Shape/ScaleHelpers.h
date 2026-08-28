// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/PhysicsSettings.h>
#include <Jolt/Math/Vec4.h>
#include <Jolt/Math/Mat44.h>
#include <Jolt/Math/Rotor.h>

JPH_NAMESPACE_BEGIN

/// Helper functions to get properties of a 4D scaling vector
namespace ScaleHelpers
{
	/// Minimum valid scale value. This is used to prevent division by zero when scaling a shape with a zero scale.
	static constexpr float	cMinScale = 1.0e-6f;

	/// The tolerance used to check if components of the scale vector are the same
	static constexpr float	cScaleToleranceSq = 1.0e-8f;

	/// Test if a scale is identity
	inline bool				IsNotScaled(Vec4Arg inScale)									{ return inScale.IsClose(Vec4::sOne(), cScaleToleranceSq); }

	/// Test if a scale is uniform (all 4 components equal)
	inline bool				IsUniformScale(Vec4Arg inScale)
	{
		float x = inScale.GetX();
		Vec4 rep = Vec4::sReplicate(x);
		return inScale.IsClose(rep, cScaleToleranceSq);
	}

	/// Test if a scale is uniform in XZ (kept for legacy cylindrical shape helpers)
	inline bool				IsUniformScaleXZ(Vec4Arg inScale)
	{
		return fabsf(inScale.GetX() - inScale.GetZ()) * fabsf(inScale.GetX() - inScale.GetZ()) < cScaleToleranceSq;
	}

	/// Scale the convex radius of an object (uses the smallest absolute scale component)
	inline float			ScaleConvexRadius(float inConvexRadius, Vec4Arg inScale)		{ return min(inConvexRadius * inScale.Abs().ReduceMin(), cDefaultConvexRadius); }

	/// Test if a scale flips an object inside out (which requires flipping all normals and winding).
	/// In 4D, this is true if an odd number of the 4 components are negative.
	inline bool				IsInsideOut(Vec4Arg inScale)									{ return (CountBits(Vec4::sLess(inScale, Vec4::sZero()).GetTrues() & 0xF) & 1) != 0; }

	/// Test if any of the components of the scale have a value below cMinScale
	inline bool				IsZeroScale(Vec4Arg inScale)									{ return Vec4::sLess(inScale.Abs(), Vec4::sReplicate(cMinScale)).TestAnyTrue(); }

	/// Ensure that the scale for each component is at least cMinScale
	inline Vec4				MakeNonZeroScale(Vec4Arg inScale)								{ return inScale.GetSign() * Vec4::sMax(inScale.Abs(), Vec4::sReplicate(cMinScale)); }

	/// Get the average scale, used to make the scale uniform when a shape doesn't support non-uniform scale
	inline Vec4				MakeUniformScale(Vec4Arg inScale)
	{
		float avg = 0.25f * (inScale.GetX() + inScale.GetY() + inScale.GetZ() + inScale.GetW());
		return Vec4::sReplicate(avg);
	}

	/// Average the scale in XZ, used to make the scale uniform when a shape doesn't support non-uniform scale in the XZ plane
	inline Vec4				MakeUniformScaleXZ(Vec4Arg inScale)
	{
		float avg_xz = 0.5f * (inScale.GetX() + inScale.GetZ());
		return Vec4(avg_xz, inScale.GetY(), avg_xz, inScale.GetW());
	}

	/// Test if a (non-uniform) scale can be passed through a rotation and remain axis-aligned.
	/// inScale is a scale in local space of the shape, so the transform (ignoring translation) is T = diag(inScale) * R.
	/// To push the scale to the child we need T = R * diag(ChildScale), i.e. ChildScale = R^T * diag(inScale) * R.
	/// This is only a pure (axis-aligned) scale if the off-diagonal elements of R^T * diag(inScale) * R vanish.
	inline bool				CanScaleBeRotated(RotorArg inRotation, Vec4Arg inScale)
	{
		Mat44 r = Mat44::sRotation(inRotation);
		Lane4 s(inScale);
		Lane4 c0 = r.GetColumn4(0), c1 = r.GetColumn4(1), c2 = r.GetColumn4(2), c3 = r.GetColumn4(3);
		// child_scale(i,j) = dot(c_i, diag(s) * c_j); require all off-diagonal entries ~ 0.
		Lane4 sc1 = s * c1, sc2 = s * c2, sc3 = s * c3;
		const float epsilon = 1.0e-6f;
		return abs(c0.Dot(sc1)) < epsilon && abs(c0.Dot(sc2)) < epsilon && abs(c0.Dot(sc3)) < epsilon
			&& abs(c1.Dot(sc2)) < epsilon && abs(c1.Dot(sc3)) < epsilon
			&& abs(c2.Dot(sc3)) < epsilon;
	}

	/// Get the diagonal of R^T * diag(inScale) * R (see comment at CanScaleBeRotated); the scale expressed in child space.
	inline Vec4				RotateScale(RotorArg inRotation, Vec4Arg inScale)
	{
		Mat44 r = Mat44::sRotation(inRotation);
		Lane4 s(inScale);
		Lane4 c0 = r.GetColumn4(0), c1 = r.GetColumn4(1), c2 = r.GetColumn4(2), c3 = r.GetColumn4(3);
		// diag(i) = dot(c_i, diag(s) * c_i) = sum_k s[k] * c_i[k]^2
		return Vec4(c0.Dot(s * c0), c1.Dot(s * c1), c2.Dot(s * c2), c3.Dot(s * c3));
	}
}

JPH_NAMESPACE_END
