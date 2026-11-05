// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Mat44.h>
#include <Jolt/Math/Real.h>

JPH_NAMESPACE_BEGIN

/// 4D rigid body transform: rotation/scale (Mat44) + translation (RVec4).
///
/// In 4D, the transform decomposes into a 4x4 linear part (rotation, optionally with scale)
/// and a 4D translation vector. Applying to a point p: result = mRotation * p + mTranslation.
///
/// This replaces the original 3D Jolt's 4x4 homogeneous matrix (old Mat44/DMat44 with translation
/// in column 3). Method names like "Multiply3x3" are kept for backward compatibility even though
/// the rotation is now a full 4x4 matrix.
class alignas(max(JPH_VECTOR_ALIGNMENT, JPH_RVECTOR_ALIGNMENT)) RMat44
{
public:
	JPH_OVERRIDE_NEW_DELETE

	using ArgType = const RMat44 &;

	/// Constructors
								RMat44() = default;
	JPH_INLINE					RMat44(Mat44Arg inRotation, RVec4Arg inTranslation)				: mRotation(inRotation), mTranslation(inTranslation) { }
	JPH_INLINE explicit			RMat44(Mat44Arg inRotation)											: mRotation(inRotation), mTranslation(RVec4::sZero()) { }

#ifdef JPH_DOUBLE_PRECISION
	/// Accept float local translations without making DVec4's widening constructors implicit.
	JPH_INLINE RMat44(Mat44Arg inRotation, Lane4Arg inTranslation) : RMat44(inRotation, DVec4(inTranslation)) { }
	JPH_INLINE RMat44 PostTranslated(Lane4Arg inTranslation) const { return PostTranslated(DVec4(inTranslation)); }
	static JPH_INLINE RMat44 sTranslation(Lane4Arg inTranslation) { return sTranslation(DVec4(inTranslation)); }
	static JPH_INLINE RMat44 sRotationTranslation(RotorArg inRotor, Lane4Arg inTranslation) { return sRotationTranslation(inRotor, DVec4(inTranslation)); }
	static JPH_INLINE RMat44 sRotationTranslation(QuatArg inQuat, Lane4Arg inTranslation) { return sRotationTranslation(inQuat, DVec4(inTranslation)); }
	static JPH_INLINE RMat44 sInverseRotationTranslation(RotorArg inRotor, Lane4Arg inTranslation) { return sInverseRotationTranslation(inRotor, DVec4(inTranslation)); }
	static JPH_INLINE RMat44 sInverseRotationTranslation(QuatArg inQuat, Lane4Arg inTranslation) { return sInverseRotationTranslation(inQuat, DVec4(inTranslation)); }
#endif

	/// Zero transform
	static JPH_INLINE RMat44	sZero()
	{
		return RMat44(Mat44::sZero(), RVec4::sZero());
	}

	/// Identity transform (no rotation, no translation)
	static JPH_INLINE RMat44	sIdentity()
	{
		return RMat44(Mat44::sIdentity(), RVec4::sZero());
	}

	/// Translation-only transform (identity rotation)
	static JPH_INLINE RMat44	sTranslation(RVec4Arg inTranslation)
	{
		return RMat44(Mat44::sIdentity(), inTranslation);
	}

	/// Rotation-only transform from rotor (no translation)
	static JPH_INLINE RMat44	sRotation(RotorArg inRotor)
	{
		return RMat44(Mat44::sRotation(inRotor), RVec4::sZero());
	}

	/// Rotation-only transform from quaternion (no translation) — backward compat
	static JPH_INLINE RMat44	sRotation(QuatArg inQuat)
	{
		return RMat44(Mat44::sRotation(inQuat), RVec4::sZero());
	}

	/// Create from rotation (via rotor) and translation
	static JPH_INLINE RMat44	sRotationTranslation(RotorArg inRotor, RVec4Arg inTranslation)
	{
		return RMat44(Mat44::sRotation(inRotor), inTranslation);
	}

	/// Create from rotation (via quaternion) and translation — backward compat
	static JPH_INLINE RMat44	sRotationTranslation(QuatArg inQuat, RVec4Arg inTranslation)
	{
		return RMat44(Mat44::sRotation(inQuat), inTranslation);
	}

	/// Create inverse of a rotation-translation transform
	/// For orthogonal R and translation t: inverse maps p → R^T * (p - t) = R^T * p - R^T * t
	static JPH_INLINE RMat44	sInverseRotationTranslation(RotorArg inRotor, RVec4Arg inTranslation)
	{
		Mat44 rot_t = Mat44::sRotation(inRotor).Transposed();
		return RMat44(rot_t, -sMultiply(rot_t, inTranslation));
	}

	/// Create inverse of a rotation-translation transform from quaternion — backward compat
	static JPH_INLINE RMat44	sInverseRotationTranslation(QuatArg inQuat, RVec4Arg inTranslation)
	{
		Mat44 rot_t = Mat44::sRotation(inQuat).Transposed();
		return RMat44(rot_t, -sMultiply(rot_t, inTranslation));
	}

	/// Diagonal scale transform (no translation)
	static JPH_INLINE RMat44	sScale(Vec4Arg inScale)
	{
		Mat44 s(
			Lane4(inScale.GetX(), 0, 0, 0),
			Lane4(0, inScale.GetY(), 0, 0),
			Lane4(0, 0, inScale.GetZ(), 0),
			Lane4(0, 0, 0, inScale.GetW()));
		return RMat44(s, RVec4::sZero());
	}

	/// Get/set translation
	JPH_INLINE RVec4			GetTranslation() const												{ return mTranslation; }
	JPH_INLINE void				SetTranslation(RVec4Arg inTranslation)								{ mTranslation = inTranslation; }

	/// Get rotation/linear part as Mat44
	JPH_INLINE Mat44			GetRotation() const													{ return mRotation; }

	/// Get rotation column (0-3)
	JPH_INLINE Lane4			GetColumn4(int inCol) const											{ return mRotation.GetColumn4(inCol); }

	/// Get diagonal of the rotation/linear part
	JPH_INLINE Lane4			GetDiagonal4() const												{ return mRotation.GetDiagonal4(); }

	/// Extract rotor from the rotation matrix
	JPH_INLINE Rotor			GetRotor() const													{ return mRotation.GetRotor(); }

	/// Transform a single-precision point: R * inV + T
	JPH_INLINE RVec4			operator * (Vec4Arg inV) const
	{
		return RVec4(mRotation * Lane4(inV)) + mTranslation;
	}

#ifdef JPH_DOUBLE_PRECISION
	/// Transform a double-precision point: R * inV + T
	JPH_INLINE DVec4			operator * (DVec4Arg inV) const
	{
		return sMultiply(mRotation, inV) + mTranslation;
	}
#endif

	/// Compose two transforms: (R1, T1) * (R2, T2) = (R1*R2, R1*T2 + T1)
	JPH_INLINE RMat44			operator * (const RMat44 &inRHS) const
	{
		return RMat44(
			mRotation * inRHS.mRotation,
			sMultiply(mRotation, inRHS.mTranslation) + mTranslation);
	}

	/// Compose with Mat44 (rotation/scale only, no translation): (R1, T1) * R2 = (R1*R2, T1)
	JPH_INLINE RMat44			operator * (Mat44Arg inRHS) const
	{
		return RMat44(mRotation * inRHS, mTranslation);
	}

	/// Apply only the rotation/linear part to a vector (no translation)
	/// Named "Multiply3x3" for backward compatibility even though it's now a 4x4 rotation
	JPH_INLINE Vec4				Multiply3x3(Vec4Arg inV) const
	{
		return Vec4(mRotation * Lane4(inV));
	}

	/// Apply the transposed rotation/linear part to a vector (no translation)
	JPH_INLINE Vec4				Multiply3x3Transposed(Vec4Arg inV) const
	{
		return Vec4(mRotation.Transposed() * Lane4(inV));
	}

	/// Pre-translate: result transforms p as R*(p + inT) + T = R*p + (R*inT + T)
	JPH_INLINE RMat44			PreTranslated(Vec4Arg inTranslation) const
	{
		return RMat44(mRotation, RVec4(mRotation * Lane4(inTranslation)) + mTranslation);
	}

	/// Post-translate: result transforms p as R*p + (T + inT)
	JPH_INLINE RMat44			PostTranslated(RVec4Arg inTranslation) const
	{
		return RMat44(mRotation, mTranslation + inTranslation);
	}

	/// Pre-scale: scale the rotation columns. result = (R * diag(S), T)
	/// Each column i of R is multiplied by inScale[i]
	JPH_INLINE RMat44			PreScaled(Vec4Arg inScale) const
	{
		Mat44 scaled(
			Lane4(Vec4(mRotation.GetColumn4(0)) * Vec4::sReplicate(inScale.GetX())),
			Lane4(Vec4(mRotation.GetColumn4(1)) * Vec4::sReplicate(inScale.GetY())),
			Lane4(Vec4(mRotation.GetColumn4(2)) * Vec4::sReplicate(inScale.GetZ())),
			Lane4(Vec4(mRotation.GetColumn4(3)) * Vec4::sReplicate(inScale.GetW())));
		return RMat44(scaled, mTranslation);
	}

	/// Post-scale: each row i of the result is multiplied by inScale[i]
	/// result = (diag(S) * R, diag(S) * T)
	JPH_INLINE RMat44			PostScaled(Vec4Arg inScale) const
	{
		Mat44 scaled(
			Lane4(Vec4(mRotation.GetColumn4(0)) * inScale),
			Lane4(Vec4(mRotation.GetColumn4(1)) * inScale),
			Lane4(Vec4(mRotation.GetColumn4(2)) * inScale),
			Lane4(Vec4(mRotation.GetColumn4(3)) * inScale));
		return RMat44(scaled, mTranslation * RVec4(Lane4(inScale)));
	}

	/// Inverse of a rotation-translation transform (orthogonal R): (R^T, -R^T * T)
	JPH_INLINE RMat44			InversedRotationTranslation() const
	{
		Mat44 rot_t = mRotation.Transposed();
		return RMat44(rot_t, -sMultiply(rot_t, mTranslation));
	}

	/// General inverse (handles rotation with scale). Uses Mat44::Inversed for the linear part.
	JPH_INLINE RMat44			Inversed() const
	{
		Mat44 inv = mRotation.Inversed();
		return RMat44(inv, -sMultiply(inv, mTranslation));
	}

	/// Get the transposed rotation part as Mat44
	JPH_INLINE Mat44			Transposed3x3() const
	{
		return mRotation.Transposed();
	}

	/// Decompose into a rotation-translation (returned) and scale (written to outScale).
	/// The linear part is decomposed as R * diag(S).
	JPH_INLINE RMat44			Decompose(Vec4 &outScale) const
	{
		Vec4 col0(mRotation.GetColumn4(0));
		Vec4 col1(mRotation.GetColumn4(1));
		Vec4 col2(mRotation.GetColumn4(2));
		Vec4 col3(mRotation.GetColumn4(3));
		float s0 = col0.Length();
		float s1 = col1.Length();
		float s2 = col2.Length();
		float s3 = col3.Length();
		outScale = Vec4(s0, s1, s2, s3);

		Vec4 inv_scale = Vec4(
			s0 > 0.0f ? 1.0f / s0 : 0.0f,
			s1 > 0.0f ? 1.0f / s1 : 0.0f,
			s2 > 0.0f ? 1.0f / s2 : 0.0f,
			s3 > 0.0f ? 1.0f / s3 : 0.0f);
		Mat44 rot(
			Lane4(col0 * Vec4::sReplicate(inv_scale.GetX())),
			Lane4(col1 * Vec4::sReplicate(inv_scale.GetY())),
			Lane4(col2 * Vec4::sReplicate(inv_scale.GetZ())),
			Lane4(col3 * Vec4::sReplicate(inv_scale.GetW())));
		return RMat44(rot, mTranslation);
	}

	/// Convert to Mat44 — returns the rotation/linear part only (translation is lost)
	JPH_INLINE Mat44			ToMat44() const
	{
		return mRotation;
	}

	/// Equality check
	JPH_INLINE bool				operator == (const RMat44 &inRHS) const
	{
		return mRotation == inRHS.mRotation && mTranslation == inRHS.mTranslation;
	}

	/// Inequality check
	JPH_INLINE bool				operator != (const RMat44 &inRHS) const
	{
		return !(*this == inRHS);
	}

private:
	/// Apply the float linear map without narrowing double world coordinates.
	/// Matrix coefficients are floats; products and accumulation use position precision.
	static JPH_INLINE RVec4 sMultiply(Mat44Arg inMatrix, RVec4Arg inV)
	{
#ifdef JPH_DOUBLE_PRECISION
		return DVec4(inMatrix.GetColumn4(0)) * inV.GetX()
			+ DVec4(inMatrix.GetColumn4(1)) * inV.GetY()
			+ DVec4(inMatrix.GetColumn4(2)) * inV.GetZ()
			+ DVec4(inMatrix.GetColumn4(3)) * inV.GetW();
#else
		return Vec4(inMatrix * Lane4(inV));
#endif
	}

	Mat44						mRotation;															///< 4x4 rotation/scale matrix
	RVec4						mTranslation;														///< 4D translation (single or double precision)
};

/// Argument type for passing RMat44 by reference
using RMat44Arg = const RMat44 &;

JPH_NAMESPACE_END
