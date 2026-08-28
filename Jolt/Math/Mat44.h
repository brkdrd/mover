// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/MathTypes.h>

JPH_NAMESPACE_BEGIN

/// 4x4 orthogonal rotation matrix for SO(4) rotations.
class [[nodiscard]] alignas(JPH_VECTOR_ALIGNMENT) Mat44
{
public:
	JPH_OVERRIDE_NEW_DELETE

	// Underlying column type
	using Type = Lane4::Type;

	// Argument type
	using ArgType = Mat44Arg;

	/// Constructor
								Mat44() = default; ///< Intentionally not initialized for performance reasons
	JPH_INLINE					Mat44(Lane4Arg inC1, Lane4Arg inC2, Lane4Arg inC3, Lane4Arg inC4);
								Mat44(const Mat44 &inM2) = default;
	Mat44 &						operator = (const Mat44 &inM2) = default;
	JPH_INLINE					Mat44(Type inC1, Type inC2, Type inC3, Type inC4);

	/// Zero matrix
	static JPH_INLINE Mat44		sZero();

	/// Identity matrix
	static JPH_INLINE Mat44		sIdentity();

	/// Matrix filled with NaN's
	static JPH_INLINE Mat44		sNaN();

	/// Load 16 floats from memory
	static JPH_INLINE Mat44		sLoadFloat4x4(const Float4 *inV);

	/// Load 16 floats from memory, 16 bytes aligned
	static JPH_INLINE Mat44		sLoadFloat4x4Aligned(const Float4 *inV);

	/// Rotate from quaternion
	static JPH_INLINE Mat44		sRotation(QuatArg inQuat);

	/// Rotate from 4D rotor (builds full 4x4 orthogonal rotation matrix)
	static JPH_INLINE Mat44		sRotation(RotorArg inRotor);

	/// 4D rotation plane factories — simple rotations in each of the 6 coordinate planes
	static JPH_INLINE Mat44		sRotationXY(float inAngle);								///< Rotation in the e1-e2 plane
	static JPH_INLINE Mat44		sRotationXZ(float inAngle);								///< Rotation in the e1-e3 plane
	static JPH_INLINE Mat44		sRotationXW(float inAngle);								///< Rotation in the e1-e4 plane
	static JPH_INLINE Mat44		sRotationYZ(float inAngle);								///< Rotation in the e2-e3 plane
	static JPH_INLINE Mat44		sRotationYW(float inAngle);								///< Rotation in the e2-e4 plane
	static JPH_INLINE Mat44		sRotationZW(float inAngle);								///< Rotation in the e3-e4 plane

	/// Get float component by element index
	JPH_INLINE float			operator () (uint inRow, uint inColumn) const			{ JPH_ASSERT(inRow < 4); JPH_ASSERT(inColumn < 4); return mCol[inColumn].mF32[inRow]; }
	JPH_INLINE float &			operator () (uint inRow, uint inColumn)					{ JPH_ASSERT(inRow < 4); JPH_ASSERT(inColumn < 4); return mCol[inColumn].mF32[inRow]; }

	/// Comparison
	JPH_INLINE bool				operator == (Mat44Arg inM2) const;
	JPH_INLINE bool				operator != (Mat44Arg inM2) const						{ return !(*this == inM2); }

	/// Test if two matrices are close
	JPH_INLINE bool				IsClose(Mat44Arg inM2, float inMaxDistSq = 1.0e-12f) const;

	/// Multiply matrix by matrix
	JPH_INLINE Mat44			operator * (Mat44Arg inM) const;

	/// Multiply vector by matrix
	JPH_INLINE Lane4			operator * (Lane4Arg inV) const;

	/// Multiply matrix with float
	JPH_INLINE Mat44			operator * (float inV) const;
	friend JPH_INLINE Mat44		operator * (float inV, Mat44Arg inM)					{ return inM * inV; }

	/// Multiply matrix with float
	JPH_INLINE Mat44 &			operator *= (float inV);

	/// Per element addition of matrix
	JPH_INLINE Mat44			operator + (Mat44Arg inM) const;

	/// Negate
	JPH_INLINE Mat44			operator - () const;

	/// Per element subtraction of matrix
	JPH_INLINE Mat44			operator - (Mat44Arg inM) const;

	/// Per element addition of matrix
	JPH_INLINE Mat44 &			operator += (Mat44Arg inM);

	/// Access to the columns
	JPH_INLINE Lane4			GetColumn4(uint inCol) const							{ JPH_ASSERT(inCol < 4); return mCol[inCol]; }
	JPH_INLINE void				SetColumn4(uint inCol, Lane4Arg inV)					{ JPH_ASSERT(inCol < 4); mCol[inCol] = inV; }
	JPH_INLINE Lane4			GetDiagonal4() const									{ return Lane4(mCol[0][0], mCol[1][1], mCol[2][2], mCol[3][3]); }
	JPH_INLINE void				SetDiagonal4(Lane4Arg inV)								{ mCol[0][0] = inV.GetX(); mCol[1][1] = inV.GetY(); mCol[2][2] = inV.GetZ(); mCol[3][3] = inV.GetW(); }

	/// Store matrix to memory
	JPH_INLINE void				StoreFloat4x4(Float4 *outV) const;

	/// Transpose matrix
	JPH_INLINE Mat44			Transposed() const;

	/// Determinant of the 4x4 matrix. Zero (within tolerance) means Inversed() is not usable.
	JPH_INLINE float			GetDeterminant() const;

	/// Inverse 4x4 matrix
	JPH_INLINE Mat44			Inversed() const;

	/// Convert to 4D rotor (matrix must be orthogonal)
	JPH_INLINE Rotor			GetRotor() const;

	/// To String
	friend ostream &			operator << (ostream &inStream, Mat44Arg inM)
	{
		inStream << inM.mCol[0] << ", " << inM.mCol[1] << ", " << inM.mCol[2] << ", " << inM.mCol[3];
		return inStream;
	}

private:
	Lane4						mCol[4];												///< Column
};

static_assert(std::is_trivial<Mat44>(), "Is supposed to be a trivial type!");

JPH_NAMESPACE_END

#include "Mat44.inl"
