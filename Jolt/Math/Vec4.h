// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Lane4.h>

JPH_NAMESPACE_BEGIN

/// 4 component spatial vector (x, y, z, w). All 4 components are meaningful spatial dimensions.
/// Inherits from Lane4 (SIMD utility) for backward compatibility — all Lane4 methods are available.
/// Vec4-specific methods (4D cross product, 4D perpendicular, etc.) are added here.
class [[nodiscard]] alignas(JPH_VECTOR_ALIGNMENT) Vec4 : public Lane4
{
public:
	// Argument type
	using ArgType = Vec4Arg;

	/// Constructor
								Vec4() = default; ///< Intentionally not initialized for performance reasons
								Vec4(const Vec4 &inRHS) = default;
	Vec4 &						operator = (const Vec4 &inRHS) = default;
	JPH_INLINE					Vec4(Lane4Arg inRHS) : Lane4(inRHS)				{ }
	JPH_INLINE					Vec4(Vec3Arg inRHS) : Lane4(inRHS.mValue)			{ }
	JPH_INLINE					Vec4(Vec3Arg inRHS, float inW);
	JPH_INLINE					Vec4(Type inRHS) : Lane4(inRHS)						{ }

	/// Load 4 floats from memory
	explicit JPH_INLINE			Vec4(const Float4 &inV);

	/// Create a vector from 4 components
	JPH_INLINE					Vec4(float inX, float inY, float inZ, float inW) : Lane4(inX, inY, inZ, inW) { }

	/// Vectors with the principal axes (overrides to return Vec4)
	static JPH_INLINE Vec4		sAxisX()											{ return Vec4(1, 0, 0, 0); }
	static JPH_INLINE Vec4		sAxisY()											{ return Vec4(0, 1, 0, 0); }
	static JPH_INLINE Vec4		sAxisZ()											{ return Vec4(0, 0, 1, 0); }
	static JPH_INLINE Vec4		sAxisW()											{ return Vec4(0, 0, 0, 1); }

	/// 4D Cross product: returns a vector perpendicular to all three input vectors.
	/// Computed as the determinant expansion of a 4x4 matrix with basis vectors in the first row.
	static JPH_INLINE Vec4		sCross(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC);

	/// Get normalized vector that is perpendicular to this vector (4D)
	JPH_INLINE Vec4				GetNormalizedPerpendicular() const;

	/// Get an orthonormal basis of the 3D subspace perpendicular to this vector. In 3D the orthogonal
	/// complement of a direction is a plane (2 tangents); in 4D it is a hyperplane, so there are 3.
	/// Used for contact friction (the tangent space of a contact is 3D) and for planar shape patches.
	/// This vector must be normalized. The result is deterministic for a given input.
	JPH_INLINE void				GetNormalizedPerpendicularBasis(Vec4 &outPerp1, Vec4 &outPerp2, Vec4 &outPerp3) const;

	/// Continuous orthonormal tangent frame on S^3, for persistent friction directions.
	/// Requires a normalized vector; uses signed permutations without branches or square roots.
	JPH_INLINE void				GetContinuousPerpendicularBasis(Vec4 &outPerp1, Vec4 &outPerp2, Vec4 &outPerp3) const;

	/// To String
	friend ostream &			operator << (ostream &inStream, Vec4Arg inV)
	{
		inStream << inV.mF32[0] << ", " << inV.mF32[1] << ", " << inV.mF32[2] << ", " << inV.mF32[3];
		return inStream;
	}
};

static_assert(std::is_trivial<Vec4>(), "Is supposed to be a trivial type!");

JPH_NAMESPACE_END

#include "Vec4.inl"
