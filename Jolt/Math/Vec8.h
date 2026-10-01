// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Float8.h>
#include <Jolt/Math/MathTypes.h>

JPH_NAMESPACE_BEGIN

/// 8-component SIMD vector. Used internally by Rotor for 4D rotation algebra.
class [[nodiscard]] alignas(JPH_VECTOR_ALIGNMENT) Vec8
{
public:
	JPH_OVERRIDE_NEW_DELETE

	// Underlying vector type
#if defined(JPH_USE_AVX)
	using Type = __m256;
#elif defined(JPH_USE_SSE)
	using Type = struct { __m128 mLow; __m128 mHigh; };
#elif defined(JPH_USE_NEON)
	using Type = struct { float32x4_t mLow; float32x4_t mHigh; };
#else
	using Type = struct { float mData[8]; };
#endif

	// Argument type
	using ArgType = Vec8Arg;

	/// Constructor
								Vec8() = default; ///< Intentionally not initialized for performance reasons
								Vec8(const Vec8 &inRHS) = default;
	Vec8 &						operator = (const Vec8 &inRHS) = default;

	/// Construct from underlying type
	JPH_INLINE					Vec8(Type inRHS) : mValue(inRHS) { }

	/// Construct from two Lane4s (low = components 0-3, high = components 4-7)
	JPH_INLINE					Vec8(Lane4Arg inLow, Lane4Arg inHigh);

	/// Construct from 8 floats
	JPH_INLINE					Vec8(float in0, float in1, float in2, float in3,
									 float in4, float in5, float in6, float in7);

	/// Load from Float8 POD storage
	explicit JPH_INLINE			Vec8(const Float8 &inV);

	/// Vector with all zeros
	static JPH_INLINE Vec8		sZero();

	/// Replicate inV across all 8 components
	static JPH_INLINE Vec8		sReplicate(float inV);

	/// Get the low 4 components as Lane4
	JPH_INLINE Lane4			GetLow() const;

	/// Get the high 4 components as Lane4
	JPH_INLINE Lane4			GetHigh() const;

	/// Get component by index
	JPH_INLINE float			operator [] (uint inCoordinate) const			{ JPH_ASSERT(inCoordinate < 8); return mF32[inCoordinate]; }

	/// Set component by index
	JPH_INLINE void				SetComponent(uint inCoordinate, float inValue)	{ JPH_ASSERT(inCoordinate < 8); mF32[inCoordinate] = inValue; }

	/// Comparison
	JPH_INLINE bool				operator == (Vec8Arg inV2) const;
	JPH_INLINE bool				operator != (Vec8Arg inV2) const				{ return !(*this == inV2); }

	/// Test if two vectors are close
	JPH_INLINE bool				IsClose(Vec8Arg inV2, float inMaxDistSq = 1.0e-12f) const;

	/// Test if vector contains NaN elements
	JPH_INLINE bool				IsNaN() const;

	/// Add (component wise)
	JPH_INLINE Vec8				operator + (Vec8Arg inV2) const;

	/// Subtract (component wise)
	JPH_INLINE Vec8				operator - (Vec8Arg inV2) const;

	/// Negate
	JPH_INLINE Vec8				operator - () const;

	/// Multiply (component wise)
	/// Bitwise AND of two vectors, for masking. TODO(4D): added for Bivec DOF masking, which is on
	/// the hot path of every constraint solve; the scalar version it replaces went through memory.
	static JPH_INLINE Vec8		sAnd(Vec8Arg inV1, Vec8Arg inV2);

	JPH_INLINE Vec8				operator * (Vec8Arg inV2) const;

	/// Multiply by scalar
	JPH_INLINE Vec8				operator * (float inV2) const;

	/// Multiply by scalar (commutative)
	friend JPH_INLINE Vec8		operator * (float inV1, Vec8Arg inV2);

	/// Store to Float8 POD storage
	JPH_INLINE void				StoreFloat8(Float8 *outV) const;

	/// 8-component dot product
	JPH_INLINE float			Dot(Vec8Arg inV2) const;

	/// Squared norm (sum of squares of all 8 components)
	JPH_INLINE float			LengthSq() const;

	/// Length (norm of all 8 components)
	JPH_INLINE float			Length() const;

	/// Normalize
	JPH_INLINE Vec8				Normalized() const;

	/// Test if normalized
	JPH_INLINE bool				IsNormalized(float inTolerance = 1.0e-6f) const;

	/// To String
	friend ostream &			operator << (ostream &inStream, Vec8Arg inV)
	{
		inStream << inV.mF32[0] << ", " << inV.mF32[1] << ", " << inV.mF32[2] << ", " << inV.mF32[3]
				 << ", " << inV.mF32[4] << ", " << inV.mF32[5] << ", " << inV.mF32[6] << ", " << inV.mF32[7];
		return inStream;
	}

	union
	{
		Type					mValue;
		float					mF32[8];
	};
};

static_assert(std::is_trivial<Vec8>(), "Is supposed to be a trivial type!");

JPH_NAMESPACE_END

#include "Vec8.inl"
