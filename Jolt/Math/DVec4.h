// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Double3.h>
#include <Jolt/Math/Double4.h>

JPH_NAMESPACE_BEGIN

/// Four independent spatial components, stored in double precision.
/// Legacy XYZ constructors embed in R4 with W = 0; no lane is padding.
class [[nodiscard]] alignas(JPH_DVECTOR_ALIGNMENT) DVec4
{
public:
	JPH_OVERRIDE_NEW_DELETE

	// Underlying vector type
#if defined(JPH_USE_AVX)
	using Type = __m256d;
	using TypeArg = __m256d;
#elif defined(JPH_USE_SSE)
	using Type = struct { __m128d mLow, mHigh; };
	using TypeArg = const Type &;
#elif defined(JPH_USE_NEON)
	using Type = float64x2x2_t;
	using TypeArg = const Type &;
#else
	using Type = struct { double mData[4]; };
	using TypeArg = const Type &;
#endif

	// Argument type
	using ArgType = DVec4Arg;

	/// Constructor
								DVec4() = default; ///< Intentionally not initialized for performance reasons
								DVec4(const DVec4 &inRHS) = default;
	DVec4 &						operator = (const DVec4 &inRHS) = default;
	JPH_INLINE explicit			DVec4(Vec3Arg inRHS);
	JPH_INLINE explicit			DVec4(Lane4Arg inRHS);
	JPH_INLINE					DVec4(TypeArg inRHS) : mValue(inRHS)			{ }

	/// Embed 3 components with W = 0 (legacy compatibility)
	JPH_INLINE					DVec4(double inX, double inY, double inZ);

	/// Create a vector from 4 components (4D spatial)
	JPH_INLINE					DVec4(double inX, double inY, double inZ, double inW);

	/// Load 3 doubles from memory
	explicit JPH_INLINE			DVec4(const Double3 &inV);

	/// Load 4 doubles from memory
	explicit JPH_INLINE			DVec4(const Double4 &inV);

	/// Vector with all zeros
	static JPH_INLINE DVec4		sZero();

	/// Vector with all ones
	static JPH_INLINE DVec4		sOne();

	/// Vectors with the principal axis
	static JPH_INLINE DVec4		sAxisX()										{ return DVec4(1, 0, 0); }
	static JPH_INLINE DVec4		sAxisY()										{ return DVec4(0, 1, 0); }
	static JPH_INLINE DVec4		sAxisZ()										{ return DVec4(0, 0, 1); }
	static JPH_INLINE DVec4		sAxisW()										{ return DVec4(0, 0, 0, 1); }

	/// Replicate inV across all components
	static JPH_INLINE DVec4		sReplicate(double inV);

	/// Vector with all NaN's
	static JPH_INLINE DVec4		sNaN();

	/// Load 3 doubles from memory (reads 64 bits extra which it doesn't use)
	static JPH_INLINE DVec4		sLoadDouble3Unsafe(const Double3 &inV);

	/// Store 3 doubles to memory
	JPH_INLINE void				StoreDouble3(Double3 *outV) const;

	/// Store 4 doubles to memory
	JPH_INLINE void				StoreDouble4(Double4 *outV) const;

	/// Explicit legacy XYZ projection (discards W), rounding to nearest
	JPH_INLINE explicit			operator Vec3() const;

	/// Convert to float vector 4 rounding to nearest (all 4 components are preserved, just narrowed to float)
	JPH_INLINE explicit			operator Vec4() const;

	/// Prepare all four coordinates for float conversion rounding towards zero
	JPH_INLINE DVec4			PrepareRoundToZero() const;

	/// Prepare all four coordinates for float conversion rounding away from zero
	JPH_INLINE DVec4			PrepareRoundToInf() const;

	/// Compatibility alias: round all four coordinates down to Vec4
	JPH_INLINE Vec4				ToVec3RoundDown() const;

	/// Compatibility alias: round all four coordinates up to Vec4
	JPH_INLINE Vec4				ToVec3RoundUp() const;

	/// Directed rounding of all four spatial components.
	JPH_INLINE Vec4 ToVec4RoundDown() const;
	JPH_INLINE Vec4 ToVec4RoundUp() const;

	/// Return the minimum value of each of the components
	static JPH_INLINE DVec4		sMin(DVec4Arg inV1, DVec4Arg inV2);

	/// Return the maximum of each of the components
	static JPH_INLINE DVec4		sMax(DVec4Arg inV1, DVec4Arg inV2);

	/// Clamp a vector between min and max (component wise)
	static JPH_INLINE DVec4		sClamp(DVec4Arg inV, DVec4Arg inMin, DVec4Arg inMax);

	/// Equals (component wise)
	static JPH_INLINE DVec4		sEquals(DVec4Arg inV1, DVec4Arg inV2);

	/// Less than (component wise)
	static JPH_INLINE DVec4		sLess(DVec4Arg inV1, DVec4Arg inV2);

	/// Less than or equal (component wise)
	static JPH_INLINE DVec4		sLessOrEqual(DVec4Arg inV1, DVec4Arg inV2);

	/// Greater than (component wise)
	static JPH_INLINE DVec4		sGreater(DVec4Arg inV1, DVec4Arg inV2);

	/// Greater than or equal (component wise)
	static JPH_INLINE DVec4		sGreaterOrEqual(DVec4Arg inV1, DVec4Arg inV2);

	/// Calculates inMul1 * inMul2 + inAdd
	static JPH_INLINE DVec4		sFusedMultiplyAdd(DVec4Arg inMul1, DVec4Arg inMul2, DVec4Arg inAdd);

	/// Component wise select, returns inNotSet when highest bit of inControl = 0 and inSet when highest bit of inControl = 1
	static JPH_INLINE DVec4		sSelect(DVec4Arg inNotSet, DVec4Arg inSet, DVec4Arg inControl);

	/// Logical or (component wise)
	static JPH_INLINE DVec4		sOr(DVec4Arg inV1, DVec4Arg inV2);

	/// Logical xor (component wise)
	static JPH_INLINE DVec4		sXor(DVec4Arg inV1, DVec4Arg inV2);

	/// Logical and (component wise)
	static JPH_INLINE DVec4		sAnd(DVec4Arg inV1, DVec4Arg inV2);

	/// Store if X is true in bit 0, Y in bit 1, Z in bit 2 and W in bit 3 (true is when highest bit of component is set)
	JPH_INLINE int				GetTrues() const;

	/// Test if any of the components are true (true is when highest bit of component is set)
	JPH_INLINE bool				TestAnyTrue() const;

	/// Test if all components are true (true is when highest bit of component is set)
	JPH_INLINE bool				TestAllTrue() const;

	/// Get individual components
#if defined(JPH_USE_AVX)
	JPH_INLINE double			GetX() const									{ return _mm_cvtsd_f64(_mm256_castpd256_pd128(mValue)); }
	JPH_INLINE double			GetY() const									{ return mF64[1]; }
	JPH_INLINE double			GetZ() const									{ return mF64[2]; }
	JPH_INLINE double			GetW() const									{ return mF64[3]; }
#elif defined(JPH_USE_SSE)
	JPH_INLINE double			GetX() const									{ return _mm_cvtsd_f64(mValue.mLow); }
	JPH_INLINE double			GetY() const									{ return mF64[1]; }
	JPH_INLINE double			GetZ() const									{ return _mm_cvtsd_f64(mValue.mHigh); }
	JPH_INLINE double			GetW() const									{ return mF64[3]; }
#elif defined(JPH_USE_NEON)
	JPH_INLINE double			GetX() const									{ return vgetq_lane_f64(mValue.val[0], 0); }
	JPH_INLINE double			GetY() const									{ return vgetq_lane_f64(mValue.val[0], 1); }
	JPH_INLINE double			GetZ() const									{ return vgetq_lane_f64(mValue.val[1], 0); }
	JPH_INLINE double			GetW() const									{ return vgetq_lane_f64(mValue.val[1], 1); }
#else
	JPH_INLINE double			GetX() const									{ return mF64[0]; }
	JPH_INLINE double			GetY() const									{ return mF64[1]; }
	JPH_INLINE double			GetZ() const									{ return mF64[2]; }
	JPH_INLINE double			GetW() const									{ return mF64[3]; }
#endif

	/// Set individual components
	JPH_INLINE void				SetX(double inX)								{ mF64[0] = inX; }
	JPH_INLINE void				SetY(double inY)								{ mF64[1] = inY; }
	JPH_INLINE void				SetZ(double inZ)								{ mF64[2] = inZ; }
	JPH_INLINE void				SetW(double inW)								{ mF64[3] = inW; }

	JPH_INLINE void Set(double inX, double inY, double inZ, double inW) { *this = DVec4(inX, inY, inZ, inW); }

	/// Set XYZ components, embedding with W = 0
	JPH_INLINE void				Set(double inX, double inY, double inZ)			{ *this = DVec4(inX, inY, inZ); }

	/// Get double component by index
	JPH_INLINE double			operator [] (uint inCoordinate) const			{ JPH_ASSERT(inCoordinate < 4); return mF64[inCoordinate]; }

	/// Set double component by index
	JPH_INLINE void				SetComponent(uint inCoordinate, double inValue)	{ JPH_ASSERT(inCoordinate < 4); mF64[inCoordinate] = inValue; }

	/// Comparison
	JPH_INLINE bool				operator == (DVec4Arg inV2) const;
	JPH_INLINE bool				operator != (DVec4Arg inV2) const				{ return !(*this == inV2); }

	/// Test if two vectors are close
	JPH_INLINE bool				IsClose(DVec4Arg inV2, double inMaxDistSq = 1.0e-24) const;

	/// Test if vector is near zero
	JPH_INLINE bool				IsNearZero(double inMaxDistSq = 1.0e-24) const;

	/// Test if vector is normalized
	JPH_INLINE bool				IsNormalized(double inTolerance = 1.0e-12) const;

	/// Test if vector contains NaN elements
	JPH_INLINE bool				IsNaN() const;

	/// Multiply two double vectors (component wise)
	JPH_INLINE DVec4			operator * (DVec4Arg inV2) const;

	/// Multiply vector with double
	JPH_INLINE DVec4			operator * (double inV2) const;

	/// Multiply vector with double
	friend JPH_INLINE DVec4		operator * (double inV1, DVec4Arg inV2);

	/// Divide vector by double
	JPH_INLINE DVec4			operator / (double inV2) const;

	/// Multiply vector with double
	JPH_INLINE DVec4 &			operator *= (double inV2);

	/// Multiply vector with vector
	JPH_INLINE DVec4 &			operator *= (DVec4Arg inV2);

	/// Divide vector by double
	JPH_INLINE DVec4 &			operator /= (double inV2);

	/// Add two vectors (component wise)
	JPH_INLINE DVec4			operator + (Lane4Arg inV2) const;

	/// Add two double vectors (component wise)
	JPH_INLINE DVec4			operator + (DVec4Arg inV2) const;

	/// Add two vectors (component wise)
	JPH_INLINE DVec4 &			operator += (Lane4Arg inV2);

	/// Add two double vectors (component wise)
	JPH_INLINE DVec4 &			operator += (DVec4Arg inV2);

	/// Negate
	JPH_INLINE DVec4			operator - () const;

	/// Subtract two vectors (component wise)
	JPH_INLINE DVec4			operator - (Lane4Arg inV2) const;

	/// Subtract two double vectors (component wise)
	JPH_INLINE DVec4			operator - (DVec4Arg inV2) const;

	/// Subtract two vectors (component wise)
	JPH_INLINE DVec4 &			operator -= (Lane4Arg inV2);

	/// Subtract two vectors (component wise)
	JPH_INLINE DVec4 &			operator -= (DVec4Arg inV2);

	/// Explicit legacy XYZ operand embedding (W = 0).
	JPH_INLINE DVec4 operator + (Vec3Arg inV2) const { return *this + DVec4(inV2); }

	/// Explicit legacy XYZ operand embedding (W = 0).
	JPH_INLINE DVec4 operator - (Vec3Arg inV2) const { return *this - DVec4(inV2); }

	/// Explicit legacy XYZ operand embedding (W = 0).
	JPH_INLINE DVec4 & operator += (Vec3Arg inV2) { return *this += DVec4(inV2); }

	/// Explicit legacy XYZ operand embedding (W = 0).
	JPH_INLINE DVec4 & operator -= (Vec3Arg inV2) { return *this -= DVec4(inV2); }

	/// Divide (component wise)
	JPH_INLINE DVec4			operator / (DVec4Arg inV2) const;

	/// Return the absolute value of each of the components
	JPH_INLINE DVec4			Abs() const;

	/// Reciprocal vector (1 / value) for each of the components
	JPH_INLINE DVec4			Reciprocal() const;

	/// Cross product (3D, ignores W)
	JPH_INLINE DVec4			Cross(DVec4Arg inV2) const;

	/// 4D Cross product: returns a vector perpendicular to all three input vectors.
	static JPH_INLINE DVec4		sCross(DVec4Arg inA, DVec4Arg inB, DVec4Arg inC);

	/// Dot product
	JPH_INLINE double			Dot(DVec4Arg inV2) const;

	/// Squared length of vector
	JPH_INLINE double			LengthSq() const;

	/// Length of vector
	JPH_INLINE double			Length() const;

	/// Normalize vector
	JPH_INLINE DVec4			Normalized() const;

	/// Component wise square root
	JPH_INLINE DVec4			Sqrt() const;

	/// Get vector that contains the sign of each element (returns 1 if positive, -1 if negative)
	JPH_INLINE DVec4			GetSign() const;

	/// To String
	friend ostream &			operator << (ostream &inStream, DVec4Arg inV)
	{
		inStream << inV.mF64[0] << ", " << inV.mF64[1] << ", " << inV.mF64[2] << ", " << inV.mF64[3];
		return inStream;
	}

	/// Legacy compatibility check: all four lanes are independent (no invariant).
	JPH_INLINE void				CheckW() const;

	/// Legacy compatibility helper: returns all four lanes unchanged.
	static JPH_INLINE Type		sFixW(TypeArg inValue);

	/// Representations of true and false for boolean operations
	inline static const double	cTrue = BitCast<double>(~uint64(0));
	inline static const double	cFalse = 0.0;

	union
	{
		Type					mValue;
		double					mF64[4];
	};
};

static_assert(std::is_trivial<DVec4>(), "Is supposed to be a trivial type!");

JPH_NAMESPACE_END

#include "DVec4.inl"
