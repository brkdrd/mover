// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Vec4.h>
#include <Jolt/Math/Vec8.h>

JPH_NAMESPACE_BEGIN

/// Bivector class for 4D physics: angular velocity, torque, angular momentum.
///
/// A bivector in R^4 has 6 components corresponding to the 6 coordinate planes:
///   Omega = w12*e12 + w13*e13 + w14*e14 + w23*e23 + w24*e24 + w34*e34
///
/// Component layout in Vec8: [e12, e13, e14, e23, e24, e34, 0, 0]
///   - e_ij represents the oriented plane spanned by basis vectors e_i and e_j
///   - Components 6 and 7 are always zero (padding for Vec8 SIMD)
///
/// In 4D rigid body dynamics:
///   - Angular velocity is a bivector (replaces Vec3 in 3D)
///   - Torque is a bivector: tau = F ^ r (wedge product, replaces cross product). Note the operand
///     order: it is the reverse of the 3D r x F, because Contract (below) is the negative adjoint
///     of the wedge. See UnitTests/Physics/DynamicsConventionsTest.cpp.
///   - Angular momentum is a bivector: L = I * Omega
///   - The inertia tensor maps bivectors to bivectors (6x6 matrix)
class [[nodiscard]] alignas(JPH_VECTOR_ALIGNMENT) Bivec
{
public:
	JPH_OVERRIDE_NEW_DELETE

	///@name Constructors
	///@{
	inline						Bivec() = default; ///< Intentionally not initialized for performance reasons
								Bivec(const Bivec &inRHS) = default;
	Bivec &						operator = (const Bivec &inRHS) = default;

	/// Construct from 6 explicit components [e12, e13, e14, e23, e24, e34]
	JPH_INLINE					Bivec(float inE12, float inE13, float inE14,
									  float inE23, float inE24, float inE34)
									: mValue(inE12, inE13, inE14, inE23, inE24, inE34, 0.0f, 0.0f) { }

	/// Construct from Vec8 (components 6,7 should be zero)
	JPH_INLINE explicit			Bivec(Vec8Arg inV) : mValue(inV) { }

	/// Construct from Float8 (components 6,7 should be zero)
	JPH_INLINE explicit			Bivec(const Float8 &inV) : mValue(inV) { }
	///@}

	///@name Get components
	///@{
	JPH_INLINE float			GetE12() const			{ return mValue[0]; }
	JPH_INLINE float			GetE13() const			{ return mValue[1]; }
	JPH_INLINE float			GetE14() const			{ return mValue[2]; }
	JPH_INLINE float			GetE23() const			{ return mValue[3]; }
	JPH_INLINE float			GetE24() const			{ return mValue[4]; }
	JPH_INLINE float			GetE34() const			{ return mValue[5]; }

	/// Get component by index (0=e12, 1=e13, 2=e14, 3=e23, 4=e24, 5=e34)
	JPH_INLINE float			operator [] (uint inIdx) const { JPH_ASSERT(inIdx < 6); return mValue[inIdx]; }

	/// Get the underlying Vec8
	JPH_INLINE Vec8				GetVec8() const			{ return mValue; }
	///@}

	///@name Default bivectors
	///@{
	JPH_INLINE static Bivec		sZero()					{ return Bivec(Vec8::sZero()); }
	///@}

	///@name Tests
	///@{
	JPH_INLINE bool				operator == (const Bivec &inRHS) const	{ return mValue == inRHS.mValue; }
	JPH_INLINE bool				operator != (const Bivec &inRHS) const	{ return mValue != inRHS.mValue; }

	/// If this bivector is close to zero
	JPH_INLINE bool				IsNearZero(float inMaxDistSq = 1.0e-12f) const { return LengthSq() <= inMaxDistSq; }

	/// If any component is NaN
	JPH_INLINE bool				IsNaN() const			{ return mValue.IsNaN(); }
	///@}

	///@name Arithmetic
	///@{
	JPH_INLINE Bivec			operator - () const		{ return Bivec(-mValue); }
	JPH_INLINE Bivec			operator + (const Bivec &inRHS) const { return Bivec(mValue + inRHS.mValue); }
	JPH_INLINE Bivec			operator - (const Bivec &inRHS) const { return Bivec(mValue - inRHS.mValue); }
	JPH_INLINE Bivec &			operator += (const Bivec &inRHS)	{ mValue = mValue + inRHS.mValue; return *this; }
	JPH_INLINE Bivec &			operator -= (const Bivec &inRHS)	{ mValue = mValue - inRHS.mValue; return *this; }

	/// Multiply by scalar
	JPH_INLINE Bivec			operator * (float inV) const		{ return Bivec(mValue * inV); }
	friend JPH_INLINE Bivec		operator * (float inV, const Bivec &inRHS) { return Bivec(inRHS.mValue * inV); }
	JPH_INLINE Bivec &			operator *= (float inV)				{ mValue = mValue * inV; return *this; }

	/// Component-wise multiply (useful for diagonal inertia: I_diag * omega)
	JPH_INLINE Bivec			operator * (const Bivec &inRHS) const { return Bivec(mValue * inRHS.mValue); }

	/// Dot product (sum of component products)
	JPH_INLINE float			Dot(const Bivec &inRHS) const		{ return mValue.Dot(inRHS.mValue); }

	/// Squared length
	JPH_INLINE float			LengthSq() const					{ return mValue.LengthSq(); }

	/// Length
	JPH_INLINE float			Length() const						{ return mValue.Length(); }

	/// Normalize
	JPH_INLINE Bivec			Normalized() const					{ return Bivec(mValue.Normalized()); }
	///@}

	///@name Bitwise operations (for DOF masking)
	///@{

	/// Bitwise AND (for masking with DOF flags)
	/// When a mask component = 0xffffffff (as float bits), result = inV1 component; when the mask
	/// component = 0x00000000, result = 0.
	JPH_INLINE static Bivec		sAnd(const Bivec &inV1, const Bivec &inV2)
	{
		// Bit-and per component through memcpy. A reinterpret_cast<uint32 *>(float *) here is
		// strict-aliasing UB and is miscompiled at -O2/-O3 (the masked write is dropped).
		Float8 f1, f2, out;
		inV1.mValue.StoreFloat8(&f1);
		inV2.mValue.StoreFloat8(&f2);
		for (int i = 0; i < 8; ++i)
		{
			uint32 a, b;
			memcpy(&a, &f1.mValue[i], sizeof(uint32));
			memcpy(&b, &f2.mValue[i], sizeof(uint32));
			uint32 r = a & b;
			memcpy(&out.mValue[i], &r, sizeof(uint32));
		}
		return Bivec(out);
	}

	///@}

	///@name 4D Physics operations
	///@{

	/// Wedge product of two vectors: a ^ b (creates a bivector)
	/// This is the 4D analog of the cross product for computing torque: tau = F ^ r. The 3D
	/// expression "a x b" in a torque / Jacobian role becomes "b ^ a" in 4D (operands swapped),
	/// because Contract is the negative adjoint of the wedge: <B _| r, v> == -<B, r ^ v>.
	JPH_INLINE static Bivec		sWedge(Vec4Arg inA, Vec4Arg inB)
	{
		float a1 = inA.GetX(), a2 = inA.GetY(), a3 = inA.GetZ(), a4 = inA.GetW();
		float b1 = inB.GetX(), b2 = inB.GetY(), b3 = inB.GetZ(), b4 = inB.GetW();
		return Bivec(
			a1 * b2 - a2 * b1,  // e12
			a1 * b3 - a3 * b1,  // e13
			a1 * b4 - a4 * b1,  // e14
			a2 * b3 - a3 * b2,  // e23
			a2 * b4 - a4 * b2,  // e24
			a3 * b4 - a4 * b3   // e34
		);
	}

	/// Left contraction with a vector: Omega _| v (returns a vector)
	/// In 4D, this replaces omega.Cross(r) for computing point velocity:
	///   v_point = v_linear + Omega._|_r
	///
	/// (Omega _| v)_i = sum_j Omega_ij * v_j   where Omega_ij is antisymmetric
	JPH_INLINE Vec4				Contract(Vec4Arg inV) const
	{
		float v1 = inV.GetX(), v2 = inV.GetY(), v3 = inV.GetZ(), v4 = inV.GetW();
		float e12 = mValue[0], e13 = mValue[1], e14 = mValue[2];
		float e23 = mValue[3], e24 = mValue[4], e34 = mValue[5];
		return Vec4(
			 e12 * v2 + e13 * v3 + e14 * v4,    // (Omega _| v).x
			-e12 * v1 + e23 * v3 + e24 * v4,    // (Omega _| v).y
			-e13 * v1 - e23 * v2 + e34 * v4,    // (Omega _| v).z
			-e14 * v1 - e24 * v2 - e34 * v3     // (Omega _| v).w
		);
	}

	///@}

	/// Store as Float8
	JPH_INLINE void				StoreFloat8(Float8 *outV) const	{ mValue.StoreFloat8(outV); }

	/// To String
	friend ostream &			operator << (ostream &inStream, const Bivec &inB)
	{
		inStream << inB.mValue[0] << "*e12 + " << inB.mValue[1] << "*e13 + " << inB.mValue[2] << "*e14 + "
				 << inB.mValue[3] << "*e23 + " << inB.mValue[4] << "*e24 + " << inB.mValue[5] << "*e34";
		return inStream;
	}

	/// 8 component storage: [e12, e13, e14, e23, e24, e34, 0, 0]
	Vec8						mValue;
};

static_assert(std::is_trivial<Bivec>(), "Is supposed to be a trivial type!");

JPH_NAMESPACE_END
