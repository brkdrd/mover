// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Math/Vec3.h>
#include <Jolt/Core/HashCombine.h>

// Create a std::hash/JPH::Hash for Vec4
JPH_MAKE_HASHABLE(JPH::Vec4, t.GetX(), t.GetY(), t.GetZ(), t.GetW())

JPH_NAMESPACE_BEGIN

Vec4::Vec4(const Float4 &inV)
{
#if defined(JPH_USE_SSE)
	mValue = _mm_loadu_ps(&inV.x);
#elif defined(JPH_USE_NEON)
	mValue = vld1q_f32(&inV.x);
#elif defined(JPH_USE_RVV)
	const vfloat32m1_t v = __riscv_vle32_v_f32m1(&inV.x, 4);
	__riscv_vse32_v_f32m1(mF32, v, 4);
#else
	mF32[0] = inV.x;
	mF32[1] = inV.y;
	mF32[2] = inV.z;
	mF32[3] = inV.w;
#endif
}

Vec4::Vec4(Vec3Arg inRHS, float inW)
{
#if defined(JPH_USE_SSE4_1)
	mValue = _mm_blend_ps(inRHS.mValue, _mm_set1_ps(inW), 8);
#elif defined(JPH_USE_NEON)
	mValue = vsetq_lane_f32(inW, inRHS.mValue, 3);
#elif defined(JPH_USE_RVV)
	const vfloat32m1_t v = __riscv_vle32_v_f32m1(inRHS.mF32, 4);
	__riscv_vse32_v_f32m1(mF32, v, 4);
	mF32[3] = inW;
#else
	for (int i = 0; i < 3; i++)
		mF32[i] = inRHS.mF32[i];
	mF32[3] = inW;
#endif
}

Vec4 Vec4::sCross(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC)
{
	// 4D cross product: given 3 vectors A, B, C in R^4, returns a vector perpendicular to all three.
	// Computed as cofactor expansion of:
	// | e1   e2   e3   e4  |
	// | a.x  a.y  a.z  a.w |
	// | b.x  b.y  b.z  b.w |
	// | c.x  c.y  c.z  c.w |
	// Evaluate the four 3x3 cofactors in parallel. Each lane keeps the scalar
	// expression's multiplication and addition order; only Y and W change sign.
	Lane4 aj = inA.Swizzle<1, 0, 0, 0>(), ak = inA.Swizzle<2, 2, 1, 1>(), al = inA.Swizzle<3, 3, 3, 2>();
	Lane4 bj = inB.Swizzle<1, 0, 0, 0>(), bk = inB.Swizzle<2, 2, 1, 1>(), bl = inB.Swizzle<3, 3, 3, 2>();
	Lane4 cj = inC.Swizzle<1, 0, 0, 0>(), ck = inC.Swizzle<2, 2, 1, 1>(), cl = inC.Swizzle<3, 3, 3, 2>();
	Lane4 cofactors = aj * (bk * cl - bl * ck) - ak * (bj * cl - bl * cj) + al * (bj * ck - bk * cj);
	return Vec4(Lane4::sXor(cofactors, UVec4(0, 0x80000000U, 0, 0x80000000U).ReinterpretAsFloat()));
}

Vec4 Vec4::GetNormalizedPerpendicular() const
{
	// Find the component with the smallest absolute value
	int min_idx = Abs().GetLowestComponentIndex();

	// Create a vector with 1 in that component position
	Vec4 other = sZero();
	other.mF32[min_idx] = 1.0f;

	// Use Gram-Schmidt: subtract the projection of other onto this vector
	Vec4 result = other - *this * (Dot(other) / LengthSq());
	return result.Normalized();
}

void Vec4::GetNormalizedPerpendicularBasis(Vec4 &outPerp1, Vec4 &outPerp2, Vec4 &outPerp3) const
{
	outPerp1 = GetNormalizedPerpendicular();

	// Pick the coordinate axis least aligned with {*this, outPerp1} to seed the next perpendicular, so
	// that the generalized cross product below is well conditioned
	Vec4 align = Abs() + outPerp1.Abs();
	int best = 0;
	float best_v = align.GetX();
	if (align.GetY() < best_v) { best = 1; best_v = align.GetY(); }
	if (align.GetZ() < best_v) { best = 2; best_v = align.GetZ(); }
	if (align.GetW() < best_v) { best = 3; }
	Vec4 seed = best == 0? sAxisX() : best == 1? sAxisY() : best == 2? sAxisZ() : sAxisW();

	// sCross(a, b, c) is perpendicular to all three of its arguments
	outPerp2 = sCross(*this, outPerp1, seed).Normalized();
	outPerp3 = sCross(*this, outPerp1, outPerp2).Normalized();
}

void Vec4::GetContinuousPerpendicularBasis(Vec4 &outPerp1, Vec4 &outPerp2, Vec4 &outPerp3) const
{
	JPH_ASSERT(IsNormalized(1.0e-5f));
	// S^3 admits a continuous tangent frame. These three signed permutations are
	// mutually orthogonal, perpendicular to n, and have |n|, without normalization.
	// Unlike choosing a least-aligned coordinate axis, this does not abruptly rotate
	// the friction frame when a nearly axis-aligned contact normal changes slightly.
	outPerp1 = Lane4::sXor(Swizzle<1, 0, 3, 2>(), UVec4(0x80000000U, 0, 0, 0x80000000U).ReinterpretAsFloat());
	outPerp2 = Lane4::sXor(Swizzle<2, 3, 0, 1>(), UVec4(0x80000000U, 0x80000000U, 0, 0).ReinterpretAsFloat());
	outPerp3 = Lane4::sXor(Swizzle<3, 2, 1, 0>(), UVec4(0x80000000U, 0, 0x80000000U, 0).ReinterpretAsFloat());
}

JPH_NAMESPACE_END
