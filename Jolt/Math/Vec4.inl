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
	float ax = inA.mF32[0], ay = inA.mF32[1], az = inA.mF32[2], aw = inA.mF32[3];
	float bx = inB.mF32[0], by = inB.mF32[1], bz = inB.mF32[2], bw = inB.mF32[3];
	float cx = inC.mF32[0], cy = inC.mF32[1], cz = inC.mF32[2], cw = inC.mF32[3];

	// Cofactor for e1 (x): det of 3x3 minor excluding column 0
	float x =  ay * (bz * cw - bw * cz) - az * (by * cw - bw * cy) + aw * (by * cz - bz * cy);
	// Cofactor for e2 (y): det of 3x3 minor excluding column 1 (negated)
	float y = -(ax * (bz * cw - bw * cz) - az * (bx * cw - bw * cx) + aw * (bx * cz - bz * cx));
	// Cofactor for e3 (z): det of 3x3 minor excluding column 2
	float z =  ax * (by * cw - bw * cy) - ay * (bx * cw - bw * cx) + aw * (bx * cy - by * cx);
	// Cofactor for e4 (w): det of 3x3 minor excluding column 3 (negated)
	float w = -(ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx));

	return Vec4(x, y, z, w);
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

JPH_NAMESPACE_END
