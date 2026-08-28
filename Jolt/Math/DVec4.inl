// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Core/HashCombine.h>

// Create a std::hash/JPH::Hash for DVec4
JPH_MAKE_HASHABLE(JPH::DVec4, t.GetX(), t.GetY(), t.GetZ())

JPH_NAMESPACE_BEGIN

DVec4::DVec4(Vec3Arg inRHS)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_cvtps_pd(inRHS.mValue);
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_cvtps_pd(inRHS.mValue);
	mValue.mHigh = _mm_cvtps_pd(_mm_shuffle_ps(inRHS.mValue, inRHS.mValue, _MM_SHUFFLE(2, 2, 2, 2)));
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vcvt_f64_f32(vget_low_f32(inRHS.mValue));
	mValue.val[1] = vcvt_high_f64_f32(inRHS.mValue);
#elif defined(JPH_USE_RVV)
	const vfloat32m1_t src = __riscv_vle32_v_f32m1(inRHS.mF32, 3);
	const vfloat64m2_t widened = __riscv_vfwcvt_f_f_v_f64m2(src, 3);
	__riscv_vse64_v_f64m2(mF64, widened, 3);
#else
	mF64[0] = (double)inRHS.GetX();
	mF64[1] = (double)inRHS.GetY();
	mF64[2] = (double)inRHS.GetZ();
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
}

DVec4::DVec4(Lane4Arg inRHS) :
	DVec4(Vec3(inRHS))
{
}

DVec4::DVec4(double inX, double inY, double inZ)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_set_pd(inZ, inZ, inY, inX); // Assure Z and W are the same
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_set_pd(inY, inX);
	mValue.mHigh = _mm_set1_pd(inZ);
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vcombine_f64(vcreate_f64(BitCast<uint64>(inX)), vcreate_f64(BitCast<uint64>(inY)));
	mValue.val[1] = vdupq_n_f64(inZ);
#elif defined(JPH_USE_RVV)
	vfloat64m2_t v = __riscv_vfmv_v_f_f64m2(inZ, 4);
	v = __riscv_vfslide1up_vf_f64m2(v, inY, 4);
	v = __riscv_vfslide1up_vf_f64m2(v, inX, 4);
	__riscv_vse64_v_f64m2(mF64, v, 4);
#else
	mF64[0] = inX;
	mF64[1] = inY;
	mF64[2] = inZ;
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
}

DVec4::DVec4(double inX, double inY, double inZ, double inW)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_set_pd(inW, inZ, inY, inX);
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_set_pd(inY, inX);
	mValue.mHigh = _mm_set_pd(inW, inZ);
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vcombine_f64(vcreate_f64(BitCast<uint64>(inX)), vcreate_f64(BitCast<uint64>(inY)));
	mValue.val[1] = vcombine_f64(vcreate_f64(BitCast<uint64>(inZ)), vcreate_f64(BitCast<uint64>(inW)));
#else
	mF64[0] = inX;
	mF64[1] = inY;
	mF64[2] = inZ;
	mF64[3] = inW;
#endif
}

DVec4::DVec4(const Double3 &inV)
{
#if defined(JPH_USE_AVX)
	Type x = _mm256_castpd128_pd256(_mm_load_sd(&inV.x));
	Type y = _mm256_castpd128_pd256(_mm_load_sd(&inV.y));
	Type z = _mm256_broadcast_sd(&inV.z);
	Type xy = _mm256_unpacklo_pd(x, y);
	mValue = _mm256_blend_pd(xy, z, 0b1100); // Assure Z and W are the same
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_loadu_pd(&inV.x);
	mValue.mHigh = _mm_set1_pd(inV.z);
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vld1q_f64(&inV.x);
	mValue.val[1] = vdupq_n_f64(inV.z);
#elif defined(JPH_USE_RVV)
	vfloat64m2_t v = __riscv_vle64_v_f64m2(&inV.x, 3);
	__riscv_vse64_v_f64m2(mF64, v, 3);
#else
	mF64[0] = inV.x;
	mF64[1] = inV.y;
	mF64[2] = inV.z;
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
}

void DVec4::CheckW() const
{
#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
	// Avoid asserts when both components are NaN
	JPH_ASSERT(reinterpret_cast<const uint64 *>(mF64)[2] == reinterpret_cast<const uint64 *>(mF64)[3]);
#endif // JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
}

/// Internal helper function that ensures that the Z component is replicated to the W component to prevent divisions by zero
DVec4::Type DVec4::sFixW(TypeArg inValue)
{
#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
	#if defined(JPH_USE_AVX)
		return _mm256_shuffle_pd(inValue, inValue, 2);
	#elif defined(JPH_USE_SSE)
		Type value;
		value.mLow = inValue.mLow;
		value.mHigh = _mm_shuffle_pd(inValue.mHigh, inValue.mHigh, 0);
		return value;
	#elif defined(JPH_USE_NEON)
		Type value;
		value.val[0] = inValue.val[0];
		value.val[1] = vdupq_laneq_f64(inValue.val[1], 0);
		return value;
	#elif defined(JPH_USE_RVV)
		Type value;
		const vfloat64m2_t buffer = __riscv_vle64_v_f64m2(inValue.mData, 3);
		__riscv_vse64_v_f64m2(value.mData, buffer, 3);
		value.mData[3] = value.mData[2];
		return value;
	#else
		Type value;
		value.mData[0] = inValue.mData[0];
		value.mData[1] = inValue.mData[1];
		value.mData[2] = inValue.mData[2];
		value.mData[3] = inValue.mData[2];
		return value;
	#endif
#else
	return inValue;
#endif // JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
}

DVec4 DVec4::sZero()
{
#if defined(JPH_USE_AVX)
	return _mm256_setzero_pd();
#elif defined(JPH_USE_SSE)
	__m128d zero = _mm_setzero_pd();
	return DVec4({ zero, zero });
#elif defined(JPH_USE_NEON)
	float64x2_t zero = vdupq_n_f64(0.0);
	return DVec4({ zero, zero });
#elif defined(JPH_USE_RVV)
	DVec4 vec;
	const vfloat64m2_t v = __riscv_vfmv_v_f_f64m2(0.0, 3);
	__riscv_vse64_v_f64m2(vec.mF64, v, 3);
	return vec;
#else
	return DVec4(0, 0, 0);
#endif
}

DVec4 DVec4::sReplicate(double inV)
{
#if defined(JPH_USE_AVX)
	return _mm256_set1_pd(inV);
#elif defined(JPH_USE_SSE)
	__m128d value = _mm_set1_pd(inV);
	return DVec4({ value, value });
#elif defined(JPH_USE_NEON)
	float64x2_t value = vdupq_n_f64(inV);
	return DVec4({ value, value });
#elif defined(JPH_USE_RVV)
	DVec4 vec;
	const vfloat64m2_t v = __riscv_vfmv_v_f_f64m2(inV, 3);
	__riscv_vse64_v_f64m2(vec.mF64, v, 3);
	return vec;
#else
	return DVec4(inV, inV, inV);
#endif
}

DVec4 DVec4::sOne()
{
	return sReplicate(1.0);
}

DVec4 DVec4::sNaN()
{
	return sReplicate(numeric_limits<double>::quiet_NaN());
}

DVec4 DVec4::sLoadDouble3Unsafe(const Double3 &inV)
{
#if defined(JPH_USE_AVX)
	Type v = _mm256_loadu_pd(&inV.x);
#elif defined(JPH_USE_SSE)
	Type v;
	v.mLow = _mm_loadu_pd(&inV.x);
	v.mHigh = _mm_set1_pd(inV.z);
#elif defined(JPH_USE_NEON)
	Type v = vld1q_f64_x2(&inV.x);
#elif defined(JPH_USE_RVV)
	Type v;
	const vfloat64m2_t vec = __riscv_vle64_v_f64m2(&inV.x, 3);
	__riscv_vse64_v_f64m2(v.mData, vec, 3);
#else
	Type v = { inV.x, inV.y, inV.z };
#endif
	return sFixW(v);
}

void DVec4::StoreDouble3(Double3 *outV) const
{
	outV->x = mF64[0];
	outV->y = mF64[1];
	outV->z = mF64[2];
}

DVec4::DVec4(const Double4 &inV)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_loadu_pd(&inV.x);
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_loadu_pd(&inV.x);
	mValue.mHigh = _mm_loadu_pd(&inV.z);
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vld1q_f64(&inV.x);
	mValue.val[1] = vld1q_f64(&inV.z);
#else
	mF64[0] = inV.x;
	mF64[1] = inV.y;
	mF64[2] = inV.z;
	mF64[3] = inV.w;
#endif
}

void DVec4::StoreDouble4(Double4 *outV) const
{
	outV->x = mF64[0];
	outV->y = mF64[1];
	outV->z = mF64[2];
	outV->w = mF64[3];
}

DVec4::operator Vec3() const
{
#if defined(JPH_USE_AVX)
	return _mm256_cvtpd_ps(mValue);
#elif defined(JPH_USE_SSE)
	__m128 low = _mm_cvtpd_ps(mValue.mLow);
	__m128 high = _mm_cvtpd_ps(mValue.mHigh);
	return _mm_shuffle_ps(low, high, _MM_SHUFFLE(1, 0, 1, 0));
#elif defined(JPH_USE_NEON)
	return vcvt_high_f32_f64(vcvtx_f32_f64(mValue.val[0]), mValue.val[1]);
#elif defined(JPH_USE_RVV)
	Vec3 v;
	const vfloat64m2_t src = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat32m1_t narrowed = __riscv_vfncvt_f_f_w_f32m1(src, 3);
	__riscv_vse32_v_f32m1(v.mF32, narrowed, 3);
	return v;
#else
	return Vec3((float)GetX(), (float)GetY(), (float)GetZ());
#endif
}

DVec4::operator Vec4() const
{
#if defined(JPH_USE_AVX)
	return _mm256_cvtpd_ps(mValue);
#elif defined(JPH_USE_SSE)
	__m128 low = _mm_cvtpd_ps(mValue.mLow);
	__m128 high = _mm_cvtpd_ps(mValue.mHigh);
	return _mm_movelh_ps(low, high);
#elif defined(JPH_USE_NEON)
	return vcvt_high_f32_f64(vcvt_f32_f64(mValue.val[0]), mValue.val[1]);
#elif defined(JPH_USE_RVV)
	Vec4 v;
	const vfloat64m2_t src = __riscv_vle64_v_f64m2(mF64, 4);
	const vfloat32m1_t narrowed = __riscv_vfncvt_f_f_w_f32m1(src, 4);
	__riscv_vse32_v_f32m1(v.mF32, narrowed, 4);
	return v;
#else
	return Vec4((float)GetX(), (float)GetY(), (float)GetZ(), (float)GetW());
#endif
}

DVec4 DVec4::sMin(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_min_pd(inV1.mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_min_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_min_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vminq_f64(inV1.mValue.val[0], inV2.mValue.val[0]), vminq_f64(inV1.mValue.val[1], inV2.mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV1.mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t min = __riscv_vfmin_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, min, 3);
	return res;
#else
	return DVec4(min(inV1.mF64[0], inV2.mF64[0]),
				 min(inV1.mF64[1], inV2.mF64[1]),
				 min(inV1.mF64[2], inV2.mF64[2]));
#endif
}

DVec4 DVec4::sMax(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_max_pd(inV1.mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_max_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_max_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vmaxq_f64(inV1.mValue.val[0], inV2.mValue.val[0]), vmaxq_f64(inV1.mValue.val[1], inV2.mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV1.mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t max = __riscv_vfmax_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, max, 3);
	return res;
#else
	return DVec4(max(inV1.mF64[0], inV2.mF64[0]),
				 max(inV1.mF64[1], inV2.mF64[1]),
				 max(inV1.mF64[2], inV2.mF64[2]));
#endif
}

DVec4 DVec4::sClamp(DVec4Arg inV, DVec4Arg inMin, DVec4Arg inMax)
{
	return sMax(sMin(inV, inMax), inMin);
}

DVec4 DVec4::sEquals(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_cmp_pd(inV1.mValue, inV2.mValue, _CMP_EQ_OQ);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_cmpeq_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_cmpeq_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vreinterpretq_f64_u64(vceqq_f64(inV1.mValue.val[0], inV2.mValue.val[0])), vreinterpretq_f64_u64(vceqq_f64(inV1.mValue.val[1], inV2.mValue.val[1])) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV1.mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vbool32_t mask = __riscv_vmfeq_vv_f64m2_b32(v1, v2, 3);
	const vfloat64m2_t zeros = __riscv_vfmv_v_f_f64m2(cFalse, 3);
	const vfloat64m2_t merged = __riscv_vfmerge_vfm_f64m2(zeros, cTrue, mask, 3);
	__riscv_vse64_v_f64m2(res.mF64, merged, 3);
	return res;
#else
	return DVec4(inV1.mF64[0] == inV2.mF64[0]? cTrue : cFalse,
				 inV1.mF64[1] == inV2.mF64[1]? cTrue : cFalse,
				 inV1.mF64[2] == inV2.mF64[2]? cTrue : cFalse);
#endif
}

DVec4 DVec4::sLess(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_cmp_pd(inV1.mValue, inV2.mValue, _CMP_LT_OQ);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_cmplt_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_cmplt_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vreinterpretq_f64_u64(vcltq_f64(inV1.mValue.val[0], inV2.mValue.val[0])), vreinterpretq_f64_u64(vcltq_f64(inV1.mValue.val[1], inV2.mValue.val[1])) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV1.mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vbool32_t mask = __riscv_vmflt_vv_f64m2_b32(v1, v2, 3);
	const vfloat64m2_t zeros = __riscv_vfmv_v_f_f64m2(cFalse, 3);
	const vfloat64m2_t merged = __riscv_vfmerge_vfm_f64m2(zeros, cTrue, mask, 3);
	__riscv_vse64_v_f64m2(res.mF64, merged, 3);
	return res;
#else
	return DVec4(inV1.mF64[0] < inV2.mF64[0]? cTrue : cFalse,
				 inV1.mF64[1] < inV2.mF64[1]? cTrue : cFalse,
				 inV1.mF64[2] < inV2.mF64[2]? cTrue : cFalse);
#endif
}

DVec4 DVec4::sLessOrEqual(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_cmp_pd(inV1.mValue, inV2.mValue, _CMP_LE_OQ);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_cmple_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_cmple_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vreinterpretq_f64_u64(vcleq_f64(inV1.mValue.val[0], inV2.mValue.val[0])), vreinterpretq_f64_u64(vcleq_f64(inV1.mValue.val[1], inV2.mValue.val[1])) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV1.mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vbool32_t mask = __riscv_vmfle_vv_f64m2_b32(v1, v2, 3);
	const vfloat64m2_t zeros = __riscv_vfmv_v_f_f64m2(cFalse, 3);
	const vfloat64m2_t merged = __riscv_vfmerge_vfm_f64m2(zeros, cTrue, mask, 3);
	__riscv_vse64_v_f64m2(res.mF64, merged, 3);
	return res;
#else
	return DVec4(inV1.mF64[0] <= inV2.mF64[0]? cTrue : cFalse,
				 inV1.mF64[1] <= inV2.mF64[1]? cTrue : cFalse,
				 inV1.mF64[2] <= inV2.mF64[2]? cTrue : cFalse);
#endif
}

DVec4 DVec4::sGreater(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_cmp_pd(inV1.mValue, inV2.mValue, _CMP_GT_OQ);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_cmpgt_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_cmpgt_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vreinterpretq_f64_u64(vcgtq_f64(inV1.mValue.val[0], inV2.mValue.val[0])), vreinterpretq_f64_u64(vcgtq_f64(inV1.mValue.val[1], inV2.mValue.val[1])) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV1.mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vbool32_t mask = __riscv_vmfgt_vv_f64m2_b32(v1, v2, 3);
	const vfloat64m2_t zeros = __riscv_vfmv_v_f_f64m2(cFalse, 3);
	const vfloat64m2_t merged = __riscv_vfmerge_vfm_f64m2(zeros, cTrue, mask, 3);
	__riscv_vse64_v_f64m2(res.mF64, merged, 3);
	return res;
#else
	return DVec4(inV1.mF64[0] > inV2.mF64[0]? cTrue : cFalse,
				 inV1.mF64[1] > inV2.mF64[1]? cTrue : cFalse,
				 inV1.mF64[2] > inV2.mF64[2]? cTrue : cFalse);
#endif
}

DVec4 DVec4::sGreaterOrEqual(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_cmp_pd(inV1.mValue, inV2.mValue, _CMP_GE_OQ);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_cmpge_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_cmpge_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vreinterpretq_f64_u64(vcgeq_f64(inV1.mValue.val[0], inV2.mValue.val[0])), vreinterpretq_f64_u64(vcgeq_f64(inV1.mValue.val[1], inV2.mValue.val[1])) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV1.mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vbool32_t mask = __riscv_vmfge_vv_f64m2_b32(v1, v2, 3);
	const vfloat64m2_t zeros = __riscv_vfmv_v_f_f64m2(cFalse, 3);
	const vfloat64m2_t merged = __riscv_vfmerge_vfm_f64m2(zeros, cTrue, mask, 3);
	__riscv_vse64_v_f64m2(res.mF64, merged, 3);
	return res;
#else
	return DVec4(inV1.mF64[0] >= inV2.mF64[0]? cTrue : cFalse,
				 inV1.mF64[1] >= inV2.mF64[1]? cTrue : cFalse,
				 inV1.mF64[2] >= inV2.mF64[2]? cTrue : cFalse);
#endif
}

DVec4 DVec4::sFusedMultiplyAdd(DVec4Arg inMul1, DVec4Arg inMul2, DVec4Arg inAdd)
{
#if defined(JPH_USE_AVX)
	#ifdef JPH_USE_FMADD
		return _mm256_fmadd_pd(inMul1.mValue, inMul2.mValue, inAdd.mValue);
	#else
		return _mm256_add_pd(_mm256_mul_pd(inMul1.mValue, inMul2.mValue), inAdd.mValue);
	#endif
#elif defined(JPH_USE_NEON)
	return DVec4({ vmlaq_f64(inAdd.mValue.val[0], inMul1.mValue.val[0], inMul2.mValue.val[0]), vmlaq_f64(inAdd.mValue.val[1], inMul1.mValue.val[1], inMul2.mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inMul1.mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inMul2.mF64, 3);
	const vfloat64m2_t rvv_add = __riscv_vle64_v_f64m2(inAdd.mF64, 3);
	const vfloat64m2_t fmadd = __riscv_vfmacc_vv_f64m2(rvv_add, v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, fmadd, 3);
	return res;
#else
	return inMul1 * inMul2 + inAdd;
#endif
}

DVec4 DVec4::sSelect(DVec4Arg inNotSet, DVec4Arg inSet, DVec4Arg inControl)
{
#if defined(JPH_USE_AVX)
	return _mm256_blendv_pd(inNotSet.mValue, inSet.mValue, inControl.mValue);
#elif defined(JPH_USE_SSE4_1)
	Type v = { _mm_blendv_pd(inNotSet.mValue.mLow, inSet.mValue.mLow, inControl.mValue.mLow), _mm_blendv_pd(inNotSet.mValue.mHigh, inSet.mValue.mHigh, inControl.mValue.mHigh) };
	return sFixW(v);
#elif defined(JPH_USE_NEON)
	Type v = { vbslq_f64(vreinterpretq_u64_s64(vshrq_n_s64(vreinterpretq_s64_f64(inControl.mValue.val[0]), 63)), inSet.mValue.val[0], inNotSet.mValue.val[0]),
			   vbslq_f64(vreinterpretq_u64_s64(vshrq_n_s64(vreinterpretq_s64_f64(inControl.mValue.val[1]), 63)), inSet.mValue.val[1], inNotSet.mValue.val[1]) };
	return sFixW(v);
#elif defined(JPH_USE_RVV)
	DVec4 masked;
	const vfloat64m2_t control_double = __riscv_vle64_v_f64m2(inControl.mF64, 3);
	const vfloat64m2_t not_set = __riscv_vle64_v_f64m2(inNotSet.mF64, 3);
	const vfloat64m2_t set = __riscv_vle64_v_f64m2(inSet.mF64, 3);
	const vuint64m2_t control = __riscv_vreinterpret_v_f64m2_u64m2(control_double);

	// Generate RVV bool mask from UVec4
	const uint64 sign_bit_mask = 0x8000000000000000u;
	const vuint64m2_t r = __riscv_vand_vx_u64m2(control, sign_bit_mask, 3);
	const vbool32_t rvv_mask = __riscv_vmsne_vx_u64m2_b32(r, 0x0, 3);
	const vfloat64m2_t merged = __riscv_vmerge_vvm_f64m2(not_set, set, rvv_mask, 3);
	__riscv_vse64_v_f64m2(masked.mF64, merged, 3);
	return masked;
#else
	DVec4 result;
	for (int i = 0; i < 3; i++)
		result.mF64[i] = (BitCast<uint64>(inControl.mF64[i]) & (uint64(1) << 63))? inSet.mF64[i] : inNotSet.mF64[i];
#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
	result.mF64[3] = result.mF64[2];
#endif // JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
	return result;
#endif
}

DVec4 DVec4::sOr(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_or_pd(inV1.mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_or_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_or_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vreinterpretq_f64_u64(vorrq_u64(vreinterpretq_u64_f64(inV1.mValue.val[0]), vreinterpretq_u64_f64(inV2.mValue.val[0]))),
				   vreinterpretq_f64_u64(vorrq_u64(vreinterpretq_u64_f64(inV1.mValue.val[1]), vreinterpretq_u64_f64(inV2.mValue.val[1]))) });
#elif defined(JPH_USE_RVV)
	DVec4 or_result;
	const vuint64m2_t v1 = __riscv_vle64_v_u64m2(reinterpret_cast<const uint64 *>(inV1.mF64), 3);
	const vuint64m2_t v2 = __riscv_vle64_v_u64m2(reinterpret_cast<const uint64 *>(inV2.mF64), 3);
	const vuint64m2_t res = __riscv_vor_vv_u64m2(v1, v2, 3);
	__riscv_vse64_v_u64m2(reinterpret_cast<uint64 *>(or_result.mF64), res, 3);
	return or_result;
#else
	return DVec4(BitCast<double>(BitCast<uint64>(inV1.mF64[0]) | BitCast<uint64>(inV2.mF64[0])),
				 BitCast<double>(BitCast<uint64>(inV1.mF64[1]) | BitCast<uint64>(inV2.mF64[1])),
				 BitCast<double>(BitCast<uint64>(inV1.mF64[2]) | BitCast<uint64>(inV2.mF64[2])));
#endif
}

DVec4 DVec4::sXor(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_xor_pd(inV1.mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_xor_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_xor_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(inV1.mValue.val[0]), vreinterpretq_u64_f64(inV2.mValue.val[0]))),
				   vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(inV1.mValue.val[1]), vreinterpretq_u64_f64(inV2.mValue.val[1]))) });
#elif defined(JPH_USE_RVV)
	DVec4 xor_result;
	const vuint64m2_t v1 = __riscv_vle64_v_u64m2(reinterpret_cast<const uint64 *>(inV1.mF64), 3);
	const vuint64m2_t v2 = __riscv_vle64_v_u64m2(reinterpret_cast<const uint64 *>(inV2.mF64), 3);
	const vuint64m2_t res = __riscv_vxor_vv_u64m2(v1, v2, 3);
	__riscv_vse64_v_u64m2(reinterpret_cast<uint64 *>(xor_result.mF64), res, 3);
	return xor_result;
#else
	return DVec4(BitCast<double>(BitCast<uint64>(inV1.mF64[0]) ^ BitCast<uint64>(inV2.mF64[0])),
				 BitCast<double>(BitCast<uint64>(inV1.mF64[1]) ^ BitCast<uint64>(inV2.mF64[1])),
				 BitCast<double>(BitCast<uint64>(inV1.mF64[2]) ^ BitCast<uint64>(inV2.mF64[2])));
#endif
}

DVec4 DVec4::sAnd(DVec4Arg inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_and_pd(inV1.mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_and_pd(inV1.mValue.mLow, inV2.mValue.mLow), _mm_and_pd(inV1.mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(inV1.mValue.val[0]), vreinterpretq_u64_f64(inV2.mValue.val[0]))),
				   vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(inV1.mValue.val[1]), vreinterpretq_u64_f64(inV2.mValue.val[1]))) });
#elif defined(JPH_USE_RVV)
	DVec4 and_result;
	const vuint64m2_t v1 = __riscv_vle64_v_u64m2(reinterpret_cast<const uint64 *>(inV1.mF64), 3);
	const vuint64m2_t v2 = __riscv_vle64_v_u64m2(reinterpret_cast<const uint64 *>(inV2.mF64), 3);
	const vuint64m2_t res = __riscv_vand_vv_u64m2(v1, v2, 3);
	__riscv_vse64_v_u64m2(reinterpret_cast<uint64 *>(and_result.mF64), res, 3);
	return and_result;
#else
	return DVec4(BitCast<double>(BitCast<uint64>(inV1.mF64[0]) & BitCast<uint64>(inV2.mF64[0])),
				 BitCast<double>(BitCast<uint64>(inV1.mF64[1]) & BitCast<uint64>(inV2.mF64[1])),
				 BitCast<double>(BitCast<uint64>(inV1.mF64[2]) & BitCast<uint64>(inV2.mF64[2])));
#endif
}

int DVec4::GetTrues() const
{
#if defined(JPH_USE_AVX)
	return _mm256_movemask_pd(mValue) & 0x7;
#elif defined(JPH_USE_SSE)
	return (_mm_movemask_pd(mValue.mLow) + (_mm_movemask_pd(mValue.mHigh) << 2)) & 0x7;
#else
	return int((BitCast<uint64>(mF64[0]) >> 63) | ((BitCast<uint64>(mF64[1]) >> 63) << 1) | ((BitCast<uint64>(mF64[2]) >> 63) << 2));
#endif
}

bool DVec4::TestAnyTrue() const
{
	return GetTrues() != 0;
}

bool DVec4::TestAllTrue() const
{
	return GetTrues() == 0x7;
}

bool DVec4::operator == (DVec4Arg inV2) const
{
	return sEquals(*this, inV2).TestAllTrue();
}

bool DVec4::IsClose(DVec4Arg inV2, double inMaxDistSq) const
{
	return (inV2 - *this).LengthSq() <= inMaxDistSq;
}

bool DVec4::IsNearZero(double inMaxDistSq) const
{
	return LengthSq() <= inMaxDistSq;
}

DVec4 DVec4::operator * (DVec4Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_mul_pd(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_mul_pd(mValue.mLow, inV2.mValue.mLow), _mm_mul_pd(mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vmulq_f64(mValue.val[0], inV2.mValue.val[0]), vmulq_f64(mValue.val[1], inV2.mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t mul = __riscv_vfmul_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, mul, 3);
	return res;
#else
	return DVec4(mF64[0] * inV2.mF64[0], mF64[1] * inV2.mF64[1], mF64[2] * inV2.mF64[2]);
#endif
}

DVec4 DVec4::operator * (double inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_mul_pd(mValue, _mm256_set1_pd(inV2));
#elif defined(JPH_USE_SSE)
	__m128d v = _mm_set1_pd(inV2);
	return DVec4({ _mm_mul_pd(mValue.mLow, v), _mm_mul_pd(mValue.mHigh, v) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vmulq_n_f64(mValue.val[0], inV2), vmulq_n_f64(mValue.val[1], inV2) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t src = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t mul = __riscv_vfmul_vf_f64m2(src, inV2, 3);
	__riscv_vse64_v_f64m2(res.mF64, mul, 3);
	return res;
#else
	return DVec4(mF64[0] * inV2, mF64[1] * inV2, mF64[2] * inV2);
#endif
}

DVec4 operator * (double inV1, DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_mul_pd(_mm256_set1_pd(inV1), inV2.mValue);
#elif defined(JPH_USE_SSE)
	__m128d v = _mm_set1_pd(inV1);
	return DVec4({ _mm_mul_pd(v, inV2.mValue.mLow), _mm_mul_pd(v, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vmulq_n_f64(inV2.mValue.val[0], inV1), vmulq_n_f64(inV2.mValue.val[1], inV1) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t mul = __riscv_vfmul_vf_f64m2(v1, inV1, 3);
	__riscv_vse64_v_f64m2(res.mF64, mul, 3);
	return res;
#else
	return DVec4(inV1 * inV2.mF64[0], inV1 * inV2.mF64[1], inV1 * inV2.mF64[2]);
#endif
}

DVec4 DVec4::operator / (double inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_div_pd(mValue, _mm256_set1_pd(inV2));
#elif defined(JPH_USE_SSE)
	__m128d v = _mm_set1_pd(inV2);
	return DVec4({ _mm_div_pd(mValue.mLow, v), _mm_div_pd(mValue.mHigh, v) });
#elif defined(JPH_USE_NEON)
	float64x2_t v = vdupq_n_f64(inV2);
	return DVec4({ vdivq_f64(mValue.val[0], v), vdivq_f64(mValue.val[1], v) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t src = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t div = __riscv_vfdiv_vf_f64m2(src, inV2, 3);
	__riscv_vse64_v_f64m2(res.mF64, div, 3);
	return res;
#else
	return DVec4(mF64[0] / inV2, mF64[1] / inV2, mF64[2] / inV2);
#endif
}

DVec4 &DVec4::operator *= (double inV2)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_mul_pd(mValue, _mm256_set1_pd(inV2));
#elif defined(JPH_USE_SSE)
	__m128d v = _mm_set1_pd(inV2);
	mValue.mLow = _mm_mul_pd(mValue.mLow, v);
	mValue.mHigh = _mm_mul_pd(mValue.mHigh, v);
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vmulq_n_f64(mValue.val[0], inV2);
	mValue.val[1] = vmulq_n_f64(mValue.val[1], inV2);
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t res = __riscv_vfmul_vf_f64m2(v1, inV2, 3);
	__riscv_vse64_v_f64m2(mF64, res, 3);
#else
	for (int i = 0; i < 3; ++i)
		mF64[i] *= inV2;
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
	return *this;
}

DVec4 &DVec4::operator *= (DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_mul_pd(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_mul_pd(mValue.mLow, inV2.mValue.mLow);
	mValue.mHigh = _mm_mul_pd(mValue.mHigh, inV2.mValue.mHigh);
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vmulq_f64(mValue.val[0], inV2.mValue.val[0]);
	mValue.val[1] = vmulq_f64(mValue.val[1], inV2.mValue.val[1]);
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t rvv_res = __riscv_vfmul_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(mF64, rvv_res, 3);
#else
	for (int i = 0; i < 3; ++i)
		mF64[i] *= inV2.mF64[i];
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
	return *this;
}

DVec4 &DVec4::operator /= (double inV2)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_div_pd(mValue, _mm256_set1_pd(inV2));
#elif defined(JPH_USE_SSE)
	__m128d v = _mm_set1_pd(inV2);
	mValue.mLow = _mm_div_pd(mValue.mLow, v);
	mValue.mHigh = _mm_div_pd(mValue.mHigh, v);
#elif defined(JPH_USE_NEON)
	float64x2_t v = vdupq_n_f64(inV2);
	mValue.val[0] = vdivq_f64(mValue.val[0], v);
	mValue.val[1] = vdivq_f64(mValue.val[1], v);
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t v = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t res = __riscv_vfdiv_vf_f64m2(v, inV2, 3);
	__riscv_vse64_v_f64m2(mF64, res, 3);
#else
	for (int i = 0; i < 3; ++i)
		mF64[i] /= inV2;
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
	return *this;
}

DVec4 DVec4::operator + (Vec3Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_add_pd(mValue, _mm256_cvtps_pd(inV2.mValue));
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_add_pd(mValue.mLow, _mm_cvtps_pd(inV2.mValue)), _mm_add_pd(mValue.mHigh, _mm_cvtps_pd(_mm_shuffle_ps(inV2.mValue, inV2.mValue, _MM_SHUFFLE(2, 2, 2, 2)))) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vaddq_f64(mValue.val[0], vcvt_f64_f32(vget_low_f32(inV2.mValue))), vaddq_f64(mValue.val[1], vcvt_high_f64_f32(inV2.mValue)) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat32m1_t v2_f32 = __riscv_vle32_v_f32m1(inV2.mF32, 3);
	const vfloat64m2_t v2 = __riscv_vfwcvt_f_f_v_f64m2(v2_f32, 3);
	const vfloat64m2_t rvv_add = __riscv_vfadd_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, rvv_add, 3);
	return res;
#else
	return DVec4(mF64[0] + inV2.mF32[0], mF64[1] + inV2.mF32[1], mF64[2] + inV2.mF32[2]);
#endif
}

DVec4 DVec4::operator + (DVec4Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_add_pd(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_add_pd(mValue.mLow, inV2.mValue.mLow), _mm_add_pd(mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vaddq_f64(mValue.val[0], inV2.mValue.val[0]), vaddq_f64(mValue.val[1], inV2.mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t rvv_add = __riscv_vfadd_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, rvv_add, 3);
	return res;
#else
	return DVec4(mF64[0] + inV2.mF64[0], mF64[1] + inV2.mF64[1], mF64[2] + inV2.mF64[2]);
#endif
}

DVec4 &DVec4::operator += (Vec3Arg inV2)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_add_pd(mValue, _mm256_cvtps_pd(inV2.mValue));
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_add_pd(mValue.mLow, _mm_cvtps_pd(inV2.mValue));
	mValue.mHigh = _mm_add_pd(mValue.mHigh, _mm_cvtps_pd(_mm_shuffle_ps(inV2.mValue, inV2.mValue, _MM_SHUFFLE(2, 2, 2, 2))));
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vaddq_f64(mValue.val[0], vcvt_f64_f32(vget_low_f32(inV2.mValue)));
	mValue.val[1] = vaddq_f64(mValue.val[1], vcvt_high_f64_f32(inV2.mValue));
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat32m1_t v2_f32 = __riscv_vle32_v_f32m1(inV2.mF32, 3);
	const vfloat64m2_t v2 = __riscv_vfwcvt_f_f_v_f64m2(v2_f32, 3);
	const vfloat64m2_t rvv_add = __riscv_vfadd_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(mF64, rvv_add, 3);
#else
	for (int i = 0; i < 3; ++i)
		mF64[i] += inV2.mF32[i];
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
	return *this;
}

DVec4 &DVec4::operator += (DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_add_pd(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_add_pd(mValue.mLow, inV2.mValue.mLow);
	mValue.mHigh = _mm_add_pd(mValue.mHigh, inV2.mValue.mHigh);
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vaddq_f64(mValue.val[0], inV2.mValue.val[0]);
	mValue.val[1] = vaddq_f64(mValue.val[1], inV2.mValue.val[1]);
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t rvv_add = __riscv_vfadd_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(mF64, rvv_add, 3);
#else
	for (int i = 0; i < 3; ++i)
		mF64[i] += inV2.mF64[i];
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
	return *this;
}

DVec4 DVec4::operator - () const
{
#if defined(JPH_USE_AVX)
	return _mm256_sub_pd(_mm256_setzero_pd(), mValue);
#elif defined(JPH_USE_SSE)
	__m128d zero = _mm_setzero_pd();
	return DVec4({ _mm_sub_pd(zero, mValue.mLow), _mm_sub_pd(zero, mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	#ifdef JPH_CROSS_PLATFORM_DETERMINISTIC
		float64x2_t zero = vdupq_n_f64(0);
		return DVec4({ vsubq_f64(zero, mValue.val[0]), vsubq_f64(zero, mValue.val[1]) });
	#else
		return DVec4({ vnegq_f64(mValue.val[0]), vnegq_f64(mValue.val[1]) });
	#endif
#elif defined(JPH_USE_RVV)
	#ifdef JPH_CROSS_PLATFORM_DETERMINISTIC
		DVec4 res;
		const vfloat64m2_t rvv_zero = __riscv_vfmv_v_f_f64m2(0.0, 3);
		const vfloat64m2_t v = __riscv_vle64_v_f64m2(mF64, 3);
		const vfloat64m2_t rvv_neg = __riscv_vfsub_vv_f64m2(rvv_zero, v, 3);
		__riscv_vse64_v_f64m2(res.mF64, rvv_neg, 3);
		return res;
	#else
		DVec4 res;
		const vfloat64m2_t v = __riscv_vle64_v_f64m2(mF64, 3);
		const vfloat64m2_t rvv_neg = __riscv_vfsgnjn_vv_f64m2(v, v, 3);
		__riscv_vse64_v_f64m2(res.mF64, rvv_neg, 3);
		return res;
	#endif
#else
	#ifdef JPH_CROSS_PLATFORM_DETERMINISTIC
		return DVec4(0.0 - mF64[0], 0.0 - mF64[1], 0.0 - mF64[2]);
	#else
		return DVec4(-mF64[0], -mF64[1], -mF64[2]);
	#endif
#endif
}

DVec4 DVec4::operator - (Vec3Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_sub_pd(mValue, _mm256_cvtps_pd(inV2.mValue));
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_sub_pd(mValue.mLow, _mm_cvtps_pd(inV2.mValue)), _mm_sub_pd(mValue.mHigh, _mm_cvtps_pd(_mm_shuffle_ps(inV2.mValue, inV2.mValue, _MM_SHUFFLE(2, 2, 2, 2)))) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vsubq_f64(mValue.val[0], vcvt_f64_f32(vget_low_f32(inV2.mValue))), vsubq_f64(mValue.val[1], vcvt_high_f64_f32(inV2.mValue)) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat32m1_t v2_f32 = __riscv_vle32_v_f32m1(inV2.mF32, 3);
	const vfloat64m2_t v2 = __riscv_vfwcvt_f_f_v_f64m2(v2_f32, 3);
	const vfloat64m2_t rvv_sub = __riscv_vfsub_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, rvv_sub, 3);
	return res;
#else
	return DVec4(mF64[0] - inV2.mF32[0], mF64[1] - inV2.mF32[1], mF64[2] - inV2.mF32[2]);
#endif
}

DVec4 DVec4::operator - (DVec4Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_sub_pd(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_sub_pd(mValue.mLow, inV2.mValue.mLow), _mm_sub_pd(mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vsubq_f64(mValue.val[0], inV2.mValue.val[0]), vsubq_f64(mValue.val[1], inV2.mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t rvv_sub = __riscv_vfsub_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, rvv_sub, 3);
	return res;
#else
	return DVec4(mF64[0] - inV2.mF64[0], mF64[1] - inV2.mF64[1], mF64[2] - inV2.mF64[2]);
#endif
}

DVec4 &DVec4::operator -= (Vec3Arg inV2)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_sub_pd(mValue, _mm256_cvtps_pd(inV2.mValue));
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_sub_pd(mValue.mLow, _mm_cvtps_pd(inV2.mValue));
	mValue.mHigh = _mm_sub_pd(mValue.mHigh, _mm_cvtps_pd(_mm_shuffle_ps(inV2.mValue, inV2.mValue, _MM_SHUFFLE(2, 2, 2, 2))));
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vsubq_f64(mValue.val[0], vcvt_f64_f32(vget_low_f32(inV2.mValue)));
	mValue.val[1] = vsubq_f64(mValue.val[1], vcvt_high_f64_f32(inV2.mValue));
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat32m1_t v2_f32 = __riscv_vle32_v_f32m1(inV2.mF32, 3);
	const vfloat64m2_t v2 = __riscv_vfwcvt_f_f_v_f64m2(v2_f32, 3);
	const vfloat64m2_t rvv_sub = __riscv_vfsub_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(mF64, rvv_sub, 3);
#else
	for (int i = 0; i < 3; ++i)
		mF64[i] -= inV2.mF32[i];
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
	return *this;
}

DVec4 &DVec4::operator -= (DVec4Arg inV2)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_sub_pd(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_sub_pd(mValue.mLow, inV2.mValue.mLow);
	mValue.mHigh = _mm_sub_pd(mValue.mHigh, inV2.mValue.mHigh);
#elif defined(JPH_USE_NEON)
	mValue.val[0] = vsubq_f64(mValue.val[0], inV2.mValue.val[0]);
	mValue.val[1] = vsubq_f64(mValue.val[1], inV2.mValue.val[1]);
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t rvv_sub = __riscv_vfsub_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(mF64, rvv_sub, 3);
#else
	for (int i = 0; i < 3; ++i)
		mF64[i] -= inV2.mF64[i];
	#ifdef JPH_FLOATING_POINT_EXCEPTIONS_ENABLED
		mF64[3] = mF64[2];
	#endif
#endif
	return *this;
}

DVec4 DVec4::operator / (DVec4Arg inV2) const
{
	inV2.CheckW();
#if defined(JPH_USE_AVX)
	return _mm256_div_pd(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_div_pd(mValue.mLow, inV2.mValue.mLow), _mm_div_pd(mValue.mHigh, inV2.mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vdivq_f64(mValue.val[0], inV2.mValue.val[0]), vdivq_f64(mValue.val[1], inV2.mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t rvv_div = __riscv_vfdiv_vv_f64m2(v1, v2, 3);
	__riscv_vse64_v_f64m2(res.mF64, rvv_div, 3);
	return res;
#else
	return DVec4(mF64[0] / inV2.mF64[0], mF64[1] / inV2.mF64[1], mF64[2] / inV2.mF64[2]);
#endif
}

DVec4 DVec4::Abs() const
{
#if defined(JPH_USE_AVX512)
	return _mm256_range_pd(mValue, mValue, 0b1000);
#elif defined(JPH_USE_AVX)
	return _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), mValue), mValue);
#elif defined(JPH_USE_SSE)
	__m128d zero = _mm_setzero_pd();
	return DVec4({ _mm_max_pd(_mm_sub_pd(zero, mValue.mLow), mValue.mLow), _mm_max_pd(_mm_sub_pd(zero, mValue.mHigh), mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vabsq_f64(mValue.val[0]), vabsq_f64(mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t rvv_abs = __riscv_vfsgnj_vf_f64m2(v, 1.0, 3);
	__riscv_vse64_v_f64m2(res.mF64, rvv_abs, 3);
	return res;
#else
	return DVec4(abs(mF64[0]), abs(mF64[1]), abs(mF64[2]));
#endif
}

DVec4 DVec4::Reciprocal() const
{
	return sOne() / mValue;
}

DVec4 DVec4::Cross(DVec4Arg inV2) const
{
#if defined(JPH_USE_AVX2)
	__m256d t1 = _mm256_permute4x64_pd(inV2.mValue, _MM_SHUFFLE(0, 0, 2, 1)); // Assure Z and W are the same
	t1 = _mm256_mul_pd(t1, mValue);
	__m256d t2 = _mm256_permute4x64_pd(mValue, _MM_SHUFFLE(0, 0, 2, 1)); // Assure Z and W are the same
	t2 = _mm256_mul_pd(t2, inV2.mValue);
	__m256d t3 = _mm256_sub_pd(t1, t2);
	return _mm256_permute4x64_pd(t3, _MM_SHUFFLE(0, 0, 2, 1)); // Assure Z and W are the same
#elif defined(JPH_USE_RVV)
	const uint64 indices[3] = { 1, 2, 0 };
	const vuint64m2_t gather_indices = __riscv_vle64_v_u64m2(indices, 3);
	const vfloat64m2_t v0 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	vfloat64m2_t t0 = __riscv_vrgather_vv_f64m2(v1, gather_indices, 3);
	t0 =  __riscv_vfmul_vv_f64m2(t0, v0, 3);
	vfloat64m2_t t1 = __riscv_vrgather_vv_f64m2(v0, gather_indices, 3);
	t1 =  __riscv_vfmul_vv_f64m2(t1, v1, 3);
	const vfloat64m2_t sub = __riscv_vfsub_vv_f64m2(t0, t1, 3);
	const vfloat64m2_t cross = __riscv_vrgather_vv_f64m2(sub, gather_indices, 3);

	DVec4 cross_result;
	__riscv_vse64_v_f64m2(cross_result.mF64, cross, 3);
	return cross_result;
#else
	return DVec4(mF64[1] * inV2.mF64[2] - mF64[2] * inV2.mF64[1],
				 mF64[2] * inV2.mF64[0] - mF64[0] * inV2.mF64[2],
				 mF64[0] * inV2.mF64[1] - mF64[1] * inV2.mF64[0]);
#endif
}

DVec4 DVec4::sCross(DVec4Arg inA, DVec4Arg inB, DVec4Arg inC)
{
	// 4D cross product: given 3 vectors A, B, C in R^4, returns a vector perpendicular to all three.
	// Computed as cofactor expansion of:
	// | e1   e2   e3   e4  |
	// | a.x  a.y  a.z  a.w |
	// | b.x  b.y  b.z  b.w |
	// | c.x  c.y  c.z  c.w |
	double ax = inA.mF64[0], ay = inA.mF64[1], az = inA.mF64[2], aw = inA.mF64[3];
	double bx = inB.mF64[0], by = inB.mF64[1], bz = inB.mF64[2], bw = inB.mF64[3];
	double cx = inC.mF64[0], cy = inC.mF64[1], cz = inC.mF64[2], cw = inC.mF64[3];

	double x =  ay * (bz * cw - bw * cz) - az * (by * cw - bw * cy) + aw * (by * cz - bz * cy);
	double y = -(ax * (bz * cw - bw * cz) - az * (bx * cw - bw * cx) + aw * (bx * cz - bz * cx));
	double z =  ax * (by * cw - bw * cy) - ay * (bx * cw - bw * cx) + aw * (bx * cy - by * cx);
	double w = -(ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx));

	return DVec4(x, y, z, w);
}

double DVec4::Dot(DVec4Arg inV2) const
{
#if defined(JPH_USE_AVX)
	__m256d mul = _mm256_mul_pd(mValue, inV2.mValue);
	__m128d xy = _mm256_castpd256_pd128(mul);
	__m128d yx = _mm_shuffle_pd(xy, xy, 1);
	__m128d sum = _mm_add_pd(xy, yx);
	__m128d zw = _mm256_extractf128_pd(mul, 1);
	sum = _mm_add_pd(sum, zw);
	return _mm_cvtsd_f64(sum);
#elif defined(JPH_USE_SSE)
	__m128d xy = _mm_mul_pd(mValue.mLow, inV2.mValue.mLow);
	__m128d yx = _mm_shuffle_pd(xy, xy, 1);
	__m128d sum = _mm_add_pd(xy, yx);
	__m128d z = _mm_mul_sd(mValue.mHigh, inV2.mValue.mHigh);
	sum = _mm_add_pd(sum, z);
	return _mm_cvtsd_f64(sum);
#elif defined(JPH_USE_NEON)
	float64x2_t mul_low = vmulq_f64(mValue.val[0], inV2.mValue.val[0]);
	float64x2_t mul_high = vmulq_f64(mValue.val[1], inV2.mValue.val[1]);
	return vaddvq_f64(mul_low) + vgetq_lane_f64(mul_high, 0);
#elif defined(JPH_USE_RVV)
	const vfloat64m1_t zeros = __riscv_vfmv_v_f_f64m1(0.0, 3);
	const vfloat64m2_t v1 = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t v2 = __riscv_vle64_v_f64m2(inV2.mF64, 3);
	const vfloat64m2_t mul = __riscv_vfmul_vv_f64m2(v1, v2, 3);
	const vfloat64m1_t sum = __riscv_vfredosum_vs_f64m2_f64m1(mul, zeros, 3);
	return __riscv_vfmv_f_s_f64m1_f64(sum);
#else
	double dot = 0.0;
	for (int i = 0; i < 3; i++)
		dot += mF64[i] * inV2.mF64[i];
	return dot;
#endif
}

double DVec4::LengthSq() const
{
	return Dot(*this);
}

DVec4 DVec4::Sqrt() const
{
#if defined(JPH_USE_AVX)
	return _mm256_sqrt_pd(mValue);
#elif defined(JPH_USE_SSE)
	return DVec4({ _mm_sqrt_pd(mValue.mLow), _mm_sqrt_pd(mValue.mHigh) });
#elif defined(JPH_USE_NEON)
	return DVec4({ vsqrtq_f64(mValue.val[0]), vsqrtq_f64(mValue.val[1]) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t v = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t rvv_sqrt = __riscv_vfsqrt_v_f64m2(v, 3);
	__riscv_vse64_v_f64m2(res.mF64, rvv_sqrt, 3);
	return res;
#else
	return DVec4(sqrt(mF64[0]), sqrt(mF64[1]), sqrt(mF64[2]));
#endif
}

double DVec4::Length() const
{
	return sqrt(Dot(*this));
}

DVec4 DVec4::Normalized() const
{
	return *this / Length();
}

bool DVec4::IsNormalized(double inTolerance) const
{
	return abs(LengthSq() - 1.0) <= inTolerance;
}

bool DVec4::IsNaN() const
{
#if defined(JPH_USE_AVX512)
	return (_mm256_fpclass_pd_mask(mValue, 0b10000001) & 0x7) != 0;
#elif defined(JPH_USE_AVX)
	return (_mm256_movemask_pd(_mm256_cmp_pd(mValue, mValue, _CMP_UNORD_Q)) & 0x7) != 0;
#elif defined(JPH_USE_SSE)
	return ((_mm_movemask_pd(_mm_cmpunord_pd(mValue.mLow, mValue.mLow)) + (_mm_movemask_pd(_mm_cmpunord_pd(mValue.mHigh, mValue.mHigh)) << 2)) & 0x7) != 0;
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t v = __riscv_vle64_v_f64m2(mF64, 3);
	const vbool32_t mask = __riscv_vmfeq_vv_f64m2_b32(v, v, 3);
	const uint32 eq = __riscv_vcpop_m_b32(mask, 3);
	return eq != 3;
#else
	return isnan(mF64[0]) || isnan(mF64[1]) || isnan(mF64[2]);
#endif
}

DVec4 DVec4::GetSign() const
{
#if defined(JPH_USE_AVX512)
	return _mm256_fixupimm_pd(mValue, mValue, _mm256_set1_epi32(0xA9A90A00), 0);
#elif defined(JPH_USE_AVX)
	__m256d minus_one = _mm256_set1_pd(-1.0);
	__m256d one = _mm256_set1_pd(1.0);
	return _mm256_or_pd(_mm256_and_pd(mValue, minus_one), one);
#elif defined(JPH_USE_SSE)
	__m128d minus_one = _mm_set1_pd(-1.0);
	__m128d one = _mm_set1_pd(1.0);
	return DVec4({ _mm_or_pd(_mm_and_pd(mValue.mLow, minus_one), one), _mm_or_pd(_mm_and_pd(mValue.mHigh, minus_one), one) });
#elif defined(JPH_USE_NEON)
	uint64x2_t minus_one = vreinterpretq_u64_f64(vdupq_n_f64(-1.0f));
	uint64x2_t one = vreinterpretq_u64_f64(vdupq_n_f64(1.0f));
	return DVec4({ vreinterpretq_f64_u64(vorrq_u64(vandq_u64(vreinterpretq_u64_f64(mValue.val[0]), minus_one), one)),
				   vreinterpretq_f64_u64(vorrq_u64(vandq_u64(vreinterpretq_u64_f64(mValue.val[1]), minus_one), one)) });
#elif defined(JPH_USE_RVV)
	DVec4 res;
	const vfloat64m2_t rvv_in = __riscv_vle64_v_f64m2(mF64, 3);
	const vfloat64m2_t rvv_one = __riscv_vfmv_v_f_f64m2(1.0, 3);
	const vfloat64m2_t rvv_signs = __riscv_vfsgnj_vv_f64m2(rvv_one, rvv_in, 3);
	__riscv_vse64_v_f64m2(res.mF64, rvv_signs, 3);
	return res;
#else
	return DVec4(std::signbit(mF64[0])? -1.0 : 1.0,
				 std::signbit(mF64[1])? -1.0 : 1.0,
				 std::signbit(mF64[2])? -1.0 : 1.0);
#endif
}

DVec4 DVec4::PrepareRoundToZero() const
{
	// Float has 23 bit mantissa, double 52 bit mantissa => we lose 29 bits when converting from double to float
	constexpr uint64 cDoubleToFloatMantissaLoss = (1U << 29) - 1;

#if defined(JPH_USE_AVX)
	return _mm256_and_pd(mValue, _mm256_castsi256_pd(_mm256_set1_epi64x(int64_t(~cDoubleToFloatMantissaLoss))));
#elif defined(JPH_USE_SSE)
	__m128d mask = _mm_castsi128_pd(_mm_set1_epi64x(int64_t(~cDoubleToFloatMantissaLoss)));
	return DVec4({ _mm_and_pd(mValue.mLow, mask), _mm_and_pd(mValue.mHigh, mask) });
#elif defined(JPH_USE_NEON)
	uint64x2_t mask = vdupq_n_u64(~cDoubleToFloatMantissaLoss);
	return DVec4({ vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(mValue.val[0]), mask)),
				   vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(mValue.val[1]), mask)) });
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t dvec = __riscv_vle64_v_f64m2(mF64, 3);
	const vuint64m2_t dvec_u64 = __riscv_vreinterpret_v_f64m2_u64m2(dvec);
	const vuint64m2_t chopped = __riscv_vand_vx_u64m2(dvec_u64, ~cDoubleToFloatMantissaLoss, 3);
	const vfloat64m2_t chopped_f64 = __riscv_vreinterpret_v_u64m2_f64m2(chopped);

	DVec4 res;
	__riscv_vse64_v_f64m2(res.mF64, chopped_f64, 3);
	return res;
#else
	double x = BitCast<double>(BitCast<uint64>(mF64[0]) & ~cDoubleToFloatMantissaLoss);
	double y = BitCast<double>(BitCast<uint64>(mF64[1]) & ~cDoubleToFloatMantissaLoss);
	double z = BitCast<double>(BitCast<uint64>(mF64[2]) & ~cDoubleToFloatMantissaLoss);

	return DVec4(x, y, z);
#endif
}

DVec4 DVec4::PrepareRoundToInf() const
{
	// Float has 23 bit mantissa, double 52 bit mantissa => we lose 29 bits when converting from double to float
	constexpr uint64 cDoubleToFloatMantissaLoss = (1U << 29) - 1;

#if defined(JPH_USE_AVX512)
	__m256i mantissa_loss = _mm256_set1_epi64x(cDoubleToFloatMantissaLoss);
	__mmask8 is_zero = _mm256_testn_epi64_mask(_mm256_castpd_si256(mValue), mantissa_loss);
	__m256d value_or_mantissa_loss = _mm256_or_pd(mValue, _mm256_castsi256_pd(mantissa_loss));
	return _mm256_mask_blend_pd(is_zero, value_or_mantissa_loss, mValue);
#elif defined(JPH_USE_AVX)
	__m256i mantissa_loss = _mm256_set1_epi64x(cDoubleToFloatMantissaLoss);
	__m256d value_and_mantissa_loss = _mm256_and_pd(mValue, _mm256_castsi256_pd(mantissa_loss));
	__m256d is_zero = _mm256_cmp_pd(value_and_mantissa_loss, _mm256_setzero_pd(), _CMP_EQ_OQ);
	__m256d value_or_mantissa_loss = _mm256_or_pd(mValue, _mm256_castsi256_pd(mantissa_loss));
	return _mm256_blendv_pd(value_or_mantissa_loss, mValue, is_zero);
#elif defined(JPH_USE_SSE4_1)
	__m128i mantissa_loss = _mm_set1_epi64x(cDoubleToFloatMantissaLoss);
	__m128d zero = _mm_setzero_pd();
	__m128d value_and_mantissa_loss_low = _mm_and_pd(mValue.mLow, _mm_castsi128_pd(mantissa_loss));
	__m128d is_zero_low = _mm_cmpeq_pd(value_and_mantissa_loss_low, zero);
	__m128d value_or_mantissa_loss_low = _mm_or_pd(mValue.mLow, _mm_castsi128_pd(mantissa_loss));
	__m128d value_and_mantissa_loss_high = _mm_and_pd(mValue.mHigh, _mm_castsi128_pd(mantissa_loss));
	__m128d is_zero_high = _mm_cmpeq_pd(value_and_mantissa_loss_high, zero);
	__m128d value_or_mantissa_loss_high = _mm_or_pd(mValue.mHigh, _mm_castsi128_pd(mantissa_loss));
	return DVec4({ _mm_blendv_pd(value_or_mantissa_loss_low, mValue.mLow, is_zero_low), _mm_blendv_pd(value_or_mantissa_loss_high, mValue.mHigh, is_zero_high) });
#elif defined(JPH_USE_NEON)
	uint64x2_t mantissa_loss = vdupq_n_u64(cDoubleToFloatMantissaLoss);
	float64x2_t zero = vdupq_n_f64(0.0);
	float64x2_t value_and_mantissa_loss_low = vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(mValue.val[0]), mantissa_loss));
	uint64x2_t is_zero_low = vceqq_f64(value_and_mantissa_loss_low, zero);
	float64x2_t value_or_mantissa_loss_low = vreinterpretq_f64_u64(vorrq_u64(vreinterpretq_u64_f64(mValue.val[0]), mantissa_loss));
	float64x2_t value_and_mantissa_loss_high = vreinterpretq_f64_u64(vandq_u64(vreinterpretq_u64_f64(mValue.val[1]), mantissa_loss));
	float64x2_t value_low = vbslq_f64(is_zero_low, mValue.val[0], value_or_mantissa_loss_low);
	uint64x2_t is_zero_high = vceqq_f64(value_and_mantissa_loss_high, zero);
	float64x2_t value_or_mantissa_loss_high = vreinterpretq_f64_u64(vorrq_u64(vreinterpretq_u64_f64(mValue.val[1]), mantissa_loss));
	float64x2_t value_high = vbslq_f64(is_zero_high, mValue.val[1], value_or_mantissa_loss_high);
	return DVec4({ value_low, value_high });
#elif defined(JPH_USE_RVV)
	const vfloat64m2_t dvec = __riscv_vle64_v_f64m2(mF64, 3);
	const vuint64m2_t dvec_u64 = __riscv_vreinterpret_v_f64m2_u64m2(dvec);
	const vuint64m2_t and_loss = __riscv_vand_vx_u64m2(dvec_u64, cDoubleToFloatMantissaLoss, 3);
	const vuint64m2_t or_loss = __riscv_vor_vx_u64m2(dvec_u64, cDoubleToFloatMantissaLoss, 3);
	const vbool32_t is_zero = __riscv_vmseq_vx_u64m2_b32(and_loss, 0x0, 3);
	const vuint64m2_t select = __riscv_vmerge_vvm_u64m2(or_loss, dvec_u64, is_zero, 3);
	const vfloat64m2_t select_f64 = __riscv_vreinterpret_v_u64m2_f64m2(select);

	DVec4 res;
	__riscv_vse64_v_f64m2(res.mF64, select_f64, 3);
	return res;
#else
	uint64 ux = BitCast<uint64>(mF64[0]);
	uint64 uy = BitCast<uint64>(mF64[1]);
	uint64 uz = BitCast<uint64>(mF64[2]);

	double x = BitCast<double>((ux & cDoubleToFloatMantissaLoss) == 0? ux : (ux | cDoubleToFloatMantissaLoss));
	double y = BitCast<double>((uy & cDoubleToFloatMantissaLoss) == 0? uy : (uy | cDoubleToFloatMantissaLoss));
	double z = BitCast<double>((uz & cDoubleToFloatMantissaLoss) == 0? uz : (uz | cDoubleToFloatMantissaLoss));

	return DVec4(x, y, z);
#endif
}

Vec3 DVec4::ToVec3RoundDown() const
{
	DVec4 to_zero = PrepareRoundToZero();
	DVec4 to_inf = PrepareRoundToInf();
	return Vec3(DVec4::sSelect(to_zero, to_inf, DVec4::sLess(*this, DVec4::sZero())));
}

Vec3 DVec4::ToVec3RoundUp() const
{
	DVec4 to_zero = PrepareRoundToZero();
	DVec4 to_inf = PrepareRoundToInf();
	return Vec3(DVec4::sSelect(to_inf, to_zero, DVec4::sLess(*this, DVec4::sZero())));
}

JPH_NAMESPACE_END
