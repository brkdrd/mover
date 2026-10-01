// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Math/Lane4.h>
#include <Jolt/Core/HashCombine.h>

// Create a std::hash/JPH::Hash for Vec8
JPH_MAKE_HASHABLE(JPH::Vec8, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7])

JPH_NAMESPACE_BEGIN

Vec8::Vec8(Lane4Arg inLow, Lane4Arg inHigh)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_set_m128(inHigh.mValue, inLow.mValue);
#elif defined(JPH_USE_SSE) || defined(JPH_USE_NEON)
	mValue.mLow = inLow.mValue;
	mValue.mHigh = inHigh.mValue;
#else
	for (int i = 0; i < 4; i++)
	{
		mF32[i] = inLow.mF32[i];
		mF32[i + 4] = inHigh.mF32[i];
	}
#endif
}

Vec8::Vec8(float in0, float in1, float in2, float in3,
		   float in4, float in5, float in6, float in7)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_set_ps(in7, in6, in5, in4, in3, in2, in1, in0);
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_set_ps(in3, in2, in1, in0);
	mValue.mHigh = _mm_set_ps(in7, in6, in5, in4);
#elif defined(JPH_USE_NEON)
	uint32x2_t low01 = vcreate_u32(static_cast<uint64>(BitCast<uint32>(in0)) | (static_cast<uint64>(BitCast<uint32>(in1)) << 32));
	uint32x2_t low23 = vcreate_u32(static_cast<uint64>(BitCast<uint32>(in2)) | (static_cast<uint64>(BitCast<uint32>(in3)) << 32));
	mValue.mLow = vreinterpretq_f32_u32(vcombine_u32(low01, low23));
	uint32x2_t hi01 = vcreate_u32(static_cast<uint64>(BitCast<uint32>(in4)) | (static_cast<uint64>(BitCast<uint32>(in5)) << 32));
	uint32x2_t hi23 = vcreate_u32(static_cast<uint64>(BitCast<uint32>(in6)) | (static_cast<uint64>(BitCast<uint32>(in7)) << 32));
	mValue.mHigh = vreinterpretq_f32_u32(vcombine_u32(hi01, hi23));
#else
	mF32[0] = in0; mF32[1] = in1; mF32[2] = in2; mF32[3] = in3;
	mF32[4] = in4; mF32[5] = in5; mF32[6] = in6; mF32[7] = in7;
#endif
}

Vec8::Vec8(const Float8 &inV)
{
#if defined(JPH_USE_AVX)
	mValue = _mm256_loadu_ps(inV.mValue);
#elif defined(JPH_USE_SSE)
	mValue.mLow = _mm_loadu_ps(inV.mValue);
	mValue.mHigh = _mm_loadu_ps(inV.mValue + 4);
#elif defined(JPH_USE_NEON)
	mValue.mLow = vld1q_f32(inV.mValue);
	mValue.mHigh = vld1q_f32(inV.mValue + 4);
#else
	for (int i = 0; i < 8; i++)
		mF32[i] = inV.mValue[i];
#endif
}

Vec8 Vec8::sZero()
{
#if defined(JPH_USE_AVX)
	return _mm256_setzero_ps();
#elif defined(JPH_USE_SSE)
	Vec8 v;
	v.mValue.mLow = _mm_setzero_ps();
	v.mValue.mHigh = _mm_setzero_ps();
	return v;
#elif defined(JPH_USE_NEON)
	Vec8 v;
	v.mValue.mLow = vdupq_n_f32(0);
	v.mValue.mHigh = vdupq_n_f32(0);
	return v;
#else
	Vec8 v;
	for (int i = 0; i < 8; i++)
		v.mF32[i] = 0.0f;
	return v;
#endif
}

Vec8 Vec8::sReplicate(float inV)
{
#if defined(JPH_USE_AVX)
	return _mm256_set1_ps(inV);
#elif defined(JPH_USE_SSE)
	Vec8 v;
	v.mValue.mLow = _mm_set1_ps(inV);
	v.mValue.mHigh = _mm_set1_ps(inV);
	return v;
#elif defined(JPH_USE_NEON)
	Vec8 v;
	v.mValue.mLow = vdupq_n_f32(inV);
	v.mValue.mHigh = vdupq_n_f32(inV);
	return v;
#else
	Vec8 v;
	for (int i = 0; i < 8; i++)
		v.mF32[i] = inV;
	return v;
#endif
}

Lane4 Vec8::GetLow() const
{
#if defined(JPH_USE_AVX)
	return _mm256_castps256_ps128(mValue);
#elif defined(JPH_USE_SSE) || defined(JPH_USE_NEON)
	return mValue.mLow;
#else
	return Lane4(mF32[0], mF32[1], mF32[2], mF32[3]);
#endif
}

Lane4 Vec8::GetHigh() const
{
#if defined(JPH_USE_AVX)
	return _mm256_extractf128_ps(mValue, 1);
#elif defined(JPH_USE_SSE) || defined(JPH_USE_NEON)
	return mValue.mHigh;
#else
	return Lane4(mF32[4], mF32[5], mF32[6], mF32[7]);
#endif
}

bool Vec8::operator == (Vec8Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return (_mm256_movemask_ps(_mm256_cmp_ps(mValue, inV2.mValue, _CMP_EQ_OQ)) & 0xff) == 0xff;
#elif defined(JPH_USE_SSE)
	return (_mm_movemask_ps(_mm_cmpeq_ps(mValue.mLow, inV2.mValue.mLow)) & 0xf) == 0xf
		&& (_mm_movemask_ps(_mm_cmpeq_ps(mValue.mHigh, inV2.mValue.mHigh)) & 0xf) == 0xf;
#elif defined(JPH_USE_NEON)
	uint32x4_t eq_low = vceqq_f32(mValue.mLow, inV2.mValue.mLow);
	uint32x4_t eq_high = vceqq_f32(mValue.mHigh, inV2.mValue.mHigh);
	return vminvq_u32(vandq_u32(eq_low, eq_high)) == 0xffffffffu;
#else
	for (int i = 0; i < 8; i++)
		if (mF32[i] != inV2.mF32[i])
			return false;
	return true;
#endif
}

bool Vec8::IsClose(Vec8Arg inV2, float inMaxDistSq) const
{
	Vec8 d = *this - inV2;
	return d.Dot(d) <= inMaxDistSq;
}

bool Vec8::IsNaN() const
{
#if defined(JPH_USE_AVX)
	return (_mm256_movemask_ps(_mm256_cmp_ps(mValue, mValue, _CMP_UNORD_Q)) & 0xff) != 0;
#elif defined(JPH_USE_SSE)
	return (_mm_movemask_ps(_mm_cmpunord_ps(mValue.mLow, mValue.mLow))
		  | _mm_movemask_ps(_mm_cmpunord_ps(mValue.mHigh, mValue.mHigh))) != 0;
#elif defined(JPH_USE_NEON)
	uint32x4_t low_nan = vmvnq_u32(vceqq_f32(mValue.mLow, mValue.mLow));
	uint32x4_t high_nan = vmvnq_u32(vceqq_f32(mValue.mHigh, mValue.mHigh));
	return vmaxvq_u32(vorrq_u32(low_nan, high_nan)) != 0;
#else
	for (int i = 0; i < 8; i++)
		if (isnan(mF32[i]))
			return true;
	return false;
#endif
}

Vec8 Vec8::operator + (Vec8Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_add_ps(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	Vec8 v;
	v.mValue.mLow = _mm_add_ps(mValue.mLow, inV2.mValue.mLow);
	v.mValue.mHigh = _mm_add_ps(mValue.mHigh, inV2.mValue.mHigh);
	return v;
#elif defined(JPH_USE_NEON)
	Vec8 v;
	v.mValue.mLow = vaddq_f32(mValue.mLow, inV2.mValue.mLow);
	v.mValue.mHigh = vaddq_f32(mValue.mHigh, inV2.mValue.mHigh);
	return v;
#else
	Vec8 v;
	for (int i = 0; i < 8; i++)
		v.mF32[i] = mF32[i] + inV2.mF32[i];
	return v;
#endif
}

Vec8 Vec8::operator - (Vec8Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_sub_ps(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	Vec8 v;
	v.mValue.mLow = _mm_sub_ps(mValue.mLow, inV2.mValue.mLow);
	v.mValue.mHigh = _mm_sub_ps(mValue.mHigh, inV2.mValue.mHigh);
	return v;
#elif defined(JPH_USE_NEON)
	Vec8 v;
	v.mValue.mLow = vsubq_f32(mValue.mLow, inV2.mValue.mLow);
	v.mValue.mHigh = vsubq_f32(mValue.mHigh, inV2.mValue.mHigh);
	return v;
#else
	Vec8 v;
	for (int i = 0; i < 8; i++)
		v.mF32[i] = mF32[i] - inV2.mF32[i];
	return v;
#endif
}

Vec8 Vec8::operator - () const
{
#if defined(JPH_USE_AVX)
	return _mm256_sub_ps(_mm256_setzero_ps(), mValue);
#elif defined(JPH_USE_SSE)
	Vec8 v;
	__m128 zero = _mm_setzero_ps();
	v.mValue.mLow = _mm_sub_ps(zero, mValue.mLow);
	v.mValue.mHigh = _mm_sub_ps(zero, mValue.mHigh);
	return v;
#elif defined(JPH_USE_NEON)
	Vec8 v;
	v.mValue.mLow = vnegq_f32(mValue.mLow);
	v.mValue.mHigh = vnegq_f32(mValue.mHigh);
	return v;
#else
	Vec8 v;
	for (int i = 0; i < 8; i++)
		v.mF32[i] = -mF32[i];
	return v;
#endif
}

Vec8 Vec8::sAnd(Vec8Arg inV1, Vec8Arg inV2)
{
#if defined(JPH_USE_AVX)
	return _mm256_and_ps(inV1.mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	Vec8 v;
	v.mValue.mLow = _mm_and_ps(inV1.mValue.mLow, inV2.mValue.mLow);
	v.mValue.mHigh = _mm_and_ps(inV1.mValue.mHigh, inV2.mValue.mHigh);
	return v;
#elif defined(JPH_USE_NEON)
	Vec8 v;
	v.mValue.mLow = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(inV1.mValue.mLow), vreinterpretq_u32_f32(inV2.mValue.mLow)));
	v.mValue.mHigh = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(inV1.mValue.mHigh), vreinterpretq_u32_f32(inV2.mValue.mHigh)));
	return v;
#else
	// No vector type to bit-and, so go through integers. memcpy, not a reinterpret_cast: punning a
	// float * through a uint32 * is strict-aliasing UB and gets miscompiled at -O2.
	Vec8 v;
	for (int i = 0; i < 8; i++)
	{
		uint32 a, b;
		memcpy(&a, &inV1.mF32[i], sizeof(uint32));
		memcpy(&b, &inV2.mF32[i], sizeof(uint32));
		uint32 r = a & b;
		memcpy(&v.mF32[i], &r, sizeof(uint32));
	}
	return v;
#endif
}

Vec8 Vec8::operator * (Vec8Arg inV2) const
{
#if defined(JPH_USE_AVX)
	return _mm256_mul_ps(mValue, inV2.mValue);
#elif defined(JPH_USE_SSE)
	Vec8 v;
	v.mValue.mLow = _mm_mul_ps(mValue.mLow, inV2.mValue.mLow);
	v.mValue.mHigh = _mm_mul_ps(mValue.mHigh, inV2.mValue.mHigh);
	return v;
#elif defined(JPH_USE_NEON)
	Vec8 v;
	v.mValue.mLow = vmulq_f32(mValue.mLow, inV2.mValue.mLow);
	v.mValue.mHigh = vmulq_f32(mValue.mHigh, inV2.mValue.mHigh);
	return v;
#else
	Vec8 v;
	for (int i = 0; i < 8; i++)
		v.mF32[i] = mF32[i] * inV2.mF32[i];
	return v;
#endif
}

Vec8 Vec8::operator * (float inV2) const
{
	return *this * sReplicate(inV2);
}

Vec8 operator * (float inV1, Vec8Arg inV2)
{
	return inV2 * inV1;
}

void Vec8::StoreFloat8(Float8 *outV) const
{
#if defined(JPH_USE_AVX)
	_mm256_storeu_ps(outV->mValue, mValue);
#elif defined(JPH_USE_SSE)
	_mm_storeu_ps(outV->mValue, mValue.mLow);
	_mm_storeu_ps(outV->mValue + 4, mValue.mHigh);
#elif defined(JPH_USE_NEON)
	vst1q_f32(outV->mValue, mValue.mLow);
	vst1q_f32(outV->mValue + 4, mValue.mHigh);
#else
	for (int i = 0; i < 8; i++)
		outV->mValue[i] = mF32[i];
#endif
}

float Vec8::Dot(Vec8Arg inV2) const
{
#if defined(JPH_USE_AVX)
	__m256 mul = _mm256_mul_ps(mValue, inV2.mValue);
	__m128 low = _mm256_castps256_ps128(mul);
	__m128 high = _mm256_extractf128_ps(mul, 1);
	__m128 sum4 = _mm_add_ps(low, high);
	__m128 shuf = _mm_movehdup_ps(sum4);
	__m128 sums = _mm_add_ps(sum4, shuf);
	shuf = _mm_movehl_ps(shuf, sums);
	sums = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
#elif defined(JPH_USE_SSE)
	__m128 mul_low = _mm_mul_ps(mValue.mLow, inV2.mValue.mLow);
	__m128 mul_high = _mm_mul_ps(mValue.mHigh, inV2.mValue.mHigh);
	__m128 sum4 = _mm_add_ps(mul_low, mul_high);
	__m128 shuf = _mm_movehdup_ps(sum4);
	__m128 sums = _mm_add_ps(sum4, shuf);
	shuf = _mm_movehl_ps(shuf, sums);
	sums = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
#elif defined(JPH_USE_NEON)
	float32x4_t mul_low = vmulq_f32(mValue.mLow, inV2.mValue.mLow);
	float32x4_t mul_high = vmulq_f32(mValue.mHigh, inV2.mValue.mHigh);
	float32x4_t sum4 = vaddq_f32(mul_low, mul_high);
	return vaddvq_f32(sum4);
#else
	float sum = 0.0f;
	for (int i = 0; i < 8; i++)
		sum += mF32[i] * inV2.mF32[i];
	return sum;
#endif
}

float Vec8::LengthSq() const
{
	return Dot(*this);
}

float Vec8::Length() const
{
	return sqrt(LengthSq());
}

Vec8 Vec8::Normalized() const
{
	return *this * (1.0f / Length());
}

bool Vec8::IsNormalized(float inTolerance) const
{
	return abs(LengthSq() - 1.0f) <= inTolerance;
}

JPH_NAMESPACE_END
