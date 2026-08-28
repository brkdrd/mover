// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Vec3.h>
#include <Jolt/Math/Lane4.h>
#include <Jolt/Math/Quat.h>
#include <Jolt/Math/Rotor.h>

JPH_NAMESPACE_BEGIN

#define JPH_EL(r, c) mCol[c].mF32[r]

Mat44::Mat44(Lane4Arg inC1, Lane4Arg inC2, Lane4Arg inC3, Lane4Arg inC4) :
	mCol { inC1, inC2, inC3, inC4 }
{
}

Mat44::Mat44(Type inC1, Type inC2, Type inC3, Type inC4) :
	mCol { inC1, inC2, inC3, inC4 }
{
}

Mat44 Mat44::sZero()
{
	return Mat44(Lane4::sZero(), Lane4::sZero(), Lane4::sZero(), Lane4::sZero());
}

Mat44 Mat44::sIdentity()
{
	return Mat44(Lane4(1, 0, 0, 0), Lane4(0, 1, 0, 0), Lane4(0, 0, 1, 0), Lane4(0, 0, 0, 1));
}

Mat44 Mat44::sNaN()
{
	return Mat44(Lane4::sNaN(), Lane4::sNaN(), Lane4::sNaN(), Lane4::sNaN());
}

Mat44 Mat44::sLoadFloat4x4(const Float4 *inV)
{
	Mat44 result;
	for (int c = 0; c < 4; ++c)
		result.mCol[c] = Lane4::sLoadFloat4(inV + c);
	return result;
}

Mat44 Mat44::sLoadFloat4x4Aligned(const Float4 *inV)
{
	Mat44 result;
	for (int c = 0; c < 4; ++c)
		result.mCol[c] = Lane4::sLoadFloat4Aligned(inV + c);
	return result;
}

Mat44 Mat44::sRotation(QuatArg inQuat)
{
	JPH_ASSERT(inQuat.IsNormalized());

	// See: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation section 'Quaternion-derived rotation matrix'
#ifdef JPH_USE_SSE4_1
	__m128 xyzw = inQuat.mValue.mValue;
	__m128 two_xyzw = _mm_add_ps(xyzw, xyzw);
	__m128 yzxw = _mm_shuffle_ps(xyzw, xyzw, _MM_SHUFFLE(3, 0, 2, 1));
	__m128 two_yzxw = _mm_add_ps(yzxw, yzxw);
	__m128 zxyw = _mm_shuffle_ps(xyzw, xyzw, _MM_SHUFFLE(3, 1, 0, 2));
	__m128 two_zxyw = _mm_add_ps(zxyw, zxyw);
	__m128 wwww = _mm_shuffle_ps(xyzw, xyzw, _MM_SHUFFLE(3, 3, 3, 3));
	__m128 diagonal = _mm_sub_ps(_mm_sub_ps(_mm_set1_ps(1.0f), _mm_mul_ps(two_yzxw, yzxw)), _mm_mul_ps(two_zxyw, zxyw));	// (1 - 2 y^2 - 2 z^2, 1 - 2 x^2 - 2 z^2, 1 - 2 x^2 - 2 y^2, 1 - 4 w^2)
	__m128 plus = _mm_add_ps(_mm_mul_ps(two_xyzw, zxyw), _mm_mul_ps(two_yzxw, wwww));										// 2 * (xz + yw, xy + zw, yz + xw, ww)
	__m128 minus = _mm_sub_ps(_mm_mul_ps(two_yzxw, xyzw), _mm_mul_ps(two_zxyw, wwww));										// 2 * (xy - zw, yz - xw, xz - yw, 0)

	// Workaround for compiler changing _mm_sub_ps(_mm_mul_ps(...), ...) into a fused multiply sub instruction, resulting in w not being 0
	// There doesn't appear to be a reliable way to turn this off in Clang
	minus = _mm_insert_ps(minus, minus, 0b1000);

	__m128 col0 = _mm_blend_ps(_mm_blend_ps(plus, diagonal, 0b0001), minus, 0b1100);	// (1 - 2 y^2 - 2 z^2, 2 xy + 2 zw, 2 xz - 2 yw, 0)
	__m128 col1 = _mm_blend_ps(_mm_blend_ps(diagonal, minus, 0b1001), plus, 0b0100);	// (2 xy - 2 zw, 1 - 2 x^2 - 2 z^2, 2 yz + 2 xw, 0)
	__m128 col2 = _mm_blend_ps(_mm_blend_ps(minus, plus, 0b0001), diagonal, 0b0100);	// (2 xz + 2 yw, 2 yz - 2 xw, 1 - 2 x^2 - 2 y^2, 0)
	__m128 col3 = _mm_set_ps(1, 0, 0, 0);

	return Mat44(col0, col1, col2, col3);
#else
	float x = inQuat.GetX();
	float y = inQuat.GetY();
	float z = inQuat.GetZ();
	float w = inQuat.GetW();

	float tx = x + x; // Note: Using x + x instead of 2.0f * x to force this function to return the same value as the SSE4.1 version across platforms.
	float ty = y + y;
	float tz = z + z;

	float xx = tx * x;
	float yy = ty * y;
	float zz = tz * z;
	float xy = tx * y;
	float xz = tx * z;
	float xw = tx * w;
	float yz = ty * z;
	float yw = ty * w;
	float zw = tz * w;

	return Mat44(Lane4((1.0f - yy) - zz, xy + zw, xz - yw, 0.0f), // Note: Added extra brackets to force this function to return the same value as the SSE4.1 version across platforms.
				 Lane4(xy - zw, (1.0f - zz) - xx, yz + xw, 0.0f),
				 Lane4(xz + yw, yz - xw, (1.0f - xx) - yy, 0.0f),
				 Lane4(0.0f, 0.0f, 0.0f, 1.0f));
#endif
}

Mat44 Mat44::sRotation(RotorArg inRotor)
{
	// Build 4x4 rotation matrix by applying sandwich product to each basis vector
	Lane4 c0 = inRotor * Vec4::sAxisX();
	Lane4 c1 = inRotor * Vec4::sAxisY();
	Lane4 c2 = inRotor * Vec4::sAxisZ();
	Lane4 c3 = inRotor * Vec4::sAxisW();
	return Mat44(c0, c1, c2, c3);
}

Mat44 Mat44::sRotationXY(float inAngle)
{
	float c = cos(inAngle), s = sin(inAngle);
	return Mat44(Lane4( c, s, 0, 0),
				 Lane4(-s, c, 0, 0),
				 Lane4( 0, 0, 1, 0),
				 Lane4( 0, 0, 0, 1));
}

Mat44 Mat44::sRotationXZ(float inAngle)
{
	float c = cos(inAngle), s = sin(inAngle);
	return Mat44(Lane4( c, 0, s, 0),
				 Lane4( 0, 1, 0, 0),
				 Lane4(-s, 0, c, 0),
				 Lane4( 0, 0, 0, 1));
}

Mat44 Mat44::sRotationXW(float inAngle)
{
	float c = cos(inAngle), s = sin(inAngle);
	return Mat44(Lane4( c, 0, 0, s),
				 Lane4( 0, 1, 0, 0),
				 Lane4( 0, 0, 1, 0),
				 Lane4(-s, 0, 0, c));
}

Mat44 Mat44::sRotationYZ(float inAngle)
{
	float c = cos(inAngle), s = sin(inAngle);
	return Mat44(Lane4( 1, 0, 0, 0),
				 Lane4( 0, c, s, 0),
				 Lane4( 0,-s, c, 0),
				 Lane4( 0, 0, 0, 1));
}

Mat44 Mat44::sRotationYW(float inAngle)
{
	float c = cos(inAngle), s = sin(inAngle);
	return Mat44(Lane4( 1, 0, 0, 0),
				 Lane4( 0, c, 0, s),
				 Lane4( 0, 0, 1, 0),
				 Lane4( 0,-s, 0, c));
}

Mat44 Mat44::sRotationZW(float inAngle)
{
	float c = cos(inAngle), s = sin(inAngle);
	return Mat44(Lane4( 1, 0, 0, 0),
				 Lane4( 0, 1, 0, 0),
				 Lane4( 0, 0, c, s),
				 Lane4( 0, 0,-s, c));
}

bool Mat44::operator == (Mat44Arg inM2) const
{
	return UVec4::sAnd(
		UVec4::sAnd(Lane4::sEquals(mCol[0], inM2.mCol[0]), Lane4::sEquals(mCol[1], inM2.mCol[1])),
		UVec4::sAnd(Lane4::sEquals(mCol[2], inM2.mCol[2]), Lane4::sEquals(mCol[3], inM2.mCol[3]))
	).TestAllTrue();
}

bool Mat44::IsClose(Mat44Arg inM2, float inMaxDistSq) const
{
	for (int i = 0; i < 4; ++i)
		if (!mCol[i].IsClose(inM2.mCol[i], inMaxDistSq))
			return false;
	return true;
}

Mat44 Mat44::operator * (Mat44Arg inM) const
{
	Mat44 result;
#if defined(JPH_USE_SSE)
	for (int i = 0; i < 4; ++i)
	{
		__m128 c = inM.mCol[i].mValue;
		__m128 t = _mm_mul_ps(mCol[0].mValue, _mm_shuffle_ps(c, c, _MM_SHUFFLE(0, 0, 0, 0)));
		t = _mm_add_ps(t, _mm_mul_ps(mCol[1].mValue, _mm_shuffle_ps(c, c, _MM_SHUFFLE(1, 1, 1, 1))));
		t = _mm_add_ps(t, _mm_mul_ps(mCol[2].mValue, _mm_shuffle_ps(c, c, _MM_SHUFFLE(2, 2, 2, 2))));
		t = _mm_add_ps(t, _mm_mul_ps(mCol[3].mValue, _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 3, 3, 3))));
		result.mCol[i].mValue = t;
	}
#elif defined(JPH_USE_NEON)
	for (int i = 0; i < 4; ++i)
	{
		Type c = inM.mCol[i].mValue;
		Type t = vmulq_f32(mCol[0].mValue, vdupq_laneq_f32(c, 0));
		t = vmlaq_f32(t, mCol[1].mValue, vdupq_laneq_f32(c, 1));
		t = vmlaq_f32(t, mCol[2].mValue, vdupq_laneq_f32(c, 2));
		t = vmlaq_f32(t, mCol[3].mValue, vdupq_laneq_f32(c, 3));
		result.mCol[i].mValue = t;
	}
#elif defined(JPH_USE_RVV)
	for (int i = 0; i < 4; ++i)
	{
		const float *c = inM.mCol[i].mF32;
		const vfloat32m1_t rep_0 = __riscv_vfmv_v_f_f32m1(c[0], 4);
		const vfloat32m1_t rep_1 = __riscv_vfmv_v_f_f32m1(c[1], 4);
		const vfloat32m1_t rep_2 = __riscv_vfmv_v_f_f32m1(c[2], 4);
		const vfloat32m1_t rep_3 = __riscv_vfmv_v_f_f32m1(c[3], 4);

		const vfloat32m1_t col0 = __riscv_vle32_v_f32m1(mCol[0].mF32, 4);
		const vfloat32m1_t col1 = __riscv_vle32_v_f32m1(mCol[1].mF32, 4);
		const vfloat32m1_t col2 = __riscv_vle32_v_f32m1(mCol[2].mF32, 4);
		const vfloat32m1_t col3 = __riscv_vle32_v_f32m1(mCol[3].mF32, 4);

		const vfloat32m1_t mul1 = __riscv_vfmul_vv_f32m1(col1, rep_1, 4);
		const vfloat32m1_t mul2 = __riscv_vfmul_vv_f32m1(col2, rep_2, 4);
		const vfloat32m1_t mul3 = __riscv_vfmul_vv_f32m1(col3, rep_3, 4);

		vfloat32m1_t t = __riscv_vfmul_vv_f32m1(col0, rep_0, 4);
		t = __riscv_vfadd_vv_f32m1(t, mul1, 4);
		t = __riscv_vfadd_vv_f32m1(t, mul2, 4);
		t = __riscv_vfadd_vv_f32m1(t, mul3, 4);
		__riscv_vse32_v_f32m1(result.mCol[i].mF32, t, 4);
	}
#else
	for (int i = 0; i < 4; ++i)
		result.mCol[i] = mCol[0] * inM.mCol[i].mF32[0] + mCol[1] * inM.mCol[i].mF32[1] + mCol[2] * inM.mCol[i].mF32[2] + mCol[3] * inM.mCol[i].mF32[3];
#endif
	return result;
}

Lane4 Mat44::operator * (Lane4Arg inV) const
{
#if defined(JPH_USE_SSE)
	__m128 t = _mm_mul_ps(mCol[0].mValue, _mm_shuffle_ps(inV.mValue, inV.mValue, _MM_SHUFFLE(0, 0, 0, 0)));
	t = _mm_add_ps(t, _mm_mul_ps(mCol[1].mValue, _mm_shuffle_ps(inV.mValue, inV.mValue, _MM_SHUFFLE(1, 1, 1, 1))));
	t = _mm_add_ps(t, _mm_mul_ps(mCol[2].mValue, _mm_shuffle_ps(inV.mValue, inV.mValue, _MM_SHUFFLE(2, 2, 2, 2))));
	t = _mm_add_ps(t, _mm_mul_ps(mCol[3].mValue, _mm_shuffle_ps(inV.mValue, inV.mValue, _MM_SHUFFLE(3, 3, 3, 3))));
	return t;
#elif defined(JPH_USE_NEON)
	Type t = vmulq_f32(mCol[0].mValue, vdupq_laneq_f32(inV.mValue, 0));
	t = vmlaq_f32(t, mCol[1].mValue, vdupq_laneq_f32(inV.mValue, 1));
	t = vmlaq_f32(t, mCol[2].mValue, vdupq_laneq_f32(inV.mValue, 2));
	t = vmlaq_f32(t, mCol[3].mValue, vdupq_laneq_f32(inV.mValue, 3));
	return t;
#elif defined(JPH_USE_RVV)
	const vfloat32m1_t v0 = __riscv_vfmv_v_f_f32m1(inV.mF32[0], 4);
	const vfloat32m1_t v1 = __riscv_vfmv_v_f_f32m1(inV.mF32[1], 4);
	const vfloat32m1_t v2 = __riscv_vfmv_v_f_f32m1(inV.mF32[2], 4);
	const vfloat32m1_t v3 = __riscv_vfmv_v_f_f32m1(inV.mF32[3], 4);

	const vfloat32m1_t col0 = __riscv_vle32_v_f32m1(mCol[0].mF32, 4);
	const vfloat32m1_t col1 = __riscv_vle32_v_f32m1(mCol[1].mF32, 4);
	const vfloat32m1_t col2 = __riscv_vle32_v_f32m1(mCol[2].mF32, 4);
	const vfloat32m1_t col3 = __riscv_vle32_v_f32m1(mCol[3].mF32, 4);

	const vfloat32m1_t mul1 = __riscv_vfmul_vv_f32m1(col1, v1, 4);
	const vfloat32m1_t mul2 = __riscv_vfmul_vv_f32m1(col2, v2, 4);
	const vfloat32m1_t mul3 = __riscv_vfmul_vv_f32m1(col3, v3, 4);

	vfloat32m1_t t = __riscv_vfmul_vv_f32m1(col0, v0, 4);
	t = __riscv_vfadd_vv_f32m1(t, mul1, 4);
	t = __riscv_vfadd_vv_f32m1(t, mul2, 4);
	t = __riscv_vfadd_vv_f32m1(t, mul3, 4);

	Lane4 v;
	__riscv_vse32_v_f32m1(v.mF32, t, 4);
	return v;
#else
	return Lane4(
		mCol[0].mF32[0] * inV.mF32[0] + mCol[1].mF32[0] * inV.mF32[1] + mCol[2].mF32[0] * inV.mF32[2] + mCol[3].mF32[0] * inV.mF32[3],
		mCol[0].mF32[1] * inV.mF32[0] + mCol[1].mF32[1] * inV.mF32[1] + mCol[2].mF32[1] * inV.mF32[2] + mCol[3].mF32[1] * inV.mF32[3],
		mCol[0].mF32[2] * inV.mF32[0] + mCol[1].mF32[2] * inV.mF32[1] + mCol[2].mF32[2] * inV.mF32[2] + mCol[3].mF32[2] * inV.mF32[3],
		mCol[0].mF32[3] * inV.mF32[0] + mCol[1].mF32[3] * inV.mF32[1] + mCol[2].mF32[3] * inV.mF32[2] + mCol[3].mF32[3] * inV.mF32[3]);
#endif
}

Mat44 Mat44::operator * (float inV) const
{
	Lane4 multiplier = Lane4::sReplicate(inV);

	Mat44 result;
	for (int c = 0; c < 4; ++c)
		result.mCol[c] = mCol[c] * multiplier;
	return result;
}

Mat44 &Mat44::operator *= (float inV)
{
	for (int c = 0; c < 4; ++c)
		mCol[c] *= inV;

	return *this;
}

Mat44 Mat44::operator + (Mat44Arg inM) const
{
	Mat44 result;
	for (int i = 0; i < 4; ++i)
		result.mCol[i] = mCol[i] + inM.mCol[i];
	return result;
}

Mat44 Mat44::operator - () const
{
	Mat44 result;
	for (int i = 0; i < 4; ++i)
		result.mCol[i] = -mCol[i];
	return result;
}

Mat44 Mat44::operator - (Mat44Arg inM) const
{
	Mat44 result;
	for (int i = 0; i < 4; ++i)
		result.mCol[i] = mCol[i] - inM.mCol[i];
	return result;
}

Mat44 &Mat44::operator += (Mat44Arg inM)
{
	for (int c = 0; c < 4; ++c)
		mCol[c] += inM.mCol[c];

	return *this;
}

void Mat44::StoreFloat4x4(Float4 *outV) const
{
	for (int c = 0; c < 4; ++c)
		mCol[c].StoreFloat4(outV + c);
}

Mat44 Mat44::Transposed() const
{
#if defined(JPH_USE_SSE)
	__m128 tmp1 = _mm_shuffle_ps(mCol[0].mValue, mCol[1].mValue, _MM_SHUFFLE(1, 0, 1, 0));
	__m128 tmp3 = _mm_shuffle_ps(mCol[0].mValue, mCol[1].mValue, _MM_SHUFFLE(3, 2, 3, 2));
	__m128 tmp2 = _mm_shuffle_ps(mCol[2].mValue, mCol[3].mValue, _MM_SHUFFLE(1, 0, 1, 0));
	__m128 tmp4 = _mm_shuffle_ps(mCol[2].mValue, mCol[3].mValue, _MM_SHUFFLE(3, 2, 3, 2));

	Mat44 result;
	result.mCol[0].mValue = _mm_shuffle_ps(tmp1, tmp2, _MM_SHUFFLE(2, 0, 2, 0));
	result.mCol[1].mValue = _mm_shuffle_ps(tmp1, tmp2, _MM_SHUFFLE(3, 1, 3, 1));
	result.mCol[2].mValue = _mm_shuffle_ps(tmp3, tmp4, _MM_SHUFFLE(2, 0, 2, 0));
	result.mCol[3].mValue = _mm_shuffle_ps(tmp3, tmp4, _MM_SHUFFLE(3, 1, 3, 1));
	return result;
#elif defined(JPH_USE_NEON)
	float32x4x2_t tmp1 = vzipq_f32(mCol[0].mValue, mCol[2].mValue);
	float32x4x2_t tmp2 = vzipq_f32(mCol[1].mValue, mCol[3].mValue);
	float32x4x2_t tmp3 = vzipq_f32(tmp1.val[0], tmp2.val[0]);
	float32x4x2_t tmp4 = vzipq_f32(tmp1.val[1], tmp2.val[1]);

	Mat44 result;
	result.mCol[0].mValue = tmp3.val[0];
	result.mCol[1].mValue = tmp3.val[1];
	result.mCol[2].mValue = tmp4.val[0];
	result.mCol[3].mValue = tmp4.val[1];
	return result;
#elif defined(JPH_USE_RVV)
	const vfloat32m1_t row0 = __riscv_vlse32_v_f32m1(&mCol[0].mF32[0], sizeof(Lane4), 4);
	const vfloat32m1_t row1 = __riscv_vlse32_v_f32m1(&mCol[0].mF32[1], sizeof(Lane4), 4);
	const vfloat32m1_t row2 = __riscv_vlse32_v_f32m1(&mCol[0].mF32[2], sizeof(Lane4), 4);
	const vfloat32m1_t row3 = __riscv_vlse32_v_f32m1(&mCol[0].mF32[3], sizeof(Lane4), 4);

	Mat44 result;
	__riscv_vse32_v_f32m1(result.mCol[0].mF32, row0, 4);
	__riscv_vse32_v_f32m1(result.mCol[1].mF32, row1, 4);
	__riscv_vse32_v_f32m1(result.mCol[2].mF32, row2, 4);
	__riscv_vse32_v_f32m1(result.mCol[3].mF32, row3, 4);
	return result;
#else
	Mat44 result;
	for (int c = 0; c < 4; ++c)
		for (int r = 0; r < 4; ++r)
			result.mCol[r].mF32[c] = mCol[c].mF32[r];
	return result;
#endif
}

float Mat44::GetDeterminant() const
{
	// Laplace expansion along the first two rows, pairing each 2x2 minor of rows 0/1 with its
	// complementary 2x2 minor of rows 2/3
	float m00 = mCol[0].mF32[0], m01 = mCol[1].mF32[0], m02 = mCol[2].mF32[0], m03 = mCol[3].mF32[0];
	float m10 = mCol[0].mF32[1], m11 = mCol[1].mF32[1], m12 = mCol[2].mF32[1], m13 = mCol[3].mF32[1];
	float m20 = mCol[0].mF32[2], m21 = mCol[1].mF32[2], m22 = mCol[2].mF32[2], m23 = mCol[3].mF32[2];
	float m30 = mCol[0].mF32[3], m31 = mCol[1].mF32[3], m32 = mCol[2].mF32[3], m33 = mCol[3].mF32[3];

	float s0 = m00 * m11 - m10 * m01;
	float s1 = m00 * m12 - m10 * m02;
	float s2 = m00 * m13 - m10 * m03;
	float s3 = m01 * m12 - m11 * m02;
	float s4 = m01 * m13 - m11 * m03;
	float s5 = m02 * m13 - m12 * m03;

	float c5 = m22 * m33 - m32 * m23;
	float c4 = m21 * m33 - m31 * m23;
	float c3 = m21 * m32 - m31 * m22;
	float c2 = m20 * m33 - m30 * m23;
	float c1 = m20 * m32 - m30 * m22;
	float c0 = m20 * m31 - m30 * m21;

	return s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
}

Mat44 Mat44::Inversed() const
{
#if defined(JPH_USE_SSE)
	// Algorithm from: http://download.intel.com/design/PentiumIII/sml/24504301.pdf
	// Streaming SIMD Extensions - Inverse of 4x4 Matrix
	// Adapted to load data using _mm_shuffle_ps instead of loading from memory
	// Replaced _mm_rcp_ps with _mm_div_ps for better accuracy

	__m128 tmp1 = _mm_shuffle_ps(mCol[0].mValue, mCol[1].mValue, _MM_SHUFFLE(1, 0, 1, 0));
	__m128 row1 = _mm_shuffle_ps(mCol[2].mValue, mCol[3].mValue, _MM_SHUFFLE(1, 0, 1, 0));
	__m128 row0 = _mm_shuffle_ps(tmp1, row1, _MM_SHUFFLE(2, 0, 2, 0));
	row1 = _mm_shuffle_ps(row1, tmp1, _MM_SHUFFLE(3, 1, 3, 1));
	tmp1 = _mm_shuffle_ps(mCol[0].mValue, mCol[1].mValue, _MM_SHUFFLE(3, 2, 3, 2));
	__m128 row3 = _mm_shuffle_ps(mCol[2].mValue, mCol[3].mValue, _MM_SHUFFLE(3, 2, 3, 2));
	__m128 row2 = _mm_shuffle_ps(tmp1, row3, _MM_SHUFFLE(2, 0, 2, 0));
	row3 = _mm_shuffle_ps(row3, tmp1, _MM_SHUFFLE(3, 1, 3, 1));

	tmp1 = _mm_mul_ps(row2, row3);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2, 3, 0, 1));
	__m128 minor0 = _mm_mul_ps(row1, tmp1);
	__m128 minor1 = _mm_mul_ps(row0, tmp1);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(1, 0, 3, 2));
	minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
	minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
	minor1 = _mm_shuffle_ps(minor1, minor1, _MM_SHUFFLE(1, 0, 3, 2));

	tmp1 = _mm_mul_ps(row1, row2);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2, 3, 0, 1));
	minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
	__m128 minor3 = _mm_mul_ps(row0, tmp1);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(1, 0, 3, 2));
	minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
	minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
	minor3 = _mm_shuffle_ps(minor3, minor3, _MM_SHUFFLE(1, 0, 3, 2));

	tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, _MM_SHUFFLE(1, 0, 3, 2)), row3);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2, 3, 0, 1));
	row2 = _mm_shuffle_ps(row2, row2, _MM_SHUFFLE(1, 0, 3, 2));
	minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
	__m128 minor2 = _mm_mul_ps(row0, tmp1);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(1, 0, 3, 2));
	minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
	minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
	minor2 = _mm_shuffle_ps(minor2, minor2, _MM_SHUFFLE(1, 0, 3, 2));

	tmp1 = _mm_mul_ps(row0, row1);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2, 3, 0, 1));
	minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
	minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(1, 0, 3, 2));
	minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
	minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));

	tmp1 = _mm_mul_ps(row0, row3);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2, 3, 0, 1));
	minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
	minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(1, 0, 3, 2));
	minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
	minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));

	tmp1 = _mm_mul_ps(row0, row2);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2, 3, 0, 1));
	minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
	minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(1, 0, 3, 2));
	minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
	minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);

	__m128 det = _mm_mul_ps(row0, minor0);
	det = _mm_add_ps(_mm_shuffle_ps(det, det, _MM_SHUFFLE(2, 3, 0, 1)), det); // Original code did (x + z) + (y + w), changed to (x + y) + (z + w) to match the ARM code below and make the result cross platform deterministic
	det = _mm_add_ss(_mm_shuffle_ps(det, det, _MM_SHUFFLE(1, 0, 3, 2)), det);
	det = _mm_div_ss(_mm_set_ss(1.0f), det);
	det = _mm_shuffle_ps(det, det, _MM_SHUFFLE(0, 0, 0, 0));

	Mat44 result;
	result.mCol[0].mValue = _mm_mul_ps(det, minor0);
	result.mCol[1].mValue = _mm_mul_ps(det, minor1);
	result.mCol[2].mValue = _mm_mul_ps(det, minor2);
	result.mCol[3].mValue = _mm_mul_ps(det, minor3);
	return result;
#elif defined(JPH_USE_NEON)
	// Adapted from the SSE version, there's surprising few articles about efficient ways of calculating an inverse for ARM on the internet
	Type tmp1 = JPH_NEON_SHUFFLE_F32x4(mCol[0].mValue, mCol[1].mValue, 0, 1, 4, 5);
	Type row1 = JPH_NEON_SHUFFLE_F32x4(mCol[2].mValue, mCol[3].mValue, 0, 1, 4, 5);
	Type row0 = JPH_NEON_SHUFFLE_F32x4(tmp1, row1, 0, 2, 4, 6);
	row1 = JPH_NEON_SHUFFLE_F32x4(row1, tmp1, 1, 3, 5, 7);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(mCol[0].mValue, mCol[1].mValue, 2, 3, 6, 7);
	Type row3 = JPH_NEON_SHUFFLE_F32x4(mCol[2].mValue, mCol[3].mValue, 2, 3, 6, 7);
	Type row2 = JPH_NEON_SHUFFLE_F32x4(tmp1, row3, 0, 2, 4, 6);
	row3 = JPH_NEON_SHUFFLE_F32x4(row3, tmp1, 1, 3, 5, 7);

	tmp1 = vmulq_f32(row2, row3);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 1, 0, 3, 2);
	Type minor0 = vmulq_f32(row1, tmp1);
	Type minor1 = vmulq_f32(row0, tmp1);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 2, 3, 0, 1);
	minor0 = vsubq_f32(vmulq_f32(row1, tmp1), minor0);
	minor1 = vsubq_f32(vmulq_f32(row0, tmp1), minor1);
	minor1 = JPH_NEON_SHUFFLE_F32x4(minor1, minor1, 2, 3, 0, 1);

	tmp1 = vmulq_f32(row1, row2);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 1, 0, 3, 2);
	minor0 = vaddq_f32(vmulq_f32(row3, tmp1), minor0);
	Type minor3 = vmulq_f32(row0, tmp1);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 2, 3, 0, 1);
	minor0 = vsubq_f32(minor0, vmulq_f32(row3, tmp1));
	minor3 = vsubq_f32(vmulq_f32(row0, tmp1), minor3);
	minor3 = JPH_NEON_SHUFFLE_F32x4(minor3, minor3, 2, 3, 0, 1);

	tmp1 = JPH_NEON_SHUFFLE_F32x4(row1, row1, 2, 3, 0, 1);
	tmp1 = vmulq_f32(tmp1, row3);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 1, 0, 3, 2);
	row2 = JPH_NEON_SHUFFLE_F32x4(row2, row2, 2, 3, 0, 1);
	minor0 = vaddq_f32(vmulq_f32(row2, tmp1), minor0);
	Type minor2 = vmulq_f32(row0, tmp1);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 2, 3, 0, 1);
	minor0 = vsubq_f32(minor0, vmulq_f32(row2, tmp1));
	minor2 = vsubq_f32(vmulq_f32(row0, tmp1), minor2);
	minor2 = JPH_NEON_SHUFFLE_F32x4(minor2, minor2, 2, 3, 0, 1);

	tmp1 = vmulq_f32(row0, row1);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 1, 0, 3, 2);
	minor2 = vaddq_f32(vmulq_f32(row3, tmp1), minor2);
	minor3 = vsubq_f32(vmulq_f32(row2, tmp1), minor3);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 2, 3, 0, 1);
	minor2 = vsubq_f32(vmulq_f32(row3, tmp1), minor2);
	minor3 = vsubq_f32(minor3, vmulq_f32(row2, tmp1));

	tmp1 = vmulq_f32(row0, row3);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 1, 0, 3, 2);
	minor1 = vsubq_f32(minor1, vmulq_f32(row2, tmp1));
	minor2 = vaddq_f32(vmulq_f32(row1, tmp1), minor2);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 2, 3, 0, 1);
	minor1 = vaddq_f32(vmulq_f32(row2, tmp1), minor1);
	minor2 = vsubq_f32(minor2, vmulq_f32(row1, tmp1));

	tmp1 = vmulq_f32(row0, row2);
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 1, 0, 3, 2);
	minor1 = vaddq_f32(vmulq_f32(row3, tmp1), minor1);
	minor3 = vsubq_f32(minor3, vmulq_f32(row1, tmp1));
	tmp1 = JPH_NEON_SHUFFLE_F32x4(tmp1, tmp1, 2, 3, 0, 1);
	minor1 = vsubq_f32(minor1, vmulq_f32(row3, tmp1));
	minor3 = vaddq_f32(vmulq_f32(row1, tmp1), minor3);

	Type det = vmulq_f32(row0, minor0);
	det = vdupq_n_f32(vaddvq_f32(det));
	det = vdivq_f32(vdupq_n_f32(1.0f), det);

	Mat44 result;
	result.mCol[0].mValue = vmulq_f32(det, minor0);
	result.mCol[1].mValue = vmulq_f32(det, minor1);
	result.mCol[2].mValue = vmulq_f32(det, minor2);
	result.mCol[3].mValue = vmulq_f32(det, minor3);
	return result;
#elif defined(JPH_USE_RVV)
	// Implementation mirrored from SSE and NEON implementations
	const vfloat32m1_t c0 = __riscv_vle32_v_f32m1(mCol[0].mF32, 4);
	const vfloat32m1_t c1 = __riscv_vle32_v_f32m1(mCol[1].mF32, 4);
	const vfloat32m1_t c2 = __riscv_vle32_v_f32m1(mCol[2].mF32, 4);
	const vfloat32m1_t c3 = __riscv_vle32_v_f32m1(mCol[3].mF32, 4);

	vfloat32m1_t minor0, minor1, minor2, minor3;
	vfloat32m1_t tmp1;
	vfloat32m1_t row0, row1, row2, row3;

	tmp1 = RVVShuffleFloat32x4<0, 1, 4, 5>(c0, c1);
	row1 = RVVShuffleFloat32x4<0, 1, 4, 5>(c2, c3);
	row0 = RVVShuffleFloat32x4<0, 2, 4, 6>(tmp1, row1);
	row1 = RVVShuffleFloat32x4<1, 3, 5, 7>(row1, tmp1);
	tmp1 = RVVShuffleFloat32x4<2, 3, 6, 7>(c0, c1);
	row3 = RVVShuffleFloat32x4<2, 3, 6, 7>(c2, c3);
	row2 = RVVShuffleFloat32x4<0, 2, 4, 6>(tmp1, row3);
	row3 = RVVShuffleFloat32x4<1, 3, 5, 7>(row3, tmp1);

	tmp1 = __riscv_vfmul_vv_f32m1(row2, row3, 4);
	tmp1 = RVVShuffleFloat32x4<1, 0, 3, 2>(tmp1, tmp1);
	minor0 = __riscv_vfmul_vv_f32m1(row1, tmp1, 4);
	minor1 = __riscv_vfmul_vv_f32m1(row0, tmp1, 4);
	tmp1 = RVVShuffleFloat32x4<2, 3, 0, 1>(tmp1, tmp1);
	minor0 = __riscv_vfsub_vv_f32m1(__riscv_vfmul_vv_f32m1(row1, tmp1, 4), minor0, 4);
	minor1 = __riscv_vfsub_vv_f32m1(__riscv_vfmul_vv_f32m1(row0, tmp1, 4), minor1, 4);
	minor1 = RVVShuffleFloat32x4<2, 3, 0, 1>(minor1, minor1);

	tmp1 = __riscv_vfmul_vv_f32m1(row1, row2, 4);
	tmp1 = RVVShuffleFloat32x4<1, 0, 3, 2>(tmp1, tmp1);
	minor0 = __riscv_vfadd_vv_f32m1(__riscv_vfmul_vv_f32m1(row3, tmp1, 4), minor0, 4);
	minor3 = __riscv_vfmul_vv_f32m1(row0, tmp1, 4);
	tmp1 = RVVShuffleFloat32x4<2, 3, 0, 1>(tmp1, tmp1);
	minor0 = __riscv_vfsub_vv_f32m1(minor0, __riscv_vfmul_vv_f32m1(row3, tmp1, 4), 4);
	minor3 = __riscv_vfsub_vv_f32m1(__riscv_vfmul_vv_f32m1(row0, tmp1, 4), minor3, 4);
	minor3 = RVVShuffleFloat32x4<2, 3, 0, 1>(minor3, minor3);

	tmp1 = RVVShuffleFloat32x4<2, 3, 0, 1>(row1, row1);
	tmp1 = __riscv_vfmul_vv_f32m1(tmp1, row3, 4);
	tmp1 = RVVShuffleFloat32x4<1, 0, 3, 2>(tmp1, tmp1);
	row2 = RVVShuffleFloat32x4<2, 3, 0, 1>(row2, row2);
	minor0 = __riscv_vfadd_vv_f32m1(__riscv_vfmul_vv_f32m1(row2, tmp1, 4), minor0, 4);
	minor2 = __riscv_vfmul_vv_f32m1(row0, tmp1, 4);
	tmp1 = RVVShuffleFloat32x4<2, 3, 0, 1>(tmp1, tmp1);
	minor0 = __riscv_vfsub_vv_f32m1(minor0, __riscv_vfmul_vv_f32m1(row2, tmp1, 4), 4);
	minor2 = __riscv_vfsub_vv_f32m1(__riscv_vfmul_vv_f32m1(row0, tmp1, 4), minor2, 4);
	minor2 = RVVShuffleFloat32x4<2, 3, 0, 1>(minor2, minor2);

	tmp1 = __riscv_vfmul_vv_f32m1(row0, row1, 4);
	tmp1 = RVVShuffleFloat32x4<1, 0, 3, 2>(tmp1, tmp1);
	minor2 = __riscv_vfadd_vv_f32m1(__riscv_vfmul_vv_f32m1(row3, tmp1, 4), minor2, 4);
	minor3 = __riscv_vfsub_vv_f32m1(__riscv_vfmul_vv_f32m1(row2, tmp1, 4), minor3, 4);
	tmp1 = RVVShuffleFloat32x4<2, 3, 0, 1>(tmp1, tmp1);
	minor2 = __riscv_vfsub_vv_f32m1(__riscv_vfmul_vv_f32m1(row3, tmp1, 4), minor2, 4);
	minor3 = __riscv_vfsub_vv_f32m1(minor3, __riscv_vfmul_vv_f32m1(row2, tmp1, 4), 4);

	tmp1 = __riscv_vfmul_vv_f32m1(row0, row3, 4);
	tmp1 = RVVShuffleFloat32x4<1, 0, 3, 2>(tmp1, tmp1);
	minor1 = __riscv_vfsub_vv_f32m1(minor1, __riscv_vfmul_vv_f32m1(row2, tmp1, 4), 4);
	minor2 = __riscv_vfadd_vv_f32m1(__riscv_vfmul_vv_f32m1(row1, tmp1, 4), minor2, 4);
	tmp1 = RVVShuffleFloat32x4<2, 3, 0, 1>(tmp1, tmp1);
	minor1 = __riscv_vfadd_vv_f32m1(__riscv_vfmul_vv_f32m1(row2, tmp1, 4), minor1, 4);
	minor2 = __riscv_vfsub_vv_f32m1(minor2, __riscv_vfmul_vv_f32m1(row1, tmp1, 4), 4);

	tmp1 = __riscv_vfmul_vv_f32m1(row0, row2, 4);
	tmp1 = RVVShuffleFloat32x4<1, 0, 3, 2>(tmp1, tmp1);
	minor1 = __riscv_vfadd_vv_f32m1(__riscv_vfmul_vv_f32m1(row3, tmp1, 4), minor1, 4);
	minor3 = __riscv_vfsub_vv_f32m1(minor3, __riscv_vfmul_vv_f32m1(row1, tmp1, 4), 4);
	tmp1 = RVVShuffleFloat32x4<2, 3, 0, 1>(tmp1, tmp1);
	minor1 = __riscv_vfsub_vv_f32m1(minor1, __riscv_vfmul_vv_f32m1(row3, tmp1, 4), 4);
	minor3 = __riscv_vfadd_vv_f32m1(__riscv_vfmul_vv_f32m1(row1, tmp1, 4), minor3, 4);

	const vfloat32m1_t v_det = __riscv_vfmul_vv_f32m1(row0, minor0, 4);
	const float s_det = RVVSumElementsFloat32x4(v_det);
	const vfloat32m1_t det_inv = __riscv_vfmv_v_f_f32m1(1.0f / s_det, 4);

	minor0 = __riscv_vfmul_vv_f32m1(det_inv, minor0, 4);
	minor1 = __riscv_vfmul_vv_f32m1(det_inv, minor1, 4);
	minor2 = __riscv_vfmul_vv_f32m1(det_inv, minor2, 4);
	minor3 = __riscv_vfmul_vv_f32m1(det_inv, minor3, 4);

	Mat44 result;
	__riscv_vse32_v_f32m1(result.mCol[0].mF32, minor0, 4);
	__riscv_vse32_v_f32m1(result.mCol[1].mF32, minor1, 4);
	__riscv_vse32_v_f32m1(result.mCol[2].mF32, minor2, 4);
	__riscv_vse32_v_f32m1(result.mCol[3].mF32, minor3, 4);
	return result;
#else
	float m00 = JPH_EL(0, 0), m10 = JPH_EL(1, 0), m20 = JPH_EL(2, 0), m30 = JPH_EL(3, 0);
	float m01 = JPH_EL(0, 1), m11 = JPH_EL(1, 1), m21 = JPH_EL(2, 1), m31 = JPH_EL(3, 1);
	float m02 = JPH_EL(0, 2), m12 = JPH_EL(1, 2), m22 = JPH_EL(2, 2), m32 = JPH_EL(3, 2);
	float m03 = JPH_EL(0, 3), m13 = JPH_EL(1, 3), m23 = JPH_EL(2, 3), m33 = JPH_EL(3, 3);

	float m10211120 = m10 * m21 - m11 * m20;
	float m10221220 = m10 * m22 - m12 * m20;
	float m10231320 = m10 * m23 - m13 * m20;
	float m10311130 = m10 * m31 - m11 * m30;
	float m10321230 = m10 * m32 - m12 * m30;
	float m10331330 = m10 * m33 - m13 * m30;
	float m11221221 = m11 * m22 - m12 * m21;
	float m11231321 = m11 * m23 - m13 * m21;
	float m11321231 = m11 * m32 - m12 * m31;
	float m11331331 = m11 * m33 - m13 * m31;
	float m12231322 = m12 * m23 - m13 * m22;
	float m12331332 = m12 * m33 - m13 * m32;
	float m20312130 = m20 * m31 - m21 * m30;
	float m20322230 = m20 * m32 - m22 * m30;
	float m20332330 = m20 * m33 - m23 * m30;
	float m21322231 = m21 * m32 - m22 * m31;
	float m21332331 = m21 * m33 - m23 * m31;
	float m22332332 = m22 * m33 - m23 * m32;

	Lane4 col0(m11 * m22332332 - m12 * m21332331 + m13 * m21322231,		-m10 * m22332332 + m12 * m20332330 - m13 * m20322230,		m10 * m21332331 - m11 * m20332330 + m13 * m20312130,		-m10 * m21322231 + m11 * m20322230 - m12 * m20312130);
	Lane4 col1(-m01 * m22332332 + m02 * m21332331 - m03 * m21322231,		m00 * m22332332 - m02 * m20332330 + m03 * m20322230,		-m00 * m21332331 + m01 * m20332330 - m03 * m20312130,		m00 * m21322231 - m01 * m20322230 + m02 * m20312130);
	Lane4 col2(m01 * m12331332 - m02 * m11331331 + m03 * m11321231,		-m00 * m12331332 + m02 * m10331330 - m03 * m10321230,		m00 * m11331331 - m01 * m10331330 + m03 * m10311130,		-m00 * m11321231 + m01 * m10321230 - m02 * m10311130);
	Lane4 col3(-m01 * m12231322 + m02 * m11231321 - m03 * m11221221,		m00 * m12231322 - m02 * m10231320 + m03 * m10221220,		-m00 * m11231321 + m01 * m10231320 - m03 * m10211120,		m00 * m11221221 - m01 * m10221220 + m02 * m10211120);

	float det = m00 * col0.mF32[0] + m01 * col0.mF32[1] + m02 * col0.mF32[2] + m03 * col0.mF32[3];

	return Mat44(col0 / det, col1 / det, col2 / det, col3 / det);
#endif
}

Rotor Mat44::GetRotor() const
{
	// Extract a 4D rotor from this orthogonal rotation matrix.
	// Uses quaternion isomorphism Spin(4) ≅ (SU(2)×SU(2))/Z₂:
	//   rotation acts as v ↦ q_L * v * conj(q_R)
	// where R⁴ ↔ H via (v1,v2,v3,v4) ↔ v1 + v2*i + v3*j + v4*k.

	// Step 1: Read columns as quaternions (w,x,y,z) = (row0,row1,row2,row3)
	float c0w = mCol[0].mF32[0], c0x = mCol[0].mF32[1], c0y = mCol[0].mF32[2], c0z = mCol[0].mF32[3];
	float c1w = mCol[1].mF32[0], c1x = mCol[1].mF32[1], c1y = mCol[1].mF32[2], c1z = mCol[1].mF32[3];
	float c2w = mCol[2].mF32[0], c2x = mCol[2].mF32[1], c2y = mCol[2].mF32[2], c2z = mCol[2].mF32[3];
	float c3w = mCol[3].mF32[0], c3x = mCol[3].mF32[1], c3y = mCol[3].mF32[2], c3z = mCol[3].mF32[3];

	// Step 2: conj(C0) * Ck for k=1,2,3 — imaginary parts form q_R's 3x3 rotation matrix
	float p1x = c0w*c1x - c0x*c1w - c0y*c1z + c0z*c1y;
	float p1y = c0w*c1y + c0x*c1z - c0y*c1w - c0z*c1x;
	float p1z = c0w*c1z - c0x*c1y + c0y*c1x - c0z*c1w;

	float p2x = c0w*c2x - c0x*c2w - c0y*c2z + c0z*c2y;
	float p2y = c0w*c2y + c0x*c2z - c0y*c2w - c0z*c2x;
	float p2z = c0w*c2z - c0x*c2y + c0y*c2x - c0z*c2w;

	float p3x = c0w*c3x - c0x*c3w - c0y*c3z + c0z*c3y;
	float p3y = c0w*c3y + c0x*c3z - c0y*c3w - c0z*c3x;
	float p3z = c0w*c3z - c0x*c3y + c0y*c3x - c0z*c3w;

	// Step 3: Extract q_R from R_3x3 using Shepperd's method
	// R_3x3[row][col]: col0=(p1x,p1y,p1z), col1=(p2x,p2y,p2z), col2=(p3x,p3y,p3z)
	float tr = p1x + p2y + p3z;
	float qrw, qrx, qry, qrz;

	if (tr >= 0.0f)
	{
		float s = sqrt(tr + 1.0f);
		float is = 0.5f / s;
		qrw = 0.5f * s;
		qrx = (p2z - p3y) * is;
		qry = (p3x - p1z) * is;
		qrz = (p1y - p2x) * is;
	}
	else
	{
		int i = 0;
		if (p2y > p1x) i = 1;
		if (p3z > (i == 0 ? p1x : p2y)) i = 2;

		if (i == 0)
		{
			float s = sqrt(p1x - p2y - p3z + 1.0f);
			float is = 0.5f / s;
			qrx = 0.5f * s;
			qry = (p2x + p1y) * is;
			qrz = (p1z + p3x) * is;
			qrw = (p2z - p3y) * is;
		}
		else if (i == 1)
		{
			float s = sqrt(p2y - p3z - p1x + 1.0f);
			float is = 0.5f / s;
			qrx = (p2x + p1y) * is;
			qry = 0.5f * s;
			qrz = (p3y + p2z) * is;
			qrw = (p3x - p1z) * is;
		}
		else
		{
			JPH_ASSERT(i == 2);

			float s = sqrt(p3z - p1x - p2y + 1.0f);
			float is = 0.5f / s;
			qrx = (p1z + p3x) * is;
			qry = (p3y + p2z) * is;
			qrz = 0.5f * s;
			qrw = (p1y - p2x) * is;
		}
	}

	// Step 4: q_L = C0 * q_R (quaternion multiply)
	float qlw = c0w*qrw - c0x*qrx - c0y*qry - c0z*qrz;
	float qlx = c0w*qrx + c0x*qrw + c0y*qrz - c0z*qry;
	float qly = c0w*qry - c0x*qrz + c0y*qrw + c0z*qrx;
	float qlz = c0w*qrz + c0x*qry - c0y*qrx + c0z*qrw;

	// Step 5: Convert (q_L, q_R) pair to rotor components [s, e12, e13, e14, e23, e24, e34, e1234].
	// The (qL, qR) computed above are labeled mirror to the rotor's SU(2) split convention (see
	// Rotor.inl: qL is the +e1234 eigenspace factor), i.e. their roles are swapped relative to the
	// rotor. Swapping qL<->qR negates exactly the anti-self-dual components e23/e24/e34/e1234 while
	// leaving s/e12/e13/e14 unchanged, which is applied here. Without it GetRotor returns the
	// opposite rotation in the e23/e24/e34 planes and sRotation(GetRotor(M)) != M.
	return Rotor(
		(qlw + qrw) * 0.5f,			// s
		(qrx - qlx) * 0.5f,			// e12
		(qry - qly) * 0.5f,			// e13
		(qrz - qlz) * 0.5f,			// e14
		-(qrz + qlz) * 0.5f,		// e23
		(qry + qly) * 0.5f,			// e24
		-(qrx + qlx) * 0.5f,		// e34
		(qrw - qlw) * 0.5f			// e1234
	).Normalized();
}

#undef JPH_EL

Mat44 Rotor::ToRotationMatrix() const
{
	return Mat44::sRotation(*this);
}

JPH_NAMESPACE_END
