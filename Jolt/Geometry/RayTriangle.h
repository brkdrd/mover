// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

/// Intersect ray with triangle, returns closest point or FLT_MAX if no hit.
/// In 4D, the ray-triangle system is overdetermined (4 equations, 3 unknowns).
/// We solve via normal equations A^T A x = A^T b, then check the residual.
JPH_INLINE float RayTriangle(Vec4Arg inOrigin, Vec4Arg inDirection, Vec4Arg inV0, Vec4Arg inV1, Vec4Arg inV2)
{
	// Edge vectors
	Vec4 e1 = inV1 - inV0;
	Vec4 e2 = inV2 - inV0;
	Vec4 s = inOrigin - inV0;

	// Build symmetric 3x3 Gram matrix: A^T A where A = [-D, e1, e2]
	float a00 = inDirection.LengthSq();
	float a01 = -inDirection.Dot(e1);
	float a02 = -inDirection.Dot(e2);
	float a11 = e1.LengthSq();
	float a12 = e1.Dot(e2);
	float a22 = e2.LengthSq();

	// RHS: A^T b where b = s
	float b0 = -inDirection.Dot(s);
	float b1 = e1.Dot(s);
	float b2 = e2.Dot(s);

	// Cofactors of the symmetric 3x3 matrix
	float c00 = a11 * a22 - a12 * a12;
	float c01 = a12 * a02 - a01 * a22;
	float c02 = a01 * a12 - a11 * a02;
	float c11 = a00 * a22 - a02 * a02;
	float c12 = a01 * a02 - a00 * a12;
	float c22 = a00 * a11 - a01 * a01;

	// Determinant
	float det = a00 * c00 + a01 * c01 + a02 * c02;
	if (abs(det) < 1.0e-12f)
		return FLT_MAX;

	float inv_det = 1.0f / det;

	// Cramer's rule for [t, u, v]
	float t = (b0 * c00 + b1 * c01 + b2 * c02) * inv_det;
	float u = (b0 * c01 + b1 * c11 + b2 * c12) * inv_det;
	float v = (b0 * c02 + b1 * c12 + b2 * c22) * inv_det;

	// Validate: t >= 0, u >= 0, v >= 0, u + v <= 1
	if (t < 0.0f || u < 0.0f || v < 0.0f || u + v > 1.0f)
		return FLT_MAX;

	// Check residual: ray may miss the triangle in 4D even if the projected system has a solution
	Vec4 residual = s + inDirection * t - e1 * u - e2 * v;
	float residual_sq = residual.LengthSq();
	float scale = max(e1.LengthSq(), max(e2.LengthSq(), inDirection.LengthSq() * t * t));
	if (residual_sq > 1.0e-6f * scale)
		return FLT_MAX;

	return t;
}

/// Intersect ray with 4 triangles in SOA format, returns 4 vector of closest points or FLT_MAX if no hit.
/// Uses normal equations (A^T A x = A^T b) for 4D overdetermined system with residual check.
JPH_INLINE Vec4 RayTriangle4(Vec4Arg inOrigin, Vec4Arg inDirection,
	Vec4Arg inV0X, Vec4Arg inV0Y, Vec4Arg inV0Z, Vec4Arg inV0W,
	Vec4Arg inV1X, Vec4Arg inV1Y, Vec4Arg inV1Z, Vec4Arg inV1W,
	Vec4Arg inV2X, Vec4Arg inV2Y, Vec4Arg inV2Z, Vec4Arg inV2W)
{
	Vec4 epsilon = Vec4::sReplicate(1.0e-12f);
	Vec4 residual_tol = Vec4::sReplicate(1.0e-6f);
	Vec4 zero = Vec4::sZero();

	// Edge vectors (SOA: each lane = one triangle)
	Vec4 e1x = inV1X - inV0X;
	Vec4 e1y = inV1Y - inV0Y;
	Vec4 e1z = inV1Z - inV0Z;
	Vec4 e1w = inV1W - inV0W;
	Vec4 e2x = inV2X - inV0X;
	Vec4 e2y = inV2Y - inV0Y;
	Vec4 e2z = inV2Z - inV0Z;
	Vec4 e2w = inV2W - inV0W;

	// Ray direction components splatted
	Vec4 dx = inDirection.SplatX();
	Vec4 dy = inDirection.SplatY();
	Vec4 dz = inDirection.SplatZ();
	Vec4 dw = inDirection.SplatW();

	// Vector from V0 to origin
	Vec4 sx = inOrigin.SplatX() - inV0X;
	Vec4 sy = inOrigin.SplatY() - inV0Y;
	Vec4 sz = inOrigin.SplatZ() - inV0Z;
	Vec4 sw = inOrigin.SplatW() - inV0W;

	// Build A^T A (symmetric 3x3) where A = [-D, e1, e2]
	// 4D dot products in SOA
	Vec4 a00 = dx * dx + dy * dy + dz * dz + dw * dw;		// D.D
	Vec4 a01 = -(dx * e1x + dy * e1y + dz * e1z + dw * e1w);	// -D.e1
	Vec4 a02 = -(dx * e2x + dy * e2y + dz * e2z + dw * e2w);	// -D.e2
	Vec4 a11 = e1x * e1x + e1y * e1y + e1z * e1z + e1w * e1w;	// e1.e1
	Vec4 a12 = e1x * e2x + e1y * e2y + e1z * e2z + e1w * e2w;	// e1.e2
	Vec4 a22 = e2x * e2x + e2y * e2y + e2z * e2z + e2w * e2w;	// e2.e2

	// RHS: A^T b where b = s
	Vec4 b0 = -(dx * sx + dy * sy + dz * sz + dw * sw);		// -D.s
	Vec4 b1 = e1x * sx + e1y * sy + e1z * sz + e1w * sw;		// e1.s
	Vec4 b2 = e2x * sx + e2y * sy + e2z * sz + e2w * sw;		// e2.s

	// Cofactors of symmetric 3x3
	Vec4 c00 = a11 * a22 - a12 * a12;
	Vec4 c01 = a12 * a02 - a01 * a22;
	Vec4 c02 = a01 * a12 - a11 * a02;
	Vec4 c11 = a00 * a22 - a02 * a02;
	Vec4 c12 = a01 * a02 - a00 * a12;
	Vec4 c22 = a00 * a11 - a01 * a01;

	// Determinant
	Vec4 det = a00 * c00 + a01 * c01 + a02 * c02;

	// Get sign bit and make positive
	Vec4 det_sign = Vec4::sAnd(det, UVec4::sReplicate(0x80000000).ReinterpretAsFloat());
	Vec4 abs_det = Vec4::sXor(det, det_sign);

	// Check for near-zero determinant
	UVec4 det_near_zero = Vec4::sLess(abs_det, epsilon);
	Vec4 safe_det = Vec4::sSelect(det, Vec4::sOne(), det_near_zero);

	// Cramer's rule: unnormalized [t, u, v] (keep sign of det for flip)
	Vec4 t_num = b0 * c00 + b1 * c01 + b2 * c02;
	Vec4 u_num = b0 * c01 + b1 * c11 + b2 * c12;
	Vec4 v_num = b0 * c02 + b1 * c12 + b2 * c22;

	// Flip signs so that det is positive (multiply numerators by sign of det)
	Vec4 t_unscaled = Vec4::sXor(t_num, det_sign);
	Vec4 u_unscaled = Vec4::sXor(u_num, det_sign);
	Vec4 v_unscaled = Vec4::sXor(v_num, det_sign);

	// Barycentric and t checks (against abs_det since we flipped signs)
	UVec4 no_intersection =
		UVec4::sOr
		(
			UVec4::sOr
			(
				UVec4::sOr
				(
					det_near_zero,
					Vec4::sLess(t_unscaled, zero)
				),
				UVec4::sOr
				(
					Vec4::sLess(u_unscaled, zero),
					Vec4::sLess(v_unscaled, zero)
				)
			),
			Vec4::sGreater(u_unscaled + v_unscaled, abs_det)
		);

	// Compute actual t, u, v for residual check
	Vec4 inv_det = Vec4::sOne() / safe_det;
	Vec4 t = t_num * inv_det;
	Vec4 u = u_num * inv_det;
	Vec4 v = v_num * inv_det;

	// Residual: r = s + D*t - e1*u - e2*v
	Vec4 rx = sx + dx * t - e1x * u - e2x * v;
	Vec4 ry = sy + dy * t - e1y * u - e2y * v;
	Vec4 rz = sz + dz * t - e1z * u - e2z * v;
	Vec4 rw = sw + dw * t - e1w * u - e2w * v;
	Vec4 residual_sq = rx * rx + ry * ry + rz * rz + rw * rw;

	// Scale for relative residual tolerance
	Vec4 scale = Vec4::sMax(a11, Vec4::sMax(a22, a00 * t * t));
	no_intersection = UVec4::sOr(no_intersection, Vec4::sGreater(residual_sq, residual_tol * scale));

	return Vec4::sSelect(t, Vec4::sReplicate(FLT_MAX), no_intersection);
}

JPH_NAMESPACE_END
