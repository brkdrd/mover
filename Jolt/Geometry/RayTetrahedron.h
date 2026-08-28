// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

/// Intersect ray with tetrahedron in 4D, returns closest point or FLT_MAX if no hit.
/// The tetrahedron defines a hyperplane in 4D; we compute the ray-hyperplane intersection
/// then check barycentric coordinates.
JPH_INLINE float RayTetrahedron(Vec4Arg inOrigin, Vec4Arg inDirection, Vec4Arg inV0, Vec4Arg inV1, Vec4Arg inV2, Vec4Arg inV3)
{
	// Edge vectors from V0
	Vec4 e1 = inV1 - inV0;
	Vec4 e2 = inV2 - inV0;
	Vec4 e3 = inV3 - inV0;

	// Hyperplane normal via 4D cross product
	Vec4 n = Vec4::sCross(e1, e2, e3);

	// Ray-hyperplane intersection: t = n.(V0 - O) / n.D
	float n_dot_d = n.Dot(inDirection);
	if (abs(n_dot_d) < 1.0e-12f)
		return FLT_MAX; // Ray parallel to hyperplane

	Vec4 s = inOrigin - inV0;
	float t = -n.Dot(s) / n_dot_d;
	if (t < 0.0f)
		return FLT_MAX; // Behind ray origin

	// Intersection point relative to V0
	Vec4 p = s + inDirection * t;

	// Compute barycentric coordinates via 3x3 Gram system:
	// [e1.e1  e1.e2  e1.e3] [l1]   [e1.p]
	// [e1.e2  e2.e2  e2.e3] [l2] = [e2.p]
	// [e1.e3  e2.e3  e3.e3] [l3]   [e3.p]
	float d00 = e1.Dot(e1), d01 = e1.Dot(e2), d02 = e1.Dot(e3);
	float d11 = e2.Dot(e2), d12 = e2.Dot(e3);
	float d22 = e3.Dot(e3);
	float b0 = e1.Dot(p), b1 = e2.Dot(p), b2 = e3.Dot(p);

	// Cofactors
	float c00 = d11 * d22 - d12 * d12;
	float c01 = d02 * d12 - d01 * d22;
	float c02 = d01 * d12 - d02 * d11;
	float c11 = d00 * d22 - d02 * d02;
	float c12 = d01 * d02 - d00 * d12;
	float c22 = d00 * d11 - d01 * d01;

	float det = d00 * c00 + d01 * c01 + d02 * c02;
	if (abs(det) < 1.0e-12f)
		return FLT_MAX; // Degenerate tetrahedron

	float inv_det = 1.0f / det;
	float l1 = (b0 * c00 + b1 * c01 + b2 * c02) * inv_det;
	float l2 = (b0 * c01 + b1 * c11 + b2 * c12) * inv_det;
	float l3 = (b0 * c02 + b1 * c12 + b2 * c22) * inv_det;
	float l0 = 1.0f - l1 - l2 - l3;

	// All barycentric coordinates must be >= 0
	if (l0 < 0.0f || l1 < 0.0f || l2 < 0.0f || l3 < 0.0f)
		return FLT_MAX;

	return t;
}

/// Intersect ray with 4 tetrahedra in SOA format, returns 4 vector of closest points or FLT_MAX if no hit.
JPH_INLINE Vec4 RayTetrahedron4(Vec4Arg inOrigin, Vec4Arg inDirection,
	Vec4Arg inV0X, Vec4Arg inV0Y, Vec4Arg inV0Z, Vec4Arg inV0W,
	Vec4Arg inV1X, Vec4Arg inV1Y, Vec4Arg inV1Z, Vec4Arg inV1W,
	Vec4Arg inV2X, Vec4Arg inV2Y, Vec4Arg inV2Z, Vec4Arg inV2W,
	Vec4Arg inV3X, Vec4Arg inV3Y, Vec4Arg inV3Z, Vec4Arg inV3W)
{
	Vec4 epsilon = Vec4::sReplicate(1.0e-12f);
	Vec4 zero = Vec4::sZero();

	// Edge vectors from V0 (SOA)
	Vec4 e1x = inV1X - inV0X, e1y = inV1Y - inV0Y, e1z = inV1Z - inV0Z, e1w = inV1W - inV0W;
	Vec4 e2x = inV2X - inV0X, e2y = inV2Y - inV0Y, e2z = inV2Z - inV0Z, e2w = inV2W - inV0W;
	Vec4 e3x = inV3X - inV0X, e3y = inV3Y - inV0Y, e3z = inV3Z - inV0Z, e3w = inV3W - inV0W;

	// Direction and origin-V0 splatted
	Vec4 dx = inDirection.SplatX(), dy = inDirection.SplatY(), dz = inDirection.SplatZ(), dw = inDirection.SplatW();
	Vec4 sx = inOrigin.SplatX() - inV0X, sy = inOrigin.SplatY() - inV0Y;
	Vec4 sz = inOrigin.SplatZ() - inV0Z, sw = inOrigin.SplatW() - inV0W;

	// 4D cross product n = sCross(e1, e2, e3) in SOA
	// n_i = sum of cofactors of the 3x3 minor obtained by deleting row i from [e1,e2,e3]^T
	// nx = e1y*(e2z*e3w - e2w*e3z) - e1z*(e2y*e3w - e2w*e3y) + e1w*(e2y*e3z - e2z*e3y)
	Vec4 nx =  e1y * (e2z * e3w - e2w * e3z) - e1z * (e2y * e3w - e2w * e3y) + e1w * (e2y * e3z - e2z * e3y);
	// ny = -(e1x*(e2z*e3w - e2w*e3z) - e1z*(e2x*e3w - e2w*e3x) + e1w*(e2x*e3z - e2z*e3x))
	Vec4 ny = -(e1x * (e2z * e3w - e2w * e3z) - e1z * (e2x * e3w - e2w * e3x) + e1w * (e2x * e3z - e2z * e3x));
	// nz = e1x*(e2y*e3w - e2w*e3y) - e1y*(e2x*e3w - e2w*e3x) + e1w*(e2x*e3y - e2y*e3x)
	Vec4 nz =  e1x * (e2y * e3w - e2w * e3y) - e1y * (e2x * e3w - e2w * e3x) + e1w * (e2x * e3y - e2y * e3x);
	// nw = -(e1x*(e2y*e3z - e2z*e3y) - e1y*(e2x*e3z - e2z*e3x) + e1z*(e2x*e3y - e2y*e3x))
	Vec4 nw = -(e1x * (e2y * e3z - e2z * e3y) - e1y * (e2x * e3z - e2z * e3x) + e1z * (e2x * e3y - e2y * e3x));

	// n.D and n.s
	Vec4 n_dot_d = nx * dx + ny * dy + nz * dz + nw * dw;
	Vec4 n_dot_s = nx * sx + ny * sy + nz * sz + nw * sw;

	// Check parallel rays
	Vec4 abs_n_dot_d = Vec4::sAnd(n_dot_d, UVec4::sReplicate(0x7fffffff).ReinterpretAsFloat());
	UVec4 parallel = Vec4::sLess(abs_n_dot_d, epsilon);
	Vec4 safe_n_dot_d = Vec4::sSelect(n_dot_d, Vec4::sOne(), parallel);

	// t = -n.s / n.D
	Vec4 t = -n_dot_s / safe_n_dot_d;

	// Intersection point relative to V0: p = s + D*t
	Vec4 px = sx + dx * t, py = sy + dy * t, pz = sz + dz * t, pw = sw + dw * t;

	// Gram matrix elements (dot products in SOA)
	Vec4 d00 = e1x * e1x + e1y * e1y + e1z * e1z + e1w * e1w;
	Vec4 d01 = e1x * e2x + e1y * e2y + e1z * e2z + e1w * e2w;
	Vec4 d02 = e1x * e3x + e1y * e3y + e1z * e3z + e1w * e3w;
	Vec4 d11 = e2x * e2x + e2y * e2y + e2z * e2z + e2w * e2w;
	Vec4 d12 = e2x * e3x + e2y * e3y + e2z * e3z + e2w * e3w;
	Vec4 d22 = e3x * e3x + e3y * e3y + e3z * e3z + e3w * e3w;

	// RHS
	Vec4 b0 = e1x * px + e1y * py + e1z * pz + e1w * pw;
	Vec4 b1 = e2x * px + e2y * py + e2z * pz + e2w * pw;
	Vec4 b2 = e3x * px + e3y * py + e3z * pz + e3w * pw;

	// Cofactors
	Vec4 c00 = d11 * d22 - d12 * d12;
	Vec4 c01 = d02 * d12 - d01 * d22;
	Vec4 c02 = d01 * d12 - d02 * d11;
	Vec4 c11 = d00 * d22 - d02 * d02;
	Vec4 c12 = d01 * d02 - d00 * d12;
	Vec4 c22 = d00 * d11 - d01 * d01;

	// Gram determinant
	Vec4 det = d00 * c00 + d01 * c01 + d02 * c02;

	// Make det positive, track sign
	Vec4 det_sign = Vec4::sAnd(det, UVec4::sReplicate(0x80000000).ReinterpretAsFloat());
	Vec4 abs_det = Vec4::sXor(det, det_sign);
	UVec4 det_near_zero = Vec4::sLess(abs_det, epsilon);
	// Note: no division by the determinant is needed - the barycentrics below are kept un-normalized
	// (scaled by abs_det) and only compared against zero, and det_near_zero rejects the degenerate case.

	// Barycentric numerators (Cramer's rule)
	Vec4 l1_num = b0 * c00 + b1 * c01 + b2 * c02;
	Vec4 l2_num = b0 * c01 + b1 * c11 + b2 * c12;
	Vec4 l3_num = b0 * c02 + b1 * c12 + b2 * c22;

	// Flip signs to match positive det
	Vec4 l1_u = Vec4::sXor(l1_num, det_sign);
	Vec4 l2_u = Vec4::sXor(l2_num, det_sign);
	Vec4 l3_u = Vec4::sXor(l3_num, det_sign);
	Vec4 l0_u = abs_det - l1_u - l2_u - l3_u;

	// Check: t >= 0, all lambdas >= 0
	UVec4 no_intersection =
		UVec4::sOr
		(
			UVec4::sOr
			(
				UVec4::sOr(parallel, det_near_zero),
				UVec4::sOr(Vec4::sLess(t, zero), Vec4::sLess(l0_u, zero))
			),
			UVec4::sOr
			(
				UVec4::sOr(Vec4::sLess(l1_u, zero), Vec4::sLess(l2_u, zero)),
				Vec4::sLess(l3_u, zero)
			)
		);

	return Vec4::sSelect(t, Vec4::sReplicate(FLT_MAX), no_intersection);
}

JPH_NAMESPACE_END
