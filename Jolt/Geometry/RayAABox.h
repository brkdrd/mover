// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

/// Helper structure holding the reciprocal of a ray for Ray vs AABox testing
class RayInvDirection
{
public:
	/// Constructors
	inline			RayInvDirection() = default;
	inline explicit	RayInvDirection(Vec4Arg inDirection) { Set(inDirection); }

	/// Set reciprocal from ray direction
	inline void		Set(Vec4Arg inDirection)
	{
		// if (abs(inDirection) <= Epsilon) the ray is nearly parallel to the slab.
		mIsParallel = Vec4::sLessOrEqual(inDirection.Abs(), Vec4::sReplicate(1.0e-20f));

		// Calculate 1 / direction while avoiding division by zero
		mInvDirection = Vec4::sSelect(inDirection, Vec4::sOne(), mIsParallel).Reciprocal();
	}

	Vec4			mInvDirection;					///< 1 / ray direction
	UVec4			mIsParallel;					///< for each component if it is parallel to the coordinate axis
};

/// Intersect AABB with ray, returns minimal distance along ray or FLT_MAX if no hit
/// Note: Can return negative value if ray starts in box
JPH_INLINE float RayAABox(Vec4Arg inOrigin, const RayInvDirection &inInvDirection, Vec4Arg inBoundsMin, Vec4Arg inBoundsMax)
{
	// Constants
	Vec4 flt_min = Vec4::sReplicate(-FLT_MAX);
	Vec4 flt_max = Vec4::sReplicate(FLT_MAX);

	// Test against all four axes simultaneously.
	Vec4 t1 = (inBoundsMin - inOrigin) * inInvDirection.mInvDirection;
	Vec4 t2 = (inBoundsMax - inOrigin) * inInvDirection.mInvDirection;

	// Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
	// use the results from any directions parallel to the slab.
	Vec4 t_min = Vec4::sSelect(Vec4::sMin(t1, t2), flt_min, inInvDirection.mIsParallel);
	Vec4 t_max = Vec4::sSelect(Vec4::sMax(t1, t2), flt_max, inInvDirection.mIsParallel);

	// t_min = maximum of all 4 components (2-step reduction)
	t_min = Vec4::sMax(t_min, t_min.Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>());
	t_min = Vec4::sMax(t_min, t_min.Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>());

	// t_max = minimum of all 4 components (2-step reduction)
	t_max = Vec4::sMin(t_max, t_max.Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>());
	t_max = Vec4::sMin(t_max, t_max.Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>());

	// if (t_min > t_max) return FLT_MAX;
	UVec4 no_intersection = Vec4::sGreater(t_min, t_max);

	// if (t_max < 0.0f) return FLT_MAX;
	no_intersection = UVec4::sOr(no_intersection, Vec4::sLess(t_max, Vec4::sZero()));

	// if (inInvDirection.mIsParallel && !(Min <= inOrigin && inOrigin <= Max)) return FLT_MAX; else return t_min;
	UVec4 no_parallel_overlap = UVec4::sOr(Vec4::sLess(inOrigin, inBoundsMin), Vec4::sGreater(inOrigin, inBoundsMax));
	no_intersection = UVec4::sOr(no_intersection, UVec4::sAnd(inInvDirection.mIsParallel, no_parallel_overlap));
	// OR all 4 components together
	no_intersection = UVec4::sOr(no_intersection, no_intersection.Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>());
	no_intersection = UVec4::sOr(no_intersection, no_intersection.Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>());
	return Vec4::sSelect(t_min, flt_max, no_intersection).GetX();
}

/// Intersect 4 AABBs with ray, returns minimal distance along ray or FLT_MAX if no hit
/// Note: Can return negative value if ray starts in box
JPH_INLINE Vec4 RayAABox4(Vec4Arg inOrigin, const RayInvDirection &inInvDirection, Vec4Arg inBoundsMinX, Vec4Arg inBoundsMinY, Vec4Arg inBoundsMinZ, Vec4Arg inBoundsMinW, Vec4Arg inBoundsMaxX, Vec4Arg inBoundsMaxY, Vec4Arg inBoundsMaxZ, Vec4Arg inBoundsMaxW)
{
	// Constants
	Vec4 flt_min = Vec4::sReplicate(-FLT_MAX);
	Vec4 flt_max = Vec4::sReplicate(FLT_MAX);

	// Origin
	Vec4 originx = inOrigin.SplatX();
	Vec4 originy = inOrigin.SplatY();
	Vec4 originz = inOrigin.SplatZ();
	Vec4 originw = inOrigin.SplatW();

	// Parallel
	UVec4 parallelx = inInvDirection.mIsParallel.SplatX();
	UVec4 parallely = inInvDirection.mIsParallel.SplatY();
	UVec4 parallelz = inInvDirection.mIsParallel.SplatZ();
	UVec4 parallelw = inInvDirection.mIsParallel.SplatW();

	// Inverse direction
	Vec4 invdirx = inInvDirection.mInvDirection.SplatX();
	Vec4 invdiry = inInvDirection.mInvDirection.SplatY();
	Vec4 invdirz = inInvDirection.mInvDirection.SplatZ();
	Vec4 invdirw = inInvDirection.mInvDirection.SplatW();

	// Test against all four axes simultaneously.
	Vec4 t1x = (inBoundsMinX - originx) * invdirx;
	Vec4 t1y = (inBoundsMinY - originy) * invdiry;
	Vec4 t1z = (inBoundsMinZ - originz) * invdirz;
	Vec4 t1w = (inBoundsMinW - originw) * invdirw;
	Vec4 t2x = (inBoundsMaxX - originx) * invdirx;
	Vec4 t2y = (inBoundsMaxY - originy) * invdiry;
	Vec4 t2z = (inBoundsMaxZ - originz) * invdirz;
	Vec4 t2w = (inBoundsMaxW - originw) * invdirw;

	// Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
	// use the results from any directions parallel to the slab.
	Vec4 t_minx = Vec4::sSelect(Vec4::sMin(t1x, t2x), flt_min, parallelx);
	Vec4 t_miny = Vec4::sSelect(Vec4::sMin(t1y, t2y), flt_min, parallely);
	Vec4 t_minz = Vec4::sSelect(Vec4::sMin(t1z, t2z), flt_min, parallelz);
	Vec4 t_minw = Vec4::sSelect(Vec4::sMin(t1w, t2w), flt_min, parallelw);
	Vec4 t_maxx = Vec4::sSelect(Vec4::sMax(t1x, t2x), flt_max, parallelx);
	Vec4 t_maxy = Vec4::sSelect(Vec4::sMax(t1y, t2y), flt_max, parallely);
	Vec4 t_maxz = Vec4::sSelect(Vec4::sMax(t1z, t2z), flt_max, parallelz);
	Vec4 t_maxw = Vec4::sSelect(Vec4::sMax(t1w, t2w), flt_max, parallelw);

	// t_min = maximum of all 4 slab axes
	Vec4 t_min = Vec4::sMax(Vec4::sMax(t_minx, t_miny), Vec4::sMax(t_minz, t_minw));

	// t_max = minimum of all 4 slab axes
	Vec4 t_max = Vec4::sMin(Vec4::sMin(t_maxx, t_maxy), Vec4::sMin(t_maxz, t_maxw));

	// if (t_min > t_max) return FLT_MAX;
	UVec4 no_intersection = Vec4::sGreater(t_min, t_max);

	// if (t_max < 0.0f) return FLT_MAX;
	no_intersection = UVec4::sOr(no_intersection, Vec4::sLess(t_max, Vec4::sZero()));

	// if bounds are invalid return FLOAT_MAX;
	UVec4 bounds_invalid = UVec4::sOr(
		UVec4::sOr(Vec4::sGreater(inBoundsMinX, inBoundsMaxX), Vec4::sGreater(inBoundsMinY, inBoundsMaxY)),
		UVec4::sOr(Vec4::sGreater(inBoundsMinZ, inBoundsMaxZ), Vec4::sGreater(inBoundsMinW, inBoundsMaxW)));
	no_intersection = UVec4::sOr(no_intersection, bounds_invalid);

	// if (inInvDirection.mIsParallel && !(Min <= inOrigin && inOrigin <= Max)) return FLT_MAX; else return t_min;
	UVec4 no_parallel_overlapx = UVec4::sAnd(parallelx, UVec4::sOr(Vec4::sLess(originx, inBoundsMinX), Vec4::sGreater(originx, inBoundsMaxX)));
	UVec4 no_parallel_overlapy = UVec4::sAnd(parallely, UVec4::sOr(Vec4::sLess(originy, inBoundsMinY), Vec4::sGreater(originy, inBoundsMaxY)));
	UVec4 no_parallel_overlapz = UVec4::sAnd(parallelz, UVec4::sOr(Vec4::sLess(originz, inBoundsMinZ), Vec4::sGreater(originz, inBoundsMaxZ)));
	UVec4 no_parallel_overlapw = UVec4::sAnd(parallelw, UVec4::sOr(Vec4::sLess(originw, inBoundsMinW), Vec4::sGreater(originw, inBoundsMaxW)));
	no_intersection = UVec4::sOr(no_intersection, UVec4::sOr(UVec4::sOr(no_parallel_overlapx, no_parallel_overlapy), UVec4::sOr(no_parallel_overlapz, no_parallel_overlapw)));
	return Vec4::sSelect(t_min, flt_max, no_intersection);
}

/// Intersect AABB with ray, returns minimal and maximal distance along ray or FLT_MAX, -FLT_MAX if no hit
/// Note: Can return negative value for outMin if ray starts in box
JPH_INLINE void RayAABox(Vec4Arg inOrigin, const RayInvDirection &inInvDirection, Vec4Arg inBoundsMin, Vec4Arg inBoundsMax, float &outMin, float &outMax)
{
	// Constants
	Vec4 flt_min = Vec4::sReplicate(-FLT_MAX);
	Vec4 flt_max = Vec4::sReplicate(FLT_MAX);

	// Test against all four axes simultaneously.
	Vec4 t1 = (inBoundsMin - inOrigin) * inInvDirection.mInvDirection;
	Vec4 t2 = (inBoundsMax - inOrigin) * inInvDirection.mInvDirection;

	// Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
	// use the results from any directions parallel to the slab.
	Vec4 t_min = Vec4::sSelect(Vec4::sMin(t1, t2), flt_min, inInvDirection.mIsParallel);
	Vec4 t_max = Vec4::sSelect(Vec4::sMax(t1, t2), flt_max, inInvDirection.mIsParallel);

	// t_min = maximum of all 4 components (2-step reduction)
	t_min = Vec4::sMax(t_min, t_min.Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>());
	t_min = Vec4::sMax(t_min, t_min.Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>());

	// t_max = minimum of all 4 components (2-step reduction)
	t_max = Vec4::sMin(t_max, t_max.Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>());
	t_max = Vec4::sMin(t_max, t_max.Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>());

	// if (t_min > t_max) return FLT_MAX;
	UVec4 no_intersection = Vec4::sGreater(t_min, t_max);

	// if (t_max < 0.0f) return FLT_MAX;
	no_intersection = UVec4::sOr(no_intersection, Vec4::sLess(t_max, Vec4::sZero()));

	// if (inInvDirection.mIsParallel && !(Min <= inOrigin && inOrigin <= Max)) return FLT_MAX;
	UVec4 no_parallel_overlap = UVec4::sOr(Vec4::sLess(inOrigin, inBoundsMin), Vec4::sGreater(inOrigin, inBoundsMax));
	no_intersection = UVec4::sOr(no_intersection, UVec4::sAnd(inInvDirection.mIsParallel, no_parallel_overlap));
	// OR all 4 components
	no_intersection = UVec4::sOr(no_intersection, no_intersection.Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>());
	no_intersection = UVec4::sOr(no_intersection, no_intersection.Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>());
	outMin = Vec4::sSelect(t_min, flt_max, no_intersection).GetX();
	outMax = Vec4::sSelect(t_max, flt_min, no_intersection).GetX();
}

/// Intersect AABB with ray, returns true if there is a hit closer than inClosest
JPH_INLINE bool RayAABoxHits(Vec4Arg inOrigin, const RayInvDirection &inInvDirection, Vec4Arg inBoundsMin, Vec4Arg inBoundsMax, float inClosest)
{
	// Constants
	Vec4 flt_min = Vec4::sReplicate(-FLT_MAX);
	Vec4 flt_max = Vec4::sReplicate(FLT_MAX);

	// Test against all four axes simultaneously.
	Vec4 t1 = (inBoundsMin - inOrigin) * inInvDirection.mInvDirection;
	Vec4 t2 = (inBoundsMax - inOrigin) * inInvDirection.mInvDirection;

	// Compute the max of min(t1,t2) and the min of max(t1,t2)
	Vec4 t_min = Vec4::sSelect(Vec4::sMin(t1, t2), flt_min, inInvDirection.mIsParallel);
	Vec4 t_max = Vec4::sSelect(Vec4::sMax(t1, t2), flt_max, inInvDirection.mIsParallel);

	// 4-component reduction
	t_min = Vec4::sMax(t_min, t_min.Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>());
	t_min = Vec4::sMax(t_min, t_min.Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>());
	t_max = Vec4::sMin(t_max, t_max.Swizzle<SWIZZLE_Z, SWIZZLE_W, SWIZZLE_X, SWIZZLE_Y>());
	t_max = Vec4::sMin(t_max, t_max.Swizzle<SWIZZLE_Y, SWIZZLE_X, SWIZZLE_W, SWIZZLE_Z>());

	// if (t_min > t_max) return false;
	UVec4 no_intersection = Vec4::sGreater(t_min, t_max);

	// if (t_max < 0.0f) return false;
	no_intersection = UVec4::sOr(no_intersection, Vec4::sLess(t_max, Vec4::sZero()));

	// if (t_min > inClosest) return false;
	no_intersection = UVec4::sOr(no_intersection, Vec4::sGreater(t_min, Vec4::sReplicate(inClosest)));

	// if (inInvDirection.mIsParallel && !(Min <= inOrigin && inOrigin <= Max)) return false;
	UVec4 no_parallel_overlap = UVec4::sOr(Vec4::sLess(inOrigin, inBoundsMin), Vec4::sGreater(inOrigin, inBoundsMax));
	no_intersection = UVec4::sOr(no_intersection, UVec4::sAnd(inInvDirection.mIsParallel, no_parallel_overlap));

	return !no_intersection.TestAnyTrue();
}

/// Intersect AABB with ray without hit fraction, based on separating axis test.
/// In 4D, we test face normals (4 axes) and edge-face cross products (6 axes from sCross).
JPH_INLINE bool RayAABoxHits(Vec4Arg inOrigin, Vec4Arg inDirection, Vec4Arg inBoundsMin, Vec4Arg inBoundsMax)
{
	Vec4 extents = inBoundsMax - inBoundsMin;
	Vec4 diff = 2.0f * inOrigin - inBoundsMin - inBoundsMax;
	Vec4 abs_diff = diff.Abs();

	// Face normal axes: check each coordinate
	UVec4 no_intersection = UVec4::sAnd(Vec4::sGreater(abs_diff, extents), Vec4::sGreaterOrEqual(diff * inDirection, Vec4::sZero()));

	// OR all 4 face-normal results
	if (no_intersection.TestAnyTrue())
		return false;

	// Edge-face cross product axes (6 combinations from 4D: sCross of ray direction with pairs of AABox face normals)
	// For AABB, axes are e_i. The 4D cross product sCross(D, e_i, e_j) gives axes perpendicular to D and two box face normals.
	// There are C(4,2) = 6 such axes. Each reduces to testing a 2D projection.
	Vec4 abs_dir = inDirection.Abs();
	float dx = inDirection.GetX(), dy = inDirection.GetY(), dz = inDirection.GetZ(), dw = inDirection.GetW();
	float adx = abs_dir.GetX(), ady = abs_dir.GetY(), adz = abs_dir.GetZ(), adw = abs_dir.GetW();
	float ex = extents.GetX(), ey = extents.GetY(), ez = extents.GetZ(), ew = extents.GetW();
	float dfx = diff.GetX(), dfy = diff.GetY(), dfz = diff.GetZ(), dfw = diff.GetW();

	// sCross(D, e_x, e_y) -> test on Z,W plane: |D_z * diff_w - D_w * diff_z| <= e_z * |D_w| + e_w * |D_z|
	if (abs(dz * dfw - dw * dfz) > ez * adw + ew * adz) return false;
	// sCross(D, e_x, e_z) -> test on Y,W plane
	if (abs(dy * dfw - dw * dfy) > ey * adw + ew * ady) return false;
	// sCross(D, e_x, e_w) -> test on Y,Z plane
	if (abs(dy * dfz - dz * dfy) > ey * adz + ez * ady) return false;
	// sCross(D, e_y, e_z) -> test on X,W plane
	if (abs(dx * dfw - dw * dfx) > ex * adw + ew * adx) return false;
	// sCross(D, e_y, e_w) -> test on X,Z plane
	if (abs(dx * dfz - dz * dfx) > ex * adz + ez * adx) return false;
	// sCross(D, e_z, e_w) -> test on X,Y plane
	if (abs(dx * dfy - dy * dfx) > ex * ady + ey * adx) return false;

	return true;
}

JPH_NAMESPACE_END
