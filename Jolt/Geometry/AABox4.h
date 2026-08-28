// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Geometry/OrientedBox.h>

JPH_NAMESPACE_BEGIN

/// Helper functions that process 4 axis aligned boxes at the same time using SIMD (4D version)
/// Test if 4 bounding boxes overlap with 1 bounding box, splat 1 box
JPH_INLINE UVec4 AABox4VsBox(const AABox &inBox1, Vec4Arg inBox2MinX, Vec4Arg inBox2MinY, Vec4Arg inBox2MinZ, Vec4Arg inBox2MinW, Vec4Arg inBox2MaxX, Vec4Arg inBox2MaxY, Vec4Arg inBox2MaxZ, Vec4Arg inBox2MaxW)
{
	// Splat values of box 1
	Vec4 box1_minx = inBox1.mMin.SplatX();
	Vec4 box1_miny = inBox1.mMin.SplatY();
	Vec4 box1_minz = inBox1.mMin.SplatZ();
	Vec4 box1_minw = inBox1.mMin.SplatW();
	Vec4 box1_maxx = inBox1.mMax.SplatX();
	Vec4 box1_maxy = inBox1.mMax.SplatY();
	Vec4 box1_maxz = inBox1.mMax.SplatZ();
	Vec4 box1_maxw = inBox1.mMax.SplatW();

	// Test separation over each axis
	UVec4 nooverlapx = UVec4::sOr(Vec4::sGreater(box1_minx, inBox2MaxX), Vec4::sGreater(inBox2MinX, box1_maxx));
	UVec4 nooverlapy = UVec4::sOr(Vec4::sGreater(box1_miny, inBox2MaxY), Vec4::sGreater(inBox2MinY, box1_maxy));
	UVec4 nooverlapz = UVec4::sOr(Vec4::sGreater(box1_minz, inBox2MaxZ), Vec4::sGreater(inBox2MinZ, box1_maxz));
	UVec4 nooverlapw = UVec4::sOr(Vec4::sGreater(box1_minw, inBox2MaxW), Vec4::sGreater(inBox2MinW, box1_maxw));

	// Return overlap
	return UVec4::sNot(UVec4::sOr(UVec4::sOr(nooverlapx, nooverlapy), UVec4::sOr(nooverlapz, nooverlapw)));
}

/// Test if 4 axis aligned boxes overlap with an oriented box.
/// TODO(4D): this is a conservative test - it collides the world-space AABB of the oriented box against
/// the 4 boxes (never misses a real overlap, but may report false positives). Tighten with a full 4D
/// OBB-vs-AABB SAT if broad-phase quality requires it.
JPH_INLINE UVec4 AABox4VsBox(const OrientedBox &inBox1, Vec4Arg inBox2MinX, Vec4Arg inBox2MinY, Vec4Arg inBox2MinZ, Vec4Arg inBox2MinW, Vec4Arg inBox2MaxX, Vec4Arg inBox2MaxY, Vec4Arg inBox2MaxZ, Vec4Arg inBox2MaxW)
{
	// World-space extent of the oriented box: extent_i = sum_j |R(i,j)| * halfExtents_j
	Mat44 abs_orient;
	for (uint i = 0; i < 4; ++i)
		abs_orient.SetColumn4(i, inBox1.mOrientation.GetColumn4(i).Abs());
	Vec4 extent = abs_orient * inBox1.mHalfExtents;

	AABox aabb(inBox1.mTranslation - extent, inBox1.mTranslation + extent);
	return AABox4VsBox(aabb, inBox2MinX, inBox2MinY, inBox2MinZ, inBox2MinW, inBox2MaxX, inBox2MaxY, inBox2MaxZ, inBox2MaxW);
}

/// Scale 4 axis aligned boxes
JPH_INLINE void AABox4Scale(Vec4Arg inScale, Vec4Arg inBoxMinX, Vec4Arg inBoxMinY, Vec4Arg inBoxMinZ, Vec4Arg inBoxMinW, Vec4Arg inBoxMaxX, Vec4Arg inBoxMaxY, Vec4Arg inBoxMaxZ, Vec4Arg inBoxMaxW, Vec4 &outBoundsMinX, Vec4 &outBoundsMinY, Vec4 &outBoundsMinZ, Vec4 &outBoundsMinW, Vec4 &outBoundsMaxX, Vec4 &outBoundsMaxY, Vec4 &outBoundsMaxZ, Vec4 &outBoundsMaxW)
{
	Vec4 scale_x = inScale.SplatX();
	Vec4 scaled_min_x = scale_x * inBoxMinX;
	Vec4 scaled_max_x = scale_x * inBoxMaxX;
	outBoundsMinX = Vec4::sMin(scaled_min_x, scaled_max_x);
	outBoundsMaxX = Vec4::sMax(scaled_min_x, scaled_max_x);

	Vec4 scale_y = inScale.SplatY();
	Vec4 scaled_min_y = scale_y * inBoxMinY;
	Vec4 scaled_max_y = scale_y * inBoxMaxY;
	outBoundsMinY = Vec4::sMin(scaled_min_y, scaled_max_y);
	outBoundsMaxY = Vec4::sMax(scaled_min_y, scaled_max_y);

	Vec4 scale_z = inScale.SplatZ();
	Vec4 scaled_min_z = scale_z * inBoxMinZ;
	Vec4 scaled_max_z = scale_z * inBoxMaxZ;
	outBoundsMinZ = Vec4::sMin(scaled_min_z, scaled_max_z);
	outBoundsMaxZ = Vec4::sMax(scaled_min_z, scaled_max_z);

	Vec4 scale_w = inScale.SplatW();
	Vec4 scaled_min_w = scale_w * inBoxMinW;
	Vec4 scaled_max_w = scale_w * inBoxMaxW;
	outBoundsMinW = Vec4::sMin(scaled_min_w, scaled_max_w);
	outBoundsMaxW = Vec4::sMax(scaled_min_w, scaled_max_w);
}

/// Enlarge 4 bounding boxes with extent (add to both sides)
JPH_INLINE void AABox4EnlargeWithExtent(Vec4Arg inExtent, Vec4 &ioBoundsMinX, Vec4 &ioBoundsMinY, Vec4 &ioBoundsMinZ, Vec4 &ioBoundsMinW, Vec4 &ioBoundsMaxX, Vec4 &ioBoundsMaxY, Vec4 &ioBoundsMaxZ, Vec4 &ioBoundsMaxW)
{
	Vec4 extent_x = inExtent.SplatX();
	ioBoundsMinX -= extent_x;
	ioBoundsMaxX += extent_x;

	Vec4 extent_y = inExtent.SplatY();
	ioBoundsMinY -= extent_y;
	ioBoundsMaxY += extent_y;

	Vec4 extent_z = inExtent.SplatZ();
	ioBoundsMinZ -= extent_z;
	ioBoundsMaxZ += extent_z;

	Vec4 extent_w = inExtent.SplatW();
	ioBoundsMinW -= extent_w;
	ioBoundsMaxW += extent_w;
}

/// Test if 4 bounding boxes overlap with a point
JPH_INLINE UVec4 AABox4VsPoint(Vec4Arg inPoint, Vec4Arg inBoxMinX, Vec4Arg inBoxMinY, Vec4Arg inBoxMinZ, Vec4Arg inBoxMinW, Vec4Arg inBoxMaxX, Vec4Arg inBoxMaxY, Vec4Arg inBoxMaxZ, Vec4Arg inBoxMaxW)
{
	// Splat point to 4 component vectors
	Vec4 point_x = inPoint.SplatX();
	Vec4 point_y = inPoint.SplatY();
	Vec4 point_z = inPoint.SplatZ();
	Vec4 point_w = inPoint.SplatW();

	// Test if point overlaps with box
	UVec4 overlapx = UVec4::sAnd(Vec4::sGreaterOrEqual(point_x, inBoxMinX), Vec4::sLessOrEqual(point_x, inBoxMaxX));
	UVec4 overlapy = UVec4::sAnd(Vec4::sGreaterOrEqual(point_y, inBoxMinY), Vec4::sLessOrEqual(point_y, inBoxMaxY));
	UVec4 overlapz = UVec4::sAnd(Vec4::sGreaterOrEqual(point_z, inBoxMinZ), Vec4::sLessOrEqual(point_z, inBoxMaxZ));
	UVec4 overlapw = UVec4::sAnd(Vec4::sGreaterOrEqual(point_w, inBoxMinW), Vec4::sLessOrEqual(point_w, inBoxMaxW));

	// Test if all are overlapping
	return UVec4::sAnd(UVec4::sAnd(overlapx, overlapy), UVec4::sAnd(overlapz, overlapw));
}

/// Get the squared distance between 4 AABoxes and a point
JPH_INLINE Vec4 AABox4DistanceSqToPoint(Vec4Arg inPointX, Vec4Arg inPointY, Vec4Arg inPointZ, Vec4Arg inPointW, Vec4Arg inBoxMinX, Vec4Arg inBoxMinY, Vec4Arg inBoxMinZ, Vec4Arg inBoxMinW, Vec4Arg inBoxMaxX, Vec4Arg inBoxMaxY, Vec4Arg inBoxMaxZ, Vec4Arg inBoxMaxW)
{
	// Get closest point on box
	Vec4 closest_x = Vec4::sMin(Vec4::sMax(inPointX, inBoxMinX), inBoxMaxX);
	Vec4 closest_y = Vec4::sMin(Vec4::sMax(inPointY, inBoxMinY), inBoxMaxY);
	Vec4 closest_z = Vec4::sMin(Vec4::sMax(inPointZ, inBoxMinZ), inBoxMaxZ);
	Vec4 closest_w = Vec4::sMin(Vec4::sMax(inPointW, inBoxMinW), inBoxMaxW);

	// Return the squared distance between the box and point
	return Square(closest_x - inPointX) + Square(closest_y - inPointY) + Square(closest_z - inPointZ) + Square(closest_w - inPointW);
}

/// Get the squared distance between 4 AABoxes and a point
JPH_INLINE Vec4 AABox4DistanceSqToPoint(Vec4Arg inPoint, Vec4Arg inBoxMinX, Vec4Arg inBoxMinY, Vec4Arg inBoxMinZ, Vec4Arg inBoxMinW, Vec4Arg inBoxMaxX, Vec4Arg inBoxMaxY, Vec4Arg inBoxMaxZ, Vec4Arg inBoxMaxW)
{
	return AABox4DistanceSqToPoint(inPoint.SplatX(), inPoint.SplatY(), inPoint.SplatZ(), inPoint.SplatW(), inBoxMinX, inBoxMinY, inBoxMinZ, inBoxMinW, inBoxMaxX, inBoxMaxY, inBoxMaxZ, inBoxMaxW);
}

/// Test 4 AABoxes vs a sphere
JPH_INLINE UVec4 AABox4VsSphere(Vec4Arg inCenterX, Vec4Arg inCenterY, Vec4Arg inCenterZ, Vec4Arg inCenterW, Vec4Arg inRadiusSq, Vec4Arg inBoxMinX, Vec4Arg inBoxMinY, Vec4Arg inBoxMinZ, Vec4Arg inBoxMinW, Vec4Arg inBoxMaxX, Vec4Arg inBoxMaxY, Vec4Arg inBoxMaxZ, Vec4Arg inBoxMaxW)
{
	// Test the distance from the center of the sphere to the box is smaller than the radius
	Vec4 distance_sq = AABox4DistanceSqToPoint(inCenterX, inCenterY, inCenterZ, inCenterW, inBoxMinX, inBoxMinY, inBoxMinZ, inBoxMinW, inBoxMaxX, inBoxMaxY, inBoxMaxZ, inBoxMaxW);
	return Vec4::sLessOrEqual(distance_sq, inRadiusSq);
}

/// Test 4 AABoxes vs a sphere
JPH_INLINE UVec4 AABox4VsSphere(Vec4Arg inCenter, float inRadiusSq, Vec4Arg inBoxMinX, Vec4Arg inBoxMinY, Vec4Arg inBoxMinZ, Vec4Arg inBoxMinW, Vec4Arg inBoxMaxX, Vec4Arg inBoxMaxY, Vec4Arg inBoxMaxZ, Vec4Arg inBoxMaxW)
{
	return AABox4VsSphere(inCenter.SplatX(), inCenter.SplatY(), inCenter.SplatZ(), inCenter.SplatW(), Vec4::sReplicate(inRadiusSq), inBoxMinX, inBoxMinY, inBoxMinZ, inBoxMinW, inBoxMaxX, inBoxMaxY, inBoxMaxZ, inBoxMaxW);
}

JPH_NAMESPACE_END
