// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/FindRoot.h>

JPH_NAMESPACE_BEGIN

/// Tests a ray against an infinite cylinder centered at the origin with the given axis direction.
/// @param inRayOrigin Origin of the ray
/// @param inRayDirection Direction of the ray (not necessarily normalized)
/// @param inCylinderAxis Unit vector along the cylinder's axis
/// @param inCylinderRadius Radius of the infinite cylinder
/// @return FLT_MAX if no intersection, otherwise the fraction along the ray. Returns 0 if ray starts inside.
JPH_INLINE float RayCylinder(Vec4Arg inRayOrigin, Vec4Arg inRayDirection, Vec4Arg inCylinderAxis, float inCylinderRadius)
{
	// Decompose ray origin and direction into parallel and perpendicular to axis
	float origin_par = inRayOrigin.Dot(inCylinderAxis);
	Vec4 origin_perp = inRayOrigin - origin_par * inCylinderAxis;
	float origin_perp_len_sq = origin_perp.LengthSq();
	float r_sq = Square(inCylinderRadius);

	if (origin_perp_len_sq > r_sq)
	{
		// Ray starts outside the infinite cylinder
		float dir_par = inRayDirection.Dot(inCylinderAxis);
		Vec4 direction_perp = inRayDirection - dir_par * inCylinderAxis;
		float a = direction_perp.LengthSq();
		float b = 2.0f * origin_perp.Dot(direction_perp);
		float c = origin_perp_len_sq - r_sq;
		float fraction1, fraction2;
		if (FindRoot(a, b, c, fraction1, fraction2) == 0)
			return FLT_MAX;

		float fraction = min(fraction1, fraction2);
		if (fraction >= 0.0f)
			return fraction;
	}
	else
	{
		// Ray starts inside the infinite cylinder
		return 0.0f;
	}

	return FLT_MAX;
}

/// Backward-compatible overload: Y-axis cylinder
JPH_INLINE float RayCylinder(Vec4Arg inRayOrigin, Vec4Arg inRayDirection, float inCylinderRadius)
{
	return RayCylinder(inRayOrigin, inRayDirection, Vec4(0, 1, 0, 0), inCylinderRadius);
}

/// Tests a ray against a finite cylinder centered at the origin with the given axis and half height.
/// @param inRayOrigin Origin of the ray
/// @param inRayDirection Direction of the ray (not necessarily normalized)
/// @param inCylinderAxis Unit vector along the cylinder's axis
/// @param inCylinderHalfHeight Distance from origin to top/bottom cap
/// @param inCylinderRadius Radius of the cylinder
/// @return FLT_MAX if no intersection, otherwise the fraction along the ray. Returns 0 if ray starts inside.
JPH_INLINE float RayCylinder(Vec4Arg inRayOrigin, Vec4Arg inRayDirection, Vec4Arg inCylinderAxis, float inCylinderHalfHeight, float inCylinderRadius)
{
	// Test infinite cylinder
	float fraction = RayCylinder(inRayOrigin, inRayDirection, inCylinderAxis, inCylinderRadius);
	if (fraction == FLT_MAX)
		return FLT_MAX;

	// If this hit is in the finite cylinder we have our fraction
	float origin_par = inRayOrigin.Dot(inCylinderAxis);
	float dir_par = inRayDirection.Dot(inCylinderAxis);
	if (abs(origin_par + fraction * dir_par) <= inCylinderHalfHeight)
		return fraction;

	// Check if ray could hit the top or bottom cap plane
	if (dir_par != 0.0f)
	{
		float plane_fraction;
		if (dir_par < 0.0f)
			plane_fraction = (inCylinderHalfHeight - origin_par) / dir_par;
		else
			plane_fraction = -(inCylinderHalfHeight + origin_par) / dir_par;

		if (plane_fraction >= 0.0f)
		{
			// Test if hit point is inside the cylinder radius (perpendicular distance from axis)
			Vec4 point = inRayOrigin + plane_fraction * inRayDirection;
			float point_par = point.Dot(inCylinderAxis);
			Vec4 point_perp = point - point_par * inCylinderAxis;
			if (point_perp.LengthSq() <= Square(inCylinderRadius))
				return plane_fraction;
		}
	}

	return FLT_MAX;
}

/// Backward-compatible overload: Y-axis finite cylinder
JPH_INLINE float RayCylinder(Vec4Arg inRayOrigin, Vec4Arg inRayDirection, float inCylinderHalfHeight, float inCylinderRadius)
{
	return RayCylinder(inRayOrigin, inRayDirection, Vec4(0, 1, 0, 0), inCylinderHalfHeight, inCylinderRadius);
}

JPH_NAMESPACE_END
