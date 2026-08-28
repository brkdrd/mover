// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Geometry/RayCylinder.h>
#include <Jolt/Geometry/RaySphere.h>

JPH_NAMESPACE_BEGIN

/// Tests a ray against a capsule centered at the origin with the given axis direction.
/// In 4D, the capsule caps are 3-spheres (hyperspheres).
/// @param inRayOrigin Origin of the ray
/// @param inRayDirection Direction of the ray (not necessarily normalized)
/// @param inCapsuleAxis Unit vector along the capsule's axis
/// @param inCapsuleHalfHeight Distance from origin to the center of each cap sphere
/// @param inCapsuleRadius Radius of the capsule
/// @return FLT_MAX if no intersection, otherwise the fraction along the ray. Returns 0 if ray starts inside.
JPH_INLINE float RayCapsule(Vec4Arg inRayOrigin, Vec4Arg inRayDirection, Vec4Arg inCapsuleAxis, float inCapsuleHalfHeight, float inCapsuleRadius)
{
	// Test infinite cylinder
	float cylinder = RayCylinder(inRayOrigin, inRayDirection, inCapsuleAxis, inCapsuleRadius);
	if (cylinder == FLT_MAX)
		return FLT_MAX;

	// If this hit is in the finite cylinder we have our fraction
	float par = inRayOrigin.Dot(inCapsuleAxis) + cylinder * inRayDirection.Dot(inCapsuleAxis);
	if (abs(par) <= inCapsuleHalfHeight)
		return cylinder;

	// Test upper and lower sphere (4D hypersphere caps)
	Vec4 upper_center = inCapsuleHalfHeight * inCapsuleAxis;
	float upper = RaySphere(inRayOrigin, inRayDirection, upper_center, inCapsuleRadius);
	float lower = RaySphere(inRayOrigin, inRayDirection, -upper_center, inCapsuleRadius);
	return min(upper, lower);
}

/// Backward-compatible overload: Y-axis capsule
JPH_INLINE float RayCapsule(Vec4Arg inRayOrigin, Vec4Arg inRayDirection, float inCapsuleHalfHeight, float inCapsuleRadius)
{
	return RayCapsule(inRayOrigin, inRayDirection, Vec4(0, 1, 0, 0), inCapsuleHalfHeight, inCapsuleRadius);
}

JPH_NAMESPACE_END
