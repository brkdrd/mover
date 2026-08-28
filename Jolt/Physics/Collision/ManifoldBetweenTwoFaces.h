// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Collision/ContactListener.h>
#include <Jolt/Physics/Collision/Shape/ConvexShape.h>

JPH_NAMESPACE_BEGIN

/// Max number of contact points that a manifold is reduced to.
///
/// TODO(4D): in 3D the contact patch is a 2D polygon and upstream reduces it to 4 points (a quad).
/// In 4D the patch lives in the 3D hyperplane perpendicular to the contact normal, so a reduced set
/// has to span 3 dimensions to resist tipping about any axis: 4 points can at best form a
/// tetrahedron and 2 of them would be redundant for a box-on-hyperplane rest pose. We keep up to 6:
/// the two endpoints of the dominant axis plus the extremes along the two remaining hyperplane
/// directions, which is the direct generalization of upstream's "line + extremes to both sides".
static constexpr int cMaxManifoldPoints = 6;

/// Remove contact points if there are > cMaxManifoldPoints of them, keeping the ones that form the
/// biggest spanning simplex of the contact patch.
/// @param inPenetrationAxis is the world space penetration axis (must be normalized)
/// @param ioContactPointsOn1 The contact points on shape 1 relative to inCenterOfMass
/// @param ioContactPointsOn2 The contact points on shape 2 relative to inCenterOfMass
/// On output ioContactPointsOn1/2 are reduced to at most cMaxManifoldPoints points.
JPH_EXPORT void PruneContactPoints(Vec4Arg inPenetrationAxis, ContactPoints &ioContactPointsOn1, ContactPoints &ioContactPointsOn2
#ifdef JPH_DEBUG_RENDERER
	, RVec4Arg inCenterOfMass
#endif // JPH_DEBUG_RENDERER
	);

/// Determine contact points between 2 faces of 2 shapes and return them in outContactPoints 1 & 2.
/// @param inContactPoint1 The contact point on shape 1 relative to inCenterOfMass
/// @param inContactPoint2 The contact point on shape 2 relative to inCenterOfMass
/// @param inPenetrationAxis The local space penetration axis in world space
/// @param inMaxContactDistance After face 2 is clipped against face 1, each remaining point on face 2 is compared to the plane of face 1. If the distance is bigger than inMaxContactDistance the point is discarded.
/// @param inShape1Face The supporting face on shape 1 relative to inCenterOfMass
/// @param inShape2Face The supporting face on shape 2 relative to inCenterOfMass
/// @param outContactPoints1 Returns the contact points between the two shapes for shape 1 relative to inCenterOfMass (any existing points in the output array are left as is)
/// @param outContactPoints2 Returns the contact points between the two shapes for shape 2 relative to inCenterOfMass (any existing points in the output array are left as is)
JPH_EXPORT void ManifoldBetweenTwoFaces(Vec4Arg inContactPoint1, Vec4Arg inContactPoint2, Vec4Arg inPenetrationAxis, float inMaxContactDistance, const ConvexShape::SupportingFace &inShape1Face, const ConvexShape::SupportingFace &inShape2Face, ContactPoints &outContactPoints1, ContactPoints &outContactPoints2
#ifdef JPH_DEBUG_RENDERER
	, RVec4Arg inCenterOfMass
#endif // JPH_DEBUG_RENDERER
	);

JPH_NAMESPACE_END
