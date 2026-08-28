// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Collision/ManifoldBetweenTwoFaces.h>
#include <Jolt/Physics/Constraints/ContactConstraintManager.h>

JPH_NAMESPACE_BEGIN

/// Build an orthonormal basis {outDirA, outDirB} of the 2D subspace that is perpendicular to both the
/// contact normal and inEdge, i.e. what is left of the 3D contact hyperplane after the dominant axis
/// of the patch has been taken out. Deterministic for a given input.
static void sGetRemainingHyperplaneBasis(Vec4Arg inNormal, Vec4Arg inEdge, Vec4 &outDirA, Vec4 &outDirB)
{
	// Pick the coordinate axis least aligned with the two directions we already have so the
	// generalized cross product below is well conditioned
	Vec4 align = inNormal.Abs() + inEdge.Abs();
	int best = 0;
	float best_v = align.GetX();
	if (align.GetY() < best_v) { best = 1; best_v = align.GetY(); }
	if (align.GetZ() < best_v) { best = 2; best_v = align.GetZ(); }
	if (align.GetW() < best_v) { best = 3; }
	Vec4 seed = best == 0? Vec4::sAxisX() : best == 1? Vec4::sAxisY() : best == 2? Vec4::sAxisZ() : Vec4::sAxisW();

	// sCross(a, b, c) is perpendicular to all three of its arguments
	outDirA = Vec4::sCross(inNormal, inEdge, seed).Normalized();
	outDirB = Vec4::sCross(inNormal, inEdge, outDirA).Normalized();
}

void PruneContactPoints(Vec4Arg inPenetrationAxis, ContactPoints &ioContactPointsOn1, ContactPoints &ioContactPointsOn2 JPH_IF_DEBUG_RENDERER(, RVec4Arg inCenterOfMass))
{
	// Makes no sense to call this with cMaxManifoldPoints or fewer points
	JPH_ASSERT(ioContactPointsOn1.size() > cMaxManifoldPoints);

	// Both arrays should have the same size
	JPH_ASSERT(ioContactPointsOn1.size() == ioContactPointsOn2.size());

	// Penetration axis must be normalized
	JPH_ASSERT(inPenetrationAxis.IsNormalized());

	// We use a heuristic of (distance to center of mass) * (penetration depth) to find the contact point that we should keep
	// Neither of those two terms should ever become zero, so we clamp against this minimum value
	constexpr float cMinDistanceSq = 1.0e-6f; // 1 mm

	ContactPoints projected;
	StaticArray<float, 64> penetration_depth_sq;
	for (ContactPoints::size_type i = 0; i < ioContactPointsOn1.size(); ++i)
	{
		// Project the contact points into the hyperplane through inCenterOfMass with normal
		// inPenetrationAxis (all points are relative to inCenterOfMass, so project through the origin)
		Vec4 v1 = ioContactPointsOn1[i];
		projected.push_back(v1 - v1.Dot(inPenetrationAxis) * inPenetrationAxis);

		// Calculate penetration depth^2 of each point and clamp against the minimal distance
		Vec4 v2 = ioContactPointsOn2[i];
		penetration_depth_sq.push_back(max(cMinDistanceSq, (v2 - v1).LengthSq()));
	}

	// Find the point that is furthest away from the center of mass (its torque will have the biggest influence)
	// and the point that has the deepest penetration depth. Use the heuristic (distance to center of mass) * (penetration depth) for this.
	uint point1 = 0;
	float val = max(cMinDistanceSq, projected[0].LengthSq()) * penetration_depth_sq[0];
	for (uint i = 0; i < projected.size(); ++i)
	{
		float v = max(cMinDistanceSq, projected[i].LengthSq()) * penetration_depth_sq[i];
		if (v > val)
		{
			val = v;
			point1 = i;
		}
	}
	Vec4 point1v = projected[point1];

	// Find point furthest from the first point forming a line segment with point1. Again combine this with the heuristic
	// for deepest point as per above.
	uint point2 = uint(-1);
	val = -FLT_MAX;
	for (uint i = 0; i < projected.size(); ++i)
		if (i != point1)
		{
			float v = max(cMinDistanceSq, (projected[i] - point1v).LengthSq()) * penetration_depth_sq[i];
			if (v > val)
			{
				val = v;
				point2 = i;
			}
		}
	JPH_ASSERT(point2 != uint(-1));
	Vec4 point2v = projected[point2];

	// TODO(4D): in 3D what is left of the contact plane after taking out the point1..point2 axis is a
	// single direction, and upstream picks the extreme point to either side of it. In 4D the patch is
	// 3 dimensional, so what is left is a 2D subspace of the contact hyperplane: we take the extremes
	// along both of its basis directions, giving up to 6 points whose hull spans the patch.
	uint point3 = uint(-1), point4 = uint(-1), point5 = uint(-1), point6 = uint(-1);
	Vec4 edge = point2v - point1v;
	float edge_len_sq = edge.LengthSq();
	if (edge_len_sq > cMinDistanceSq)
	{
		Vec4 dir_a, dir_b;
		sGetRemainingHyperplaneBasis(inPenetrationAxis, edge / sqrt(edge_len_sq), dir_a, dir_b);

		float min_a = 0.0f, max_a = 0.0f, min_b = 0.0f, max_b = 0.0f;
		for (uint i = 0; i < projected.size(); ++i)
			if (i != point1 && i != point2)
			{
				Vec4 rel = projected[i] - point1v;

				float va = dir_a.Dot(rel);
				if (va < min_a) { min_a = va; point3 = i; }
				else if (va > max_a) { max_a = va; point4 = i; }

				float vb = dir_b.Dot(rel);
				if (vb < min_b) { min_b = vb; point5 = i; }
				else if (vb > max_b) { max_b = vb; point6 = i; }
			}
	}

	// Add points to array. point1 and point2 are the endpoints of the dominant axis, the rest are the
	// extremes of the two remaining hyperplane directions; duplicates are skipped (a single point can
	// be the extreme along more than one direction).
	StaticArray<Vec4, cMaxManifoldPoints> points_to_keep_on_1, points_to_keep_on_2;
	uint candidates[] = { point1, point2, point3, point4, point5, point6 };
	for (uint i = 0; i < std::size(candidates); ++i)
	{
		uint c = candidates[i];
		if (c == uint(-1))
			continue;

		// A single point can be the extreme along more than one direction, so skip indices we already took
		bool already_added = false;
		for (uint j = 0; j < i; ++j)
			if (candidates[j] == c)
			{
				already_added = true;
				break;
			}

		if (!already_added)
		{
			points_to_keep_on_1.push_back(ioContactPointsOn1[c]);
			points_to_keep_on_2.push_back(ioContactPointsOn2[c]);
		}
	}

	// Copy the points back to the input buffer
	ioContactPointsOn1 = points_to_keep_on_1;
	ioContactPointsOn2 = points_to_keep_on_2;
}

void ManifoldBetweenTwoFaces(Vec4Arg inContactPoint1, Vec4Arg inContactPoint2, Vec4Arg inPenetrationAxis, float inMaxContactDistance, const ConvexShape::SupportingFace &inShape1Face, const ConvexShape::SupportingFace &inShape2Face, ContactPoints &outContactPoints1, ContactPoints &outContactPoints2 JPH_IF_DEBUG_RENDERER(, RVec4Arg inCenterOfMass))
{
	JPH_ASSERT(inMaxContactDistance > 0.0f);

	// Remember size before adding new points, to check at the end if we added some
	ContactPoints::size_type old_size = outContactPoints1.size();

	// TODO(4D): THIS IS THE REDUCED v1 PATH -- conceptual gap 1 in CLAUDE.md.
	//
	// In 3D a supporting face is a polygon and upstream clips face 2 against face 1 with
	// Sutherland-Hodgman (ClipPolyVsPoly / ClipPolyVsEdge in Geometry/ClipPoly.h) to produce a contact
	// patch of up to 32 points, which PruneContactPoints then reduces.
	//
	// In 4D a supporting face is a 3D *polyhedron* living in the contact hyperplane (BoxShape returns a
	// cubic cell, PlaneShape a 3-cube patch), and ConvexShape::SupportingFace is only an unordered
	// vertex list -- clipping a polyhedron against a polyhedron needs face connectivity, which the
	// support function does not currently return. Until that plumbing plus a halfspace clip for
	// polyhedra exists (ClipPoly.h is still 2D), we fall back to the single deepest contact point,
	// which upstream also does whenever clipping produces nothing.
	//
	// Consequence: collisions, restitution and friction are correct, but a resting box gets one contact
	// point per manifold instead of a patch, so large flat stacks will be less stable than in 3D Jolt.
	// The reduction step above is already the real 4D algorithm and is ready for the clipped patch.
	(void)inPenetrationAxis;
	(void)inMaxContactDistance;
	(void)inShape1Face;
	(void)inShape2Face;

	// If the clipping result is empty, use the contact point itself
	if (outContactPoints1.size() == old_size)
	{
		outContactPoints1.push_back(inContactPoint1);
		outContactPoints2.push_back(inContactPoint2);
	}
}

JPH_NAMESPACE_END
