// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Geometry/AABox.h>

JPH_NAMESPACE_BEGIN

// TODO(4D): everything below the "4D polytope clipping" block is 3D Jolt's polygon clipping. It
// assumes the vertices form an ordered 2D loop (it walks consecutive pairs as edges), which a 4D
// supporting face does not: that face is a 3-dimensional polytope and GetSupportingFace returns its
// vertices unordered. Those functions are unused and kept only for reference while the shapes that
// needed them (Mesh/HeightField/Triangle) are cut.

/// A halfspace used to clip a contact patch, in the form mNormal . x <= mDistance.
///
/// The normal is perpendicular to the penetration axis, so the halfspace is a "prism" that extends
/// along the axis - the 4D counterpart of 3D Jolt building its clip planes from an edge and the
/// penetration axis. Testing a point against it is therefore automatically a test of that point
/// *projected along the axis*, which is what contact clipping wants.
struct ClipHalfSpace
{
	Vec4			mNormal;								///< Outward normal, normalized, perpendicular to the penetration axis
	float			mDistance;								///< Inside is mNormal . x <= mDistance
};

/// Upper bounds for the scratch tables in ClipPolytopeVsPolytope. A supporting face in the v1 shape
/// set has at most 8 vertices (a tesseract's cubic cell), so these are generous.
static constexpr int cMaxClipHalfSpaces = 64;
static constexpr int cMaxClipVertices = 32;

using ClipHalfSpaces = StaticArray<ClipHalfSpace, cMaxClipHalfSpaces>;

/// Derive the supporting halfspaces of the convex hull of inVertices, within the hyperplane
/// perpendicular to inAxis.
///
/// TODO(4D): a 4D supporting face is a 3D polytope given as an unordered vertex set, so there is no
/// edge loop to walk as in 3D. We enumerate candidate facets over vertex triples and keep a triple
/// when every other vertex lies on one side of it. That is O(n^3 * n), which is fine here because n
/// is small by construction: a tesseract's or plane patch's supporting cell is a cube (8 vertices)
/// and a convex hull's is a tetrahedral cell (4), against a hard cap of 32.
template <class VERTEX_ARRAY>
void GetClipHalfSpaces(const VERTEX_ARRAY &inVertices, Vec4Arg inAxis, ClipHalfSpaces &outHalfSpaces)
{
	outHalfSpaces.clear();

	int n = int(inVertices.size());
	if (n < 3)
		return;

	// Scale the tolerances to the size of the face so this is resolution independent
	Vec4 centroid = Vec4::sZero();
	for (int i = 0; i < n; ++i)
		centroid += inVertices[i];
	centroid = centroid / float(n);
	float extent_sq = 0.0f;
	for (int i = 0; i < n; ++i)
		extent_sq = max(extent_sq, (inVertices[i] - centroid).LengthSq());
	float extent = sqrt(extent_sq);
	if (extent < FLT_EPSILON)
		return;
	float dist_tolerance = 1.0e-4f * extent;

	for (int i = 0; i < n; ++i)
		for (int j = i + 1; j < n; ++j)
			for (int k = j + 1; k < n; ++k)
			{
				Vec4 origin = inVertices[i];

				// Normal of the facet within the hyperplane: perpendicular to both facet directions
				// and to the axis
				Vec4 normal = Vec4::sCross(inVertices[j] - origin, inVertices[k] - origin, inAxis);
				float normal_len_sq = normal.LengthSq();
				if (normal_len_sq < FLT_EPSILON)
					continue; // Degenerate triple (collinear, or the facet contains the axis)

				// Discard a plane we already have before doing anything expensive with it. Many
				// triples describe the same facet - a cube cell has 8 vertices and C(8,3) = 56
				// triples but at most 12 facets - so both the normalization (a square root) and the
				// support scan below are wasted on most of them. The stored normals are unit, so a
				// candidate is parallel to one when dot^2 is |candidate|^2, and then the plane's
				// distance can be read off the stored normal without normalizing the candidate.
				bool duplicate = false;
				for (const ClipHalfSpace &h : outHalfSpaces)
				{
					float dot = h.mNormal.Dot(normal);
					if (dot * dot > 0.998f * normal_len_sq && abs(h.mNormal.Dot(origin) - h.mDistance) < dist_tolerance)
					{
						duplicate = true;
						break;
					}
				}
				if (duplicate)
					continue;

				normal = normal / sqrt(normal_len_sq);
				float d = normal.Dot(origin);

				// Keep the triple only if it supports the whole set
				bool all_below = true, all_above = true;
				for (int m = 0; m < n; ++m)
				{
					float dm = normal.Dot(inVertices[m]) - d;
					if (dm > dist_tolerance)
						all_below = false;
					if (dm < -dist_tolerance)
						all_above = false;
					if (!all_below && !all_above)
						break;
				}
				if (!all_below && !all_above)
					continue;
				if (!all_below)
				{
					normal = -normal;
					d = -d;
				}

				if (outHalfSpaces.size() < outHalfSpaces.capacity())
					outHalfSpaces.push_back({ normal, d });
			}
}

/// Check if inPoint is inside all of inHalfSpaces (within inTolerance)
inline bool IsInsideClipHalfSpaces(Vec4Arg inPoint, const ClipHalfSpaces &inHalfSpaces, float inTolerance)
{
	for (const ClipHalfSpace &h : inHalfSpaces)
		if (h.mNormal.Dot(inPoint) - h.mDistance > inTolerance)
			return false;
	return true;
}

/// Clip the segment inA -> inB against inHalfSpaces, returning the surviving parameter range.
/// Returns false if the segment lies entirely outside.
inline bool ClipSegmentVsHalfSpaces(Vec4Arg inA, Vec4Arg inB, const ClipHalfSpaces &inHalfSpaces, float inTolerance, float &outT0, float &outT1)
{
	outT0 = 0.0f;
	outT1 = 1.0f;
	Vec4 delta = inB - inA;

	for (const ClipHalfSpace &h : inHalfSpaces)
	{
		float da = h.mNormal.Dot(inA) - h.mDistance;
		float dd = h.mNormal.Dot(delta);
		if (abs(dd) < FLT_EPSILON)
		{
			// Parallel to the halfspace boundary
			if (da > inTolerance)
				return false;
			continue;
		}
		float t = -da / dd;
		if (dd > 0.0f)
			outT1 = min(outT1, t);		// Leaving the halfspace
		else
			outT0 = max(outT0, t);		// Entering it
		if (outT0 > outT1)
			return false;
	}
	return true;
}

/// Clip the convex polytope inPolytopeToClip against the convex polytope inClippingPolytope, both
/// given as unordered vertex sets of a 4D supporting face, along inAxis.
///
/// TODO(4D): this replaces 3D Jolt's Sutherland-Hodgman polygon clip. In 3D the contact patch is a
/// 2D polygon and you can clip an ordered loop edge by edge; in 4D it is a 3D polyhedron and there
/// is no loop. Instead we build the intersection from its definition: the vertices of conv(P) n
/// conv(Q) are exactly the vertices of P inside Q, the vertices of Q inside P, and the points where
/// an edge of one crosses a facet of the other. The crossings come out of clipping segments against
/// the other polytope's halfspaces.
///
/// We feed in ALL vertex pairs rather than computing the true edge set. A pair that is not an edge
/// is a chord lying inside the polytope, so clipping it can only produce points that are already
/// inside the intersection - harmless, because the caller reduces the result with
/// PruneContactPoints, which selects extreme points. The output is therefore a superset of the
/// intersection's vertices, all of which lie in the intersection, rather than exactly its vertex
/// set. Measured output sizes for typical v1 faces: a tesseract resting on a plane patch gives
/// exactly its 8 corners, a half-overlapping pair of cells gives 18, against a cap of 64.
///
/// Output points lie on the affine hull of inPolytopeToClip (mirroring 3D Jolt, where the clipped
/// points lie on face 2 and the caller then projects them onto the plane of face 1).
template <class VERTEX_ARRAY, class OUT_ARRAY>
void ClipPolytopeVsPolytope(const VERTEX_ARRAY &inPolytopeToClip, const VERTEX_ARRAY &inClippingPolytope, Vec4Arg inAxis, OUT_ARRAY &outClippedPolytope)
{
	outClippedPolytope.clear();

	int n_clip = int(inPolytopeToClip.size());
	int n_clipper = int(inClippingPolytope.size());
	if (n_clip == 0 || n_clipper == 0)
		return;

	ClipHalfSpaces clipper_planes;
	GetClipHalfSpaces(inClippingPolytope, inAxis, clipper_planes);
	if (clipper_planes.empty())
		return;

	// Tolerance relative to the size of the face being clipped
	Vec4 centroid = Vec4::sZero();
	for (int i = 0; i < n_clip; ++i)
		centroid += inPolytopeToClip[i];
	centroid = centroid / float(n_clip);
	float extent = 0.0f;
	for (int i = 0; i < n_clip; ++i)
		extent = max(extent, (inPolytopeToClip[i] - centroid).Length());
	float tolerance = 1.0e-4f * max(extent, FLT_EPSILON);

	// Fast path: the polytope is entirely inside the clipper, so the intersection is the polytope
	// itself. This is the common case in a resting stack - a box face landing inside a larger face,
	// or two equal boxes flush - and it skips the clipper's facet enumeration for the second
	// polytope plus every segment test. It is exact, not an approximation: if P is contained in Q
	// then no vertex of Q lies strictly inside P (a vertex of Q is extreme in Q, so it cannot be
	// interior to a subset of it) and no segment of P crosses the boundary of Q.
	bool all_inside = true;
	for (int i = 0; i < n_clip && all_inside; ++i)
		all_inside = IsInsideClipHalfSpaces(inPolytopeToClip[i], clipper_planes, tolerance);
	if (all_inside)
	{
		for (int i = 0; i < n_clip && outClippedPolytope.size() < outClippedPolytope.capacity(); ++i)
		{
			bool duplicate = false;
			for (const Vec4 &existing : outClippedPolytope)
				if ((existing - inPolytopeToClip[i]).LengthSq() < Square(tolerance))
				{
					duplicate = true;
					break;
				}
			if (!duplicate)
				outClippedPolytope.push_back(inPolytopeToClip[i]);
		}
		return;
	}

	ClipHalfSpaces to_clip_planes;
	GetClipHalfSpaces(inPolytopeToClip, inAxis, to_clip_planes);

	// The 3-flat of the polytope being clipped, so results can be lifted back onto it along the axis
	Vec4 plane_origin = inPolytopeToClip[0];
	Vec4 plane_normal = Vec4::sZero();
	for (int j = 1; j < n_clip && plane_normal.IsNearZero(); ++j)
		for (int k = j + 1; k < n_clip && plane_normal.IsNearZero(); ++k)
			for (int l = k + 1; l < n_clip && plane_normal.IsNearZero(); ++l)
			{
				Vec4 candidate = Vec4::sCross(inPolytopeToClip[j] - plane_origin, inPolytopeToClip[k] - plane_origin, inPolytopeToClip[l] - plane_origin);
				if (candidate.LengthSq() > FLT_EPSILON)
					plane_normal = candidate;
			}
	float axis_dot_normal = inAxis.Dot(plane_normal);
	bool can_project = !plane_normal.IsNearZero() && abs(axis_dot_normal) > FLT_EPSILON;

	// Add a point, lifting it onto the 3-flat of the clipped polytope and skipping near duplicates
	auto add_point = [&](Vec4Arg inPoint)
	{
		if (outClippedPolytope.size() >= outClippedPolytope.capacity())
			return;

		Vec4 p = inPoint;
		if (can_project)
			p = p - ((p - plane_origin).Dot(plane_normal) / axis_dot_normal) * inAxis;

		for (const Vec4 &existing : outClippedPolytope)
			if ((existing - p).LengthSq() < Square(tolerance))
				return;

		outClippedPolytope.push_back(p);
	};

	// Vertices of one polytope that are inside the other, plus the crossings of all vertex pairs.
	// Each vertex's signed distance to each halfspace is computed once up front: the planes are
	// affine, so a segment's entry/exit parameters follow from its endpoints' distances alone
	// (n . (B - A) is just dist_B - dist_A), which takes every dot product out of the O(n^2)
	// segment loop below. For a cube cell against a cube cell that is 28 segments x 6 planes of
	// arithmetic instead of 28 x 6 x 2 dot products.
	auto clip_against = [&](const VERTEX_ARRAY &inVerts, int inNumVerts, const ClipHalfSpaces &inPlanes, const ClipHalfSpaces &inOwnPlanes)
	{
		if (inPlanes.empty())
			return;

		int n_planes = int(inPlanes.size());
		float dist[cMaxClipVertices][cMaxClipHalfSpaces];
		int n = min(inNumVerts, int(cMaxClipVertices));
		for (int i = 0; i < n; ++i)
			for (int h = 0; h < n_planes; ++h)
				dist[i][h] = inPlanes[h].mNormal.Dot(inVerts[i]) - inPlanes[h].mDistance;

		for (int i = 0; i < n; ++i)
		{
			bool inside = true;
			for (int h = 0; h < n_planes && inside; ++h)
				inside = dist[i][h] <= tolerance;
			if (inside)
				add_point(inVerts[i]);
		}

		// Which of this polytope's own facets each vertex lies on. Two vertices of a convex polytope
		// span an edge exactly when they share at least two facets, so this picks the real edges out
		// of the n^2 vertex pairs - 12 instead of 28 for a cube cell. It is not only cheaper: a
		// chord that is not an edge can cross the other polytope's boundary at points that are not
		// vertices of the intersection, which is where the spurious extra points came from.
		uint64 facet_mask[cMaxClipVertices];
		int n_own = min(int(inOwnPlanes.size()), 64);
		bool have_edges = n_own >= 2;
		if (have_edges)
			for (int i = 0; i < n; ++i)
			{
				uint64 m = 0;
				for (int h = 0; h < n_own; ++h)
					if (abs(inOwnPlanes[h].mNormal.Dot(inVerts[i]) - inOwnPlanes[h].mDistance) <= tolerance)
						m |= uint64(1) << h;
				facet_mask[i] = m;
			}

		for (int i = 0; i < n; ++i)
			for (int j = i + 1; j < n; ++j)
			{
				if (have_edges && CountBits(uint32(facet_mask[i] & facet_mask[j])) + CountBits(uint32((facet_mask[i] & facet_mask[j]) >> 32)) < 2)
					continue;

				float t0 = 0.0f, t1 = 1.0f;
				bool keep = true;
				for (int h = 0; h < n_planes; ++h)
				{
					float da = dist[i][h];
					float dd = dist[j][h] - da;
					if (abs(dd) < FLT_EPSILON)
					{
						if (da > tolerance)
						{
							keep = false;
							break;
						}
						continue;
					}
					float t = -da / dd;
					if (dd > 0.0f)
						t1 = min(t1, t);	// Leaving the halfspace
					else
						t0 = max(t0, t);	// Entering it
					if (t0 > t1)
					{
						keep = false;
						break;
					}
				}
				if (keep)
				{
					Vec4 d = inVerts[j] - inVerts[i];
					if (t0 > 0.0f) add_point(inVerts[i] + t0 * d);
					if (t1 < 1.0f) add_point(inVerts[i] + t1 * d);
				}
			}
	};

	clip_against(inPolytopeToClip, n_clip, clipper_planes, to_clip_planes);
	clip_against(inClippingPolytope, n_clipper, to_clip_planes, clipper_planes);
}

#ifdef JPH_ENABLE_2D_POLY_CLIPPING
// TODO(4D): never defined. See the note at the top of this file - these walk an ordered 2D
// vertex loop, which a 4D supporting face is not. Kept verbatim for reference.

/// Clip inPolygonToClip against the positive halfspace of plane defined by inPlaneOrigin and inPlaneNormal.
/// inPlaneNormal does not need to be normalized.
template <class VERTEX_ARRAY>
void ClipPolyVsPlane(const VERTEX_ARRAY &inPolygonToClip, Vec4Arg inPlaneOrigin, Vec4Arg inPlaneNormal, VERTEX_ARRAY &outClippedPolygon)
{
	JPH_ASSERT(inPolygonToClip.size() >= 2);
	JPH_ASSERT(outClippedPolygon.empty());

	// Determine state of last point
	Vec4 e1 = inPolygonToClip[inPolygonToClip.size() - 1];
	float prev_num = (inPlaneOrigin - e1).Dot(inPlaneNormal);
	bool prev_inside = prev_num < 0.0f;

	// Loop through all vertices
	for (typename VERTEX_ARRAY::size_type j = 0; j < inPolygonToClip.size(); ++j)
	{
		// Check if second point is inside
		Vec4Arg e2 = inPolygonToClip[j];
		float num = (inPlaneOrigin - e2).Dot(inPlaneNormal);
		bool cur_inside = num < 0.0f;

		// In -> Out or Out -> In: Add point on clipping plane
		if (cur_inside != prev_inside)
		{
			// Solve: (X - inPlaneOrigin) . inPlaneNormal = 0 and X = e1 + t * (e2 - e1) for X
			Vec4 e12 = e2 - e1;
			float denom = e12.Dot(inPlaneNormal);
			if (denom != 0.0f)
				outClippedPolygon.push_back(e1 + (prev_num / denom) * e12);
			else
				cur_inside = prev_inside; // Edge is parallel to plane, treat point as if it were on the same side as the last point
		}

		// Point inside, add it
		if (cur_inside)
			outClippedPolygon.push_back(e2);

		// Update previous state
		prev_num = num;
		prev_inside = cur_inside;
		e1 = e2;
	}
}

/// Clip polygon versus polygon.
/// Both polygons are assumed to be in counter clockwise order.
/// @param inClippingPolygonNormal is the hyperplane normal of the clipping polygon's plane, does not need to be normalized
/// @param inClippingPolygon is the polygon which inClippedPolygon is clipped against
/// @param inPolygonToClip is the polygon that is clipped
/// @param outClippedPolygon will contain clipped polygon when function returns
template <class VERTEX_ARRAY>
void ClipPolyVsPoly(const VERTEX_ARRAY &inPolygonToClip, const VERTEX_ARRAY &inClippingPolygon, Vec4Arg inClippingPolygonNormal, VERTEX_ARRAY &outClippedPolygon)
{
	JPH_ASSERT(inPolygonToClip.size() >= 2);
	JPH_ASSERT(inClippingPolygon.size() >= 3);

	// Compute the second normal of the clipping polygon's 2D plane in 4D.
	// The polygon lies in a 2D subspace. inClippingPolygonNormal is one normal direction.
	// We derive the second from two polygon edges and the first normal using the 4D cross product.
	Vec4 poly_e1 = inClippingPolygon[1] - inClippingPolygon[0];
	Vec4 poly_e2 = inClippingPolygon[2] - inClippingPolygon[0];
	Vec4 second_normal = Vec4::sCross(poly_e1, poly_e2, inClippingPolygonNormal);

	VERTEX_ARRAY tmp_vertices[2];
	int tmp_vertices_idx = 0;

	for (typename VERTEX_ARRAY::size_type i = 0; i < inClippingPolygon.size(); ++i)
	{
		// Get edge to clip against
		Vec4 clip_e1 = inClippingPolygon[i];
		Vec4 clip_e2 = inClippingPolygon[(i + 1) % inClippingPolygon.size()];
		Vec4 clip_edge = clip_e2 - clip_e1;

		// Compute inward-pointing normal perpendicular to the edge, lying in the polygon plane.
		// sCross(edge, N1, N2) gives a vector perpendicular to all three; since N1 and N2 span
		// the normal space, the result lies in the polygon plane and is perpendicular to the edge.
		Vec4 clip_normal = Vec4::sCross(clip_edge, inClippingPolygonNormal, second_normal);

		// Get source and target polygon
		const VERTEX_ARRAY &src_polygon = (i == 0)? inPolygonToClip : tmp_vertices[tmp_vertices_idx];
		tmp_vertices_idx ^= 1;
		VERTEX_ARRAY &tgt_polygon = (i == inClippingPolygon.size() - 1)? outClippedPolygon : tmp_vertices[tmp_vertices_idx];
		tgt_polygon.clear();

		// Clip against the edge
		ClipPolyVsPlane(src_polygon, clip_e1, clip_normal, tgt_polygon);

		// Break out if no polygon left
		if (tgt_polygon.size() < 3)
		{
			outClippedPolygon.clear();
			break;
		}
	}
}

/// Clip inPolygonToClip against an edge, the edge is projected on inPolygonToClip using inClippingEdgeNormal.
/// The positive half space (the side on the edge in the direction of inClippingEdgeNormal) is cut away.
/// @param inClippingEdgeNormal is the hyperplane normal of the polygon's plane
template <class VERTEX_ARRAY>
void ClipPolyVsEdge(const VERTEX_ARRAY &inPolygonToClip, Vec4Arg inEdgeVertex1, Vec4Arg inEdgeVertex2, Vec4Arg inClippingEdgeNormal, VERTEX_ARRAY &outClippedPolygon)
{
	JPH_ASSERT(inPolygonToClip.size() >= 3);
	JPH_ASSERT(outClippedPolygon.empty());

	// Compute the second normal of the polygon's 2D plane
	Vec4 pe1 = inPolygonToClip[1] - inPolygonToClip[0];
	Vec4 pe2 = inPolygonToClip[2] - inPolygonToClip[0];
	Vec4 second_normal = Vec4::sCross(pe1, pe2, inClippingEdgeNormal);

	// Get normal perpendicular to the edge, lying in the polygon plane
	Vec4 edge = inEdgeVertex2 - inEdgeVertex1;
	Vec4 edge_normal = Vec4::sCross(edge, inClippingEdgeNormal, second_normal);

	// Project vertices of edge onto the polygon's 2D affine plane.
	// The plane passes through inPolygonToClip[0] with normal space spanned by {N1, N2}.
	// Projection: v_proj = v - N1*(N1.(v-p0)/N1.N1) - N2*(N2.(v-p0)/N2.N2)
	// (assuming N1, N2 are not necessarily orthogonal — use Gram-Schmidt first)
	Vec4 p0 = inPolygonToClip[0];
	Vec4 n1 = inClippingEdgeNormal;
	Vec4 n2 = second_normal;
	// Gram-Schmidt: make n2 orthogonal to n1
	float n1_len_sq = n1.LengthSq();
	if (n1_len_sq > FLT_EPSILON)
		n2 = n2 - n1 * (n2.Dot(n1) / n1_len_sq);
	float n2_len_sq = n2.LengthSq();

	Vec4 d1 = inEdgeVertex1 - p0;
	Vec4 v1 = inEdgeVertex1;
	if (n1_len_sq > FLT_EPSILON)
		v1 = v1 - n1 * (n1.Dot(d1) / n1_len_sq);
	if (n2_len_sq > FLT_EPSILON)
		v1 = v1 - n2 * (n2.Dot(d1) / n2_len_sq);

	Vec4 d2 = inEdgeVertex2 - p0;
	Vec4 v2 = inEdgeVertex2;
	if (n1_len_sq > FLT_EPSILON)
		v2 = v2 - n1 * (n1.Dot(d2) / n1_len_sq);
	if (n2_len_sq > FLT_EPSILON)
		v2 = v2 - n2 * (n2.Dot(d2) / n2_len_sq);

	Vec4 v12 = v2 - v1;
	float v12_len_sq = v12.LengthSq();

	// Determine state of last point
	Vec4 e1 = inPolygonToClip[inPolygonToClip.size() - 1];
	float prev_num = (inEdgeVertex1 - e1).Dot(edge_normal);
	bool prev_inside = prev_num < 0.0f;

	// Loop through all vertices
	for (typename VERTEX_ARRAY::size_type j = 0; j < inPolygonToClip.size(); ++j)
	{
		// Check if second point is inside
		Vec4 e2 = inPolygonToClip[j];
		float num = (inEdgeVertex1 - e2).Dot(edge_normal);
		bool cur_inside = num < 0.0f;

		// In -> Out or Out -> In: Add point on clipping plane
		if (cur_inside != prev_inside)
		{
			// Solve: (inEdgeVertex1 - X) . edge_normal = 0 and X = e1 + t * (e2 - e1) for X
			Vec4 e12 = e2 - e1;
			float denom = e12.Dot(edge_normal);
			Vec4 clipped_point = denom != 0.0f? e1 + (prev_num / denom) * e12 : e1;

			// Project point on line segment v1, v2 so see if it falls outside if the edge
			float projection = (clipped_point - v1).Dot(v12);
			if (projection < 0.0f)
				outClippedPolygon.push_back(v1);
			else if (projection > v12_len_sq)
				outClippedPolygon.push_back(v2);
			else
				outClippedPolygon.push_back(clipped_point);
		}

		// Update previous state
		prev_num = num;
		prev_inside = cur_inside;
		e1 = e2;
	}
}

/// Clip polygon vs axis aligned box, inPolygonToClip is assume to be in counter clockwise order.
/// Output will be stored in outClippedPolygon. Everything inside inAABox will be kept.
template <class VERTEX_ARRAY>
void ClipPolyVsAABox(const VERTEX_ARRAY &inPolygonToClip, const AABox &inAABox, VERTEX_ARRAY &outClippedPolygon)
{
	JPH_ASSERT(inPolygonToClip.size() >= 2);

	VERTEX_ARRAY tmp_vertices[2];
	int tmp_vertices_idx = 0;

	for (int coord = 0; coord < 4; ++coord)
		for (int side = 0; side < 2; ++side)
		{
			// Get plane to clip against
			Vec4 origin = Vec4::sZero(), normal = Vec4::sZero();
			if (side == 0)
			{
				normal.SetComponent(coord, 1.0f);
				origin.SetComponent(coord, inAABox.mMin[coord]);
			}
			else
			{
				normal.SetComponent(coord, -1.0f);
				origin.SetComponent(coord, inAABox.mMax[coord]);
			}

			// Get source and target polygon
			const VERTEX_ARRAY &src_polygon = tmp_vertices_idx == 0? inPolygonToClip : tmp_vertices[tmp_vertices_idx & 1];
			tmp_vertices_idx++;
			VERTEX_ARRAY &tgt_polygon = tmp_vertices_idx == 8? outClippedPolygon : tmp_vertices[tmp_vertices_idx & 1];
			tgt_polygon.clear();

			// Clip against the edge
			ClipPolyVsPlane(src_polygon, origin, normal, tgt_polygon);

			// Break out if no polygon left
			if (tgt_polygon.size() < 3)
			{
				outClippedPolygon.clear();
				return;
			}

			// Flip normal
			normal = -normal;
		}
}

#endif // JPH_ENABLE_2D_POLY_CLIPPING

JPH_NAMESPACE_END
