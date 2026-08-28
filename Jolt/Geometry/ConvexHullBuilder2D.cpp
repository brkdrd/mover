// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Geometry/ConvexHullBuilder2D.h>

#ifdef JPH_CONVEX_BUILDER_2D_DEBUG
	#include <Jolt/Renderer/DebugRenderer.h>
#endif

JPH_NAMESPACE_BEGIN

void ConvexHullBuilder2D::Edge::CalculateNormalAndCenter(const Vec4 *inPositions)
{
	Vec4 p1 = inPositions[mStartIdx];
	Vec4 p2 = inPositions[mNextEdge->mStartIdx];

	// Center of edge
	mCenter = 0.5f * (p1 + p2);

	// Create outward pointing normal.
	// We have two choices for the normal (which satisfies normal . edge = 0):
	// normal1 = (-edge.y, edge.x, 0, 0)
	// normal2 = (edge.y, -edge.x, 0, 0)
	// We want the normal to point out of the polygon (CCW winding). Only normal2 satisfies this condition.
	Vec4 edge = p2 - p1;
	mNormal = Vec4(edge.GetY(), -edge.GetX(), 0, 0);
}

ConvexHullBuilder2D::ConvexHullBuilder2D(const Positions &inPositions) :
	mPositions(inPositions)
{
#ifdef JPH_CONVEX_BUILDER_2D_DEBUG
	// Center the drawing of the first hull around the origin and calculate the delta offset between states
	mOffset = RVec3::sZero();
	if (mPositions.empty())
	{
		// No hull will be generated
		mDelta = Vec3::sZero();
	}
	else
	{
		Vec4 maxv = Vec4::sReplicate(-FLT_MAX), minv = Vec4::sReplicate(FLT_MAX);
		for (Vec4 v : mPositions)
		{
			minv = Vec4::sMin(minv, v);
			maxv = Vec4::sMax(maxv, v);
			mOffset -= Vec3(v.GetX(), v.GetY(), v.GetZ());
		}
		mOffset /= Real(mPositions.size());
		mDelta = Vec3((maxv - minv).GetX() + 0.5f, 0, 0);
		mOffset += mDelta; // Don't start at origin, we're already drawing the final hull there
	}
#endif
}

ConvexHullBuilder2D::~ConvexHullBuilder2D()
{
	FreeEdges();
}

void ConvexHullBuilder2D::FreeEdges()
{
	if (mFirstEdge == nullptr)
		return;

	Edge *edge = mFirstEdge;
	do
	{
		Edge *next = edge->mNextEdge;
		delete edge;
		edge = next;
	} while (edge != mFirstEdge);

	mFirstEdge = nullptr;
	mNumEdges = 0;
}

#ifdef JPH_ENABLE_ASSERTS

void ConvexHullBuilder2D::ValidateEdges() const
{
	if (mFirstEdge == nullptr)
	{
		JPH_ASSERT(mNumEdges == 0);
		return;
	}

	int count = 0;

	Edge *edge = mFirstEdge;
	do
	{
		// Validate connectivity
		JPH_ASSERT(edge->mNextEdge->mPrevEdge == edge);
		JPH_ASSERT(edge->mPrevEdge->mNextEdge == edge);

		++count;
		edge = edge->mNextEdge;
	} while (edge != mFirstEdge);

	// Validate that count matches
	JPH_ASSERT(count == mNumEdges);
}

#endif // JPH_ENABLE_ASSERTS

void ConvexHullBuilder2D::AssignPointToEdge(int inPositionIdx, const Array<Edge *> &inEdges) const
{
	Vec4 point = mPositions[inPositionIdx];

	Edge *best_edge = nullptr;
	float best_dist_sq = 0.0f;

	// Test against all edges
	for (Edge *edge : inEdges)
	{
		// Determine distance to edge
		float dot = edge->mNormal.Dot(point - edge->mCenter);
		if (dot > 0.0f)
		{
			float dist_sq = dot * dot / edge->mNormal.LengthSq();
			if (dist_sq > best_dist_sq)
			{
				best_edge = edge;
				best_dist_sq = dist_sq;
			}
		}
	}

	// If this point is in front of the edge, add it to the conflict list
	if (best_edge != nullptr)
	{
		if (best_dist_sq > best_edge->mFurthestPointDistanceSq)
		{
			// This point is further away than any others, update the distance and add point as last point
			best_edge->mFurthestPointDistanceSq = best_dist_sq;
			best_edge->mConflictList.push_back(inPositionIdx);
		}
		else
		{
			// Not the furthest point, add it as the before last point
			best_edge->mConflictList.insert(best_edge->mConflictList.begin() + best_edge->mConflictList.size() - 1, inPositionIdx);
		}
	}
}

ConvexHullBuilder2D::EResult ConvexHullBuilder2D::Initialize(int inIdx1, int inIdx2, int inIdx3, int inMaxVertices, float inTolerance, Edges &outEdges)
{
	// Clear any leftovers
	FreeEdges();
	outEdges.clear();

	// Reset flag
	EResult result = EResult::Success;

	// Determine a suitable tolerance for detecting that points are colinear
	// Formula as per: Implementing Quickhull - Dirk Gregorius.
	Vec4 vmax = Vec4::sZero();
	for (Vec4 v : mPositions)
		vmax = Vec4::sMax(vmax, v.Abs());
	float colinear_tolerance_sq = Square(2.0f * FLT_EPSILON * (vmax.GetX() + vmax.GetY()));

	// Increase desired tolerance if accuracy doesn't allow it
	float tolerance_sq = max(colinear_tolerance_sq, Square(inTolerance));

	// Start with the initial indices in counter clockwise order
	// 2D cross product (determinant of XY components) to determine winding
	Vec4 edge12 = mPositions[inIdx2] - mPositions[inIdx1];
	Vec4 edge13 = mPositions[inIdx3] - mPositions[inIdx1];
	float z = edge12.GetX() * edge13.GetY() - edge12.GetY() * edge13.GetX();
	if (z < 0.0f)
		std::swap(inIdx1, inIdx2);

	// Create and link edges
	Edge *e1 = new Edge(inIdx1);
	Edge *e2 = new Edge(inIdx2);
	Edge *e3 = new Edge(inIdx3);
	e1->mNextEdge = e2;
	e1->mPrevEdge = e3;
	e2->mNextEdge = e3;
	e2->mPrevEdge = e1;
	e3->mNextEdge = e1;
	e3->mPrevEdge = e2;
	mFirstEdge = e1;
	mNumEdges = 3;

	// Build the initial conflict lists
	Array<Edge *> edges { e1, e2, e3 };
	for (Edge *edge : edges)
		edge->CalculateNormalAndCenter(mPositions.data());
	for (int idx = 0; idx < (int)mPositions.size(); ++idx)
		if (idx != inIdx1 && idx != inIdx2 && idx != inIdx3)
			AssignPointToEdge(idx, edges);

	JPH_IF_ENABLE_ASSERTS(ValidateEdges();)
#ifdef JPH_CONVEX_BUILDER_2D_DEBUG
	DrawState();
#endif

	// Add the remaining points to the hull
	for (;;)
	{
		// Check if we've reached the max amount of vertices that are allowed
		if (mNumEdges >= inMaxVertices)
		{
			result = EResult::MaxVerticesReached;
			break;
		}

		// Find the edge with the furthest point on it
		Edge *edge_with_furthest_point = nullptr;
		float furthest_dist_sq = 0.0f;
		Edge *edge = mFirstEdge;
		do
		{
			if (edge->mFurthestPointDistanceSq > furthest_dist_sq)
			{
				furthest_dist_sq = edge->mFurthestPointDistanceSq;
				edge_with_furthest_point = edge;
			}
			edge = edge->mNextEdge;
		} while (edge != mFirstEdge);

		// If there is none closer than our tolerance value, we're done
		if (edge_with_furthest_point == nullptr || furthest_dist_sq < tolerance_sq)
			break;

		// Take the furthest point
		int furthest_point_idx = edge_with_furthest_point->mConflictList.back();
		edge_with_furthest_point->mConflictList.pop_back();
		Vec4 furthest_point = mPositions[furthest_point_idx];

		// Find the horizon of edges that need to be removed
		Edge *first_edge = edge_with_furthest_point;
		do
		{
			Edge *prev = first_edge->mPrevEdge;
			if (!prev->IsFacing(furthest_point))
				break;
			first_edge = prev;
		} while (first_edge != edge_with_furthest_point);

		Edge *last_edge = edge_with_furthest_point;
		do
		{
			Edge *next = last_edge->mNextEdge;
			if (!next->IsFacing(furthest_point))
				break;
			last_edge = next;
		} while (last_edge != edge_with_furthest_point);

		// Create new edges
		e1 = new Edge(first_edge->mStartIdx);
		e2 = new Edge(furthest_point_idx);
		e1->mNextEdge = e2;
		e1->mPrevEdge = first_edge->mPrevEdge;
		e2->mPrevEdge = e1;
		e2->mNextEdge = last_edge->mNextEdge;
		e1->mPrevEdge->mNextEdge = e1;
		e2->mNextEdge->mPrevEdge = e2;
		mFirstEdge = e1; // We could delete mFirstEdge so just update it to the newly created edge
		mNumEdges += 2;

		// Calculate normals
		Array<Edge *> new_edges { e1, e2 };
		for (Edge *new_edge : new_edges)
			new_edge->CalculateNormalAndCenter(mPositions.data());

		// Delete the old edges
		for (;;)
		{
			Edge *next = first_edge->mNextEdge;

			// Redistribute points in conflict list
			for (int idx : first_edge->mConflictList)
				AssignPointToEdge(idx, new_edges);

			// Delete the old edge
			delete first_edge;
			--mNumEdges;

			if (first_edge == last_edge)
				break;
			first_edge = next;
		}

		JPH_IF_ENABLE_ASSERTS(ValidateEdges();)
	#ifdef JPH_CONVEX_BUILDER_2D_DEBUG
		DrawState();
	#endif
	}

	// Convert the edge list to a list of indices
	outEdges.reserve(mNumEdges);
	Edge *edge = mFirstEdge;
	do
	{
		outEdges.push_back(edge->mStartIdx);
		edge = edge->mNextEdge;
	} while (edge != mFirstEdge);

	return result;
}

#ifdef JPH_CONVEX_BUILDER_2D_DEBUG

void ConvexHullBuilder2D::DrawState()
{
	int color_idx = 0;

	const Edge *edge = mFirstEdge;
	do
	{
		const Edge *next = edge->mNextEdge;

		// Get unique color per edge
		Color color = Color::sGetDistinctColor(color_idx++);

		// Project Vec4 positions to Vec3 for debug drawing
		auto ToVec3 = [](Vec4 v) -> Vec3 { return Vec3(v.GetX(), v.GetY(), v.GetZ()); };

		// Draw edge and normal
		DebugRenderer::sInstance->DrawArrow(cDrawScale * (mOffset + ToVec3(mPositions[edge->mStartIdx])), cDrawScale * (mOffset + ToVec3(mPositions[next->mStartIdx])), color, 0.1f);
		Vec3 center3 = ToVec3(edge->mCenter);
		Vec3 normal3 = ToVec3(edge->mNormal);
		float normal_len = normal3.Length();
		Vec3 normalized_normal = normal_len > 0.0f ? normal3 / normal_len : Vec3::sZero();
		DebugRenderer::sInstance->DrawArrow(cDrawScale * (mOffset + center3), cDrawScale * (mOffset + center3) + normalized_normal, Color::sGreen, 0.1f);

		// Draw points that belong to this edge in the same color
		for (int idx : edge->mConflictList)
			DebugRenderer::sInstance->DrawMarker(cDrawScale * (mOffset + ToVec3(mPositions[idx])), color, 0.05f);

		edge = next;
	} while (edge != mFirstEdge);

	mOffset += mDelta;
}

#endif

JPH_NAMESPACE_END
