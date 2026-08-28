// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Geometry/ConvexHullBuilder.h>
#include <Jolt/Geometry/ConvexHullBuilder2D.h>
#include <Jolt/Geometry/ClosestPoint.h>
#include <Jolt/Core/StringTools.h>
#include <Jolt/Core/UnorderedSet.h>

#ifdef JPH_CONVEX_BUILDER_DUMP_SHAPE
JPH_SUPPRESS_WARNINGS_STD_BEGIN
#include <fstream>
JPH_SUPPRESS_WARNINGS_STD_END
#endif // JPH_CONVEX_BUILDER_DUMP_SHAPE

#ifdef JPH_CONVEX_BUILDER_DEBUG
	#include <Jolt/Renderer/DebugRenderer.h>
#endif

JPH_NAMESPACE_BEGIN

void ConvexHullBuilder::Cell::CalculateNormalAndCentroid(const Vec4 *inPositions)
{
	Vec4 v0 = inPositions[mIdx[0]];
	Vec4 v1 = inPositions[mIdx[1]];
	Vec4 v2 = inPositions[mIdx[2]];
	Vec4 v3 = inPositions[mIdx[3]];

	// 4D cross product of edge vectors gives the outward normal
	mNormal = Vec4::sCross(v1 - v0, v2 - v0, v3 - v0);

	// Centroid is the average of the 4 vertices
	mCentroid = 0.25f * (v0 + v1 + v2 + v3);
}

void ConvexHullBuilder::Cell::Initialize(int inIdx0, int inIdx1, int inIdx2, int inIdx3, const Vec4 *inPositions)
{
	JPH_ASSERT(inIdx0 != inIdx1 && inIdx0 != inIdx2 && inIdx0 != inIdx3);
	JPH_ASSERT(inIdx1 != inIdx2 && inIdx1 != inIdx3 && inIdx2 != inIdx3);

	mIdx[0] = inIdx0;
	mIdx[1] = inIdx1;
	mIdx[2] = inIdx2;
	mIdx[3] = inIdx3;

	CalculateNormalAndCentroid(inPositions);
}

ConvexHullBuilder::ConvexHullBuilder(const Positions &inPositions) :
	mPositions(inPositions)
{
#ifdef JPH_CONVEX_BUILDER_DEBUG
	mIteration = 0;
#endif
}

void ConvexHullBuilder::FreeCells()
{
	for (Cell *c : mCells)
		delete c;
	mCells.clear();
}

void ConvexHullBuilder::GetCellForPoint(Vec4Arg inPoint, const Cells &inCells, Cell *&outCell, float &outDistSq) const
{
	outCell = nullptr;
	outDistSq = 0.0f;

	for (Cell *c : inCells)
		if (!c->mRemoved)
		{
			float dot = c->mNormal.Dot(inPoint - c->mCentroid);
			if (dot > 0.0f)
			{
				float dist_sq = dot * dot / c->mNormal.LengthSq();
				if (dist_sq > outDistSq)
				{
					outCell = c;
					outDistSq = dist_sq;
				}
			}
		}
}

float ConvexHullBuilder::GetDistanceToRidgeSq(Vec4Arg inPoint, const Cell *inCell) const
{
	// Use ClosestPoint::GetClosestPointOnTetrahedron to determine if the point projects inside
	Vec4 v0 = mPositions[inCell->mIdx[0]] - inPoint;
	Vec4 v1 = mPositions[inCell->mIdx[1]] - inPoint;
	Vec4 v2 = mPositions[inCell->mIdx[2]] - inPoint;
	Vec4 v3 = mPositions[inCell->mIdx[3]] - inPoint;

	uint32 set = 0;
	Vec4 closest = ClosestPoint::GetClosestPointOnTetrahedron(v0, v1, v2, v3, set);

	// If point projects into the interior of the tetrahedron, distance to boundary is 0
	if (set == 0b1111)
		return 0.0f;

	return closest.LengthSq();
}

bool ConvexHullBuilder::AssignPointToCell(int inPositionIdx, const Cells &inCells, float inToleranceSq)
{
	Vec4 point = mPositions[inPositionIdx];

	// Find the cell for which the point is furthest away
	Cell *best_cell;
	float best_dist_sq;
	GetCellForPoint(point, inCells, best_cell, best_dist_sq);

	if (best_cell != nullptr)
	{
		// Check if this point is within tolerance of the hyperplane
		if (best_dist_sq <= inToleranceSq)
		{
			// Check distance to cell boundary
			float dist_to_ridge_sq = GetDistanceToRidgeSq(point, best_cell);
			if (dist_to_ridge_sq > inToleranceSq)
			{
				// Point is outside the cell projection and too far to discard
				mCoplanarList.push_back({ inPositionIdx, dist_to_ridge_sq });
			}
		}
		else
		{
			// Point is in front of the cell, add to conflict list
			if (best_dist_sq > best_cell->mFurthestPointDistanceSq)
			{
				best_cell->mFurthestPointDistanceSq = best_dist_sq;
				best_cell->mConflictList.push_back(inPositionIdx);
			}
			else
			{
				best_cell->mConflictList.insert(best_cell->mConflictList.begin() + best_cell->mConflictList.size() - 1, inPositionIdx);
			}

			return true;
		}
	}

	return false;
}

float ConvexHullBuilder::DetermineCoplanarDistance() const
{
	// Formula adapted for 4D: 4 * eps * (|max_x| + |max_y| + |max_z| + |max_w|)
	Vec4 vmax = Vec4::sZero();
	for (const Vec4 &v : mPositions)
		vmax = Vec4::sMax(vmax, v.Abs());
	return 4.0f * FLT_EPSILON * (vmax.GetX() + vmax.GetY() + vmax.GetZ() + vmax.GetW());
}

int ConvexHullBuilder::GetNumVerticesUsed() const
{
	UnorderedSet<int> used_verts;
	used_verts.reserve(UnorderedSet<int>::size_type(mPositions.size()));
	for (const Cell *c : mCells)
		if (!c->mRemoved)
			for (int i = 0; i < 4; ++i)
				used_verts.insert(c->mIdx[i]);
	return (int)used_verts.size();
}

bool ConvexHullBuilder::ContainsCell(const Array<int> &inIndices) const
{
	for (const Cell *c : mCells)
	{
		if (c->mRemoved)
			continue;

		// Check if the cell's vertex set matches (order-independent)
		int matches = 0;
		for (int idx : inIndices)
			for (int j = 0; j < 4; ++j)
				if (c->mIdx[j] == idx)
				{
					matches++;
					break;
				}

		if (matches == (int)inIndices.size() && (int)inIndices.size() == 4)
			return true;
	}

	return false;
}

ConvexHullBuilder::EResult ConvexHullBuilder::Initialize(int inMaxVertices, float inTolerance, const char *&outError)
{
	// Free cells from earlier hull
	FreeCells();

	// Need at least 5 points for a 4D hull (pentachoron)
	if (mPositions.size() < 5)
	{
		outError = "Need at least 5 points to make a 4D hull";
		return EResult::TooFewPoints;
	}

	// Determine coplanar tolerance
	float coplanar_tolerance_sq = Square(DetermineCoplanarDistance());
	float tolerance_sq = max(coplanar_tolerance_sq, Square(inTolerance));

	// === Step 1: Find the two most distant points ===
	// Find point furthest from origin
	int idx1 = -1;
	float max_dist_sq = -1.0f;
	for (int i = 0; i < (int)mPositions.size(); ++i)
	{
		float dist_sq = mPositions[i].LengthSq();
		if (dist_sq > max_dist_sq)
		{
			max_dist_sq = dist_sq;
			idx1 = i;
		}
	}
	JPH_ASSERT(idx1 >= 0);

	// Find point furthest from idx1
	int idx2 = -1;
	max_dist_sq = -1.0f;
	for (int i = 0; i < (int)mPositions.size(); ++i)
		if (i != idx1)
		{
			float dist_sq = (mPositions[i] - mPositions[idx1]).LengthSq();
			if (dist_sq > max_dist_sq)
			{
				max_dist_sq = dist_sq;
				idx2 = i;
			}
		}
	JPH_ASSERT(idx2 >= 0);

	// === Step 2: Find 3rd point maximizing triangle area (via Gram determinant) ===
	Vec4 e_base = mPositions[idx2] - mPositions[idx1];
	float e_base_len_sq = e_base.LengthSq();
	int idx3 = -1;
	float best_gram_sq = -1.0f;
	for (int i = 0; i < (int)mPositions.size(); ++i)
		if (i != idx1 && i != idx2)
		{
			Vec4 ei = mPositions[i] - mPositions[idx1];
			// Gram determinant = |e_base|^2 * |ei|^2 - (e_base . ei)^2  (= 4 * triangle_area^2)
			float gram = e_base_len_sq * ei.LengthSq() - Square(e_base.Dot(ei));
			if (gram > best_gram_sq)
			{
				best_gram_sq = gram;
				idx3 = i;
			}
		}
	JPH_ASSERT(idx3 >= 0);
	if (best_gram_sq < cMinTriangleGramSq)
	{
		outError = "Could not find a non-degenerate initial triangle";
		return EResult::Degenerate;
	}

	// === Step 3: Find 4th point maximizing distance from triangle's 2D span ===
	// Use Gram-Schmidt to orthogonalize the triangle's spanning vectors
	Vec4 e1 = mPositions[idx2] - mPositions[idx1];
	Vec4 e2_raw = mPositions[idx3] - mPositions[idx1];
	Vec4 e1_norm = e1 / e1.Length();
	Vec4 e2_orth = e2_raw - e1_norm * e2_raw.Dot(e1_norm);
	float e2_orth_len = e2_orth.Length();
	if (e2_orth_len < FLT_EPSILON)
	{
		outError = "Initial triangle is degenerate";
		return EResult::Degenerate;
	}
	Vec4 e2_norm = e2_orth / e2_orth_len;

	int idx4 = -1;
	float max_dist_from_plane = 0.0f;
	for (int i = 0; i < (int)mPositions.size(); ++i)
		if (i != idx1 && i != idx2 && i != idx3)
		{
			Vec4 v = mPositions[i] - mPositions[idx1];
			// Project onto triangle span, measure residual
			float proj1 = v.Dot(e1_norm);
			float proj2 = v.Dot(e2_norm);
			float dist_sq_from_span = v.LengthSq() - proj1 * proj1 - proj2 * proj2;
			if (dist_sq_from_span < 0.0f) dist_sq_from_span = 0.0f; // numerical safety
			if (dist_sq_from_span > max_dist_from_plane)
			{
				max_dist_from_plane = dist_sq_from_span;
				idx4 = i;
			}
		}

	if (idx4 < 0 || max_dist_from_plane <= Square(cCoplanarSlopFactor) * coplanar_tolerance_sq)
	{
		outError = "All points are co-planar in 4D";
		return EResult::Degenerate;
	}

	// === Step 4: Find 5th point maximizing distance from tetrahedron's hyperplane ===
	Vec4 tet_e1 = mPositions[idx2] - mPositions[idx1];
	Vec4 tet_e2 = mPositions[idx3] - mPositions[idx1];
	Vec4 tet_e3 = mPositions[idx4] - mPositions[idx1];
	Vec4 hyperplane_normal = Vec4::sCross(tet_e1, tet_e2, tet_e3);
	float hyperplane_normal_len = hyperplane_normal.Length();
	if (hyperplane_normal_len < FLT_EPSILON)
	{
		outError = "Initial tetrahedron is degenerate";
		return EResult::Degenerate;
	}

	int idx5 = -1;
	float max_hyperplane_dist = 0.0f;
	for (int i = 0; i < (int)mPositions.size(); ++i)
		if (i != idx1 && i != idx2 && i != idx3 && i != idx4)
		{
			float dist = abs((mPositions[i] - mPositions[idx1]).Dot(hyperplane_normal)) / hyperplane_normal_len;
			if (dist > max_hyperplane_dist)
			{
				max_hyperplane_dist = dist;
				idx5 = i;
			}
		}

	if (idx5 < 0 || Square(max_hyperplane_dist) <= Square(cCoplanarSlopFactor) * coplanar_tolerance_sq)
	{
		outError = "All points are co-hyperplanar";
		return EResult::Degenerate;
	}

	// === Step 5: Create pentachoron (5 tetrahedral cells) ===
	// The 5 extreme points
	int pent[5] = { idx1, idx2, idx3, idx4, idx5 };

	// Create 5 cells: cell[c] uses all vertices except pent[c]
	Cell *cells[5];
	for (int c = 0; c < 5; ++c)
	{
		int v[4];
		int k = 0;
		for (int j = 0; j < 5; ++j)
			if (j != c)
				v[k++] = pent[j];

		cells[c] = CreateTetrahedron(v[0], v[1], v[2], v[3]);
	}

	// Compute interior point as centroid of pentachoron (guaranteed inside hull at all times)
	mInteriorPoint = 0.2f * (mPositions[pent[0]] + mPositions[pent[1]] + mPositions[pent[2]] + mPositions[pent[3]] + mPositions[pent[4]]);

	// Orient each cell: normal should point away from the interior
	for (int c = 0; c < 5; ++c)
	{
		if (cells[c]->mNormal.Dot(mInteriorPoint - cells[c]->mCentroid) > 0.0f)
		{
			// Normal points toward the interior; swap two vertices to flip
			std::swap(cells[c]->mIdx[1], cells[c]->mIdx[2]);
			cells[c]->mNormal = -cells[c]->mNormal;
		}
	}

	// Link cells: cell[c] and cell[d] share a ridge (3 vertices = {pent[0..4]} \ {pent[c], pent[d]})
	for (int c = 0; c < 5; ++c)
		for (int d = c + 1; d < 5; ++d)
		{
			// Find which ridge of cell c is opposite the vertex pent[d]
			int ridge_c = -1;
			for (int j = 0; j < 4; ++j)
				if (cells[c]->mIdx[j] == pent[d])
				{
					ridge_c = j;
					break;
				}

			// Find which ridge of cell d is opposite the vertex pent[c]
			int ridge_d = -1;
			for (int j = 0; j < 4; ++j)
				if (cells[d]->mIdx[j] == pent[c])
				{
					ridge_d = j;
					break;
				}

			JPH_ASSERT(ridge_c >= 0 && ridge_d >= 0);
			sLinkCell(cells[c], ridge_c, cells[d], ridge_d);
		}

	// Build initial conflict lists
	Cells initial_cells(mCells.begin(), mCells.end());
	for (int idx = 0; idx < (int)mPositions.size(); ++idx)
	{
		bool is_pent = false;
		for (int p : pent)
			if (idx == p) { is_pent = true; break; }
		if (!is_pent)
			AssignPointToCell(idx, initial_cells, tolerance_sq);
	}

#ifdef JPH_CONVEX_BUILDER_DEBUG
	DrawState(true);
	++mIteration;
#endif

	// Overestimate of vertices used
	int num_vertices_used = 5;

	// === Quickhull loop ===
	for (;;)
	{
		// Find cell with the furthest conflict point
		Cell *cell_with_furthest = nullptr;
		float furthest_dist_sq = 0.0f;
		for (Cell *c : mCells)
			if (!c->mRemoved && c->mFurthestPointDistanceSq > furthest_dist_sq)
			{
				furthest_dist_sq = c->mFurthestPointDistanceSq;
				cell_with_furthest = c;
			}

		int furthest_point_idx;
		if (cell_with_furthest != nullptr)
		{
			furthest_point_idx = cell_with_furthest->mConflictList.back();
			cell_with_furthest->mConflictList.pop_back();
		}
		else if (!mCoplanarList.empty())
		{
			// Try to assign co-hyperplanar points to cells
			CoplanarList coplanar;
			mCoplanarList.swap(coplanar);
			bool added = false;
			for (const Coplanar &cp : coplanar)
				added |= AssignPointToCell(cp.mPositionIdx, mCells, tolerance_sq);

			if (added)
				continue;

			if (mCoplanarList.empty())
				break;

			do
			{
				// Find the vertex furthest from the hull
				CoplanarList::size_type best_idx = 0;
				float best_dist_sq = mCoplanarList.front().mDistanceSq;
				for (CoplanarList::size_type idx = 1; idx < mCoplanarList.size(); ++idx)
				{
					if (mCoplanarList[idx].mDistanceSq > best_dist_sq)
					{
						best_idx = idx;
						best_dist_sq = mCoplanarList[idx].mDistanceSq;
					}
				}

				std::swap(mCoplanarList[best_idx], mCoplanarList.back());
				furthest_point_idx = mCoplanarList.back().mPositionIdx;
				mCoplanarList.pop_back();

				GetCellForPoint(mPositions[furthest_point_idx], mCells, cell_with_furthest, best_dist_sq);
			} while (!mCoplanarList.empty() && cell_with_furthest == nullptr);

			if (cell_with_furthest == nullptr)
				break;
		}
		else
		{
			break;
		}

		// Check vertex limit
		if (num_vertices_used >= inMaxVertices)
		{
			num_vertices_used = GetNumVerticesUsed();
			if (num_vertices_used >= inMaxVertices)
				return EResult::MaxVerticesReached;
		}

		++num_vertices_used;

		// Add the point to the hull
		Cells new_cells;
		AddPoint(cell_with_furthest, furthest_point_idx, new_cells);

		// Redistribute conflict lists from removed cells
		for (const Cell *cell : mCells)
			if (cell->mRemoved)
				for (int idx : cell->mConflictList)
					AssignPointToCell(idx, new_cells, tolerance_sq);

		// Free removed cells
		GarbageCollectCells();

#ifdef JPH_CONVEX_BUILDER_DEBUG
		DrawState(true);
		++mIteration;
#endif
	}

	// Check we have enough cells
	if (mCells.size() < 5)
	{
		outError = "Too few cells in hull";
		return EResult::TooFewFaces;
	}

	return EResult::Success;
}

void ConvexHullBuilder::AddPoint(Cell *inFacingCell, int inIdx, Cells &outNewCells)
{
	Vec4 pos = mPositions[inIdx];

#ifdef JPH_ENABLE_ASSERTS
	ValidateCells();
#endif

	// Find boundary ridges via DFS
	BoundaryRidges ridges;
	FindBoundaryRidges(inFacingCell, pos, ridges);
	JPH_ASSERT(ridges.size() >= 4);

	// Create new cells from boundary ridges + new point
	outNewCells.reserve(ridges.size());
	for (const BoundaryRidge &r : ridges)
	{
		Cell *c = CreateTetrahedron(r.mVertices[0], r.mVertices[1], r.mVertices[2], inIdx);

		// If the cell is degenerate (near-zero normal from coplanar vertices),
		// inherit the normal from the non-facing cell (same hyperplane)
		float normal_len_sq = c->mNormal.LengthSq();
		if (normal_len_sq < 1.0e-20f)
		{
			c->mNormal = r.mNonFacingCell->mNormal;
		}
		else
		{
			// Orient: normal should point away from the interior (toward outside)
			if (c->mNormal.Dot(mInteriorPoint - c->mCentroid) > 0.0f)
			{
				// Normal points toward interior, flip
				std::swap(c->mIdx[1], c->mIdx[2]);
				c->mNormal = -c->mNormal;
			}
		}

		outNewCells.push_back(c);

		// Link to non-facing cell. Ridge 3 (opposite mIdx[3] = new point) connects to it.
		sLinkCell(c, 3, r.mNonFacingCell, r.mNonFacingRidgeIdx);
	}

	// Link new cells to each other
	for (Cells::size_type i = 0; i < outNewCells.size(); ++i)
		for (Cells::size_type k = i + 1; k < outNewCells.size(); ++k)
			sTryLinkCells(outNewCells[i], outNewCells[k]);

#ifdef JPH_ENABLE_ASSERTS
	// Verify all internal ridges are linked
	for (Cell *c : outNewCells)
		for (int j = 0; j < 4; ++j)
			JPH_ASSERT(c->mNeighbour[j] != nullptr);

	ValidateCells();
#endif
}

void ConvexHullBuilder::GarbageCollectCells()
{
	for (int i = (int)mCells.size() - 1; i >= 0; --i)
	{
		Cell *c = mCells[i];
		if (c->mRemoved)
		{
			FreeCell(c);
			mCells.erase(mCells.begin() + i);
		}
	}
}

ConvexHullBuilder::Cell *ConvexHullBuilder::CreateCell()
{
	Cell *c = new Cell();

#ifdef JPH_CONVEX_BUILDER_DEBUG
	c->mIteration = mIteration;
#endif

	mCells.push_back(c);
	return c;
}

ConvexHullBuilder::Cell *ConvexHullBuilder::CreateTetrahedron(int inIdx0, int inIdx1, int inIdx2, int inIdx3)
{
	Cell *c = CreateCell();
	c->Initialize(inIdx0, inIdx1, inIdx2, inIdx3, mPositions.data());
	return c;
}

void ConvexHullBuilder::FreeCell(Cell *inCell)
{
	JPH_ASSERT(inCell->mRemoved);

#ifdef JPH_ENABLE_ASSERTS
	for (int j = 0; j < 4; ++j)
		JPH_ASSERT(inCell->mNeighbour[j] == nullptr);
#endif

	delete inCell;
}

void ConvexHullBuilder::sLinkCell(Cell *inCell1, int inRidge1, Cell *inCell2, int inRidge2)
{
	JPH_ASSERT(inCell1->mNeighbour[inRidge1] == nullptr);
	JPH_ASSERT(inCell2->mNeighbour[inRidge2] == nullptr);
	JPH_ASSERT(inCell1 != inCell2);

	inCell1->mNeighbour[inRidge1] = inCell2;
	inCell1->mNeighbourRidge[inRidge1] = inRidge2;
	inCell2->mNeighbour[inRidge2] = inCell1;
	inCell2->mNeighbourRidge[inRidge2] = inRidge1;
}

void ConvexHullBuilder::sUnlinkCell(Cell *inCell)
{
	for (int j = 0; j < 4; ++j)
	{
		Cell *nb = inCell->mNeighbour[j];
		if (nb != nullptr)
		{
			int nb_ridge = inCell->mNeighbourRidge[j];
			JPH_ASSERT(nb->mNeighbour[nb_ridge] == inCell);
			nb->mNeighbour[nb_ridge] = nullptr;
			inCell->mNeighbour[j] = nullptr;
		}
	}
}

void ConvexHullBuilder::FindBoundaryRidges(Cell *inFacingCell, Vec4Arg inVertex, BoundaryRidges &outRidges) const
{
	JPH_ASSERT(outRidges.empty());
	JPH_ASSERT(inFacingCell->IsFacing(inVertex));

	inFacingCell->mRemoved = true;

	struct StackEntry
	{
		Cell *	mCell;
		int		mNextRidge;		// 0..4, which ridge to examine next
	};

	constexpr int cMaxStackSize = 256;
	StackEntry stack[cMaxStackSize];
	int stack_pos = 0;
	stack[0] = { inFacingCell, 0 };

	while (stack_pos >= 0)
	{
		StackEntry &entry = stack[stack_pos];

		if (entry.mNextRidge >= 4)
		{
			// Done with this cell, unlink it (caller will free later)
			sUnlinkCell(entry.mCell);
			--stack_pos;
			continue;
		}

		int j = entry.mNextRidge++;
		Cell *neighbour = entry.mCell->mNeighbour[j];

		if (neighbour == nullptr || neighbour->mRemoved)
			continue;

		if (neighbour->IsFacing(inVertex))
		{
			// Facing: recurse deeper
			neighbour->mRemoved = true;
			++stack_pos;
			JPH_ASSERT(stack_pos < cMaxStackSize);
			stack[stack_pos] = { neighbour, 0 };
		}
		else
		{
			// Non-facing: this is a boundary ridge
			BoundaryRidge ridge;
			ridge.mNonFacingCell = neighbour;
			ridge.mNonFacingRidgeIdx = entry.mCell->mNeighbourRidge[j];

			// Ridge vertices = cell's vertices minus the one opposite this ridge
			int k = 0;
			for (int i = 0; i < 4; ++i)
				if (i != j)
					ridge.mVertices[k++] = entry.mCell->mIdx[i];

			outRidges.push_back(ridge);
		}
	}
}

bool ConvexHullBuilder::sTryLinkCells(Cell *inCell1, Cell *inCell2)
{
	// Two new cells share a ridge iff they have 3 common vertices.
	// Only check ridges 0,1,2 (ridge 3 connects to non-facing cell and is already linked).
	for (int r1 = 0; r1 < 3; ++r1)
	{
		if (inCell1->mNeighbour[r1] != nullptr)
			continue;

		for (int r2 = 0; r2 < 3; ++r2)
		{
			if (inCell2->mNeighbour[r2] != nullptr)
				continue;

			// Count matching vertices between the two ridges
			// Ridge r1 of cell1: vertices {mIdx[0..3]} \ {mIdx[r1]}
			// Ridge r2 of cell2: vertices {mIdx[0..3]} \ {mIdx[r2]}
			int count = 0;
			for (int a = 0; a < 4; ++a)
			{
				if (a == r1) continue;
				for (int b = 0; b < 4; ++b)
				{
					if (b == r2) continue;
					if (inCell1->mIdx[a] == inCell2->mIdx[b])
					{
						++count;
						break;
					}
				}
			}

			if (count == 3) // All 3 ridge vertices match
			{
				sLinkCell(inCell1, r1, inCell2, r2);
				return true;
			}
		}
	}
	return false;
}

#ifdef JPH_ENABLE_ASSERTS

void ConvexHullBuilder::DumpCell(const Cell *inCell) const
{
	Trace("cell:0x%p idx:{%d,%d,%d,%d} removed:%d", inCell,
		inCell->mIdx[0], inCell->mIdx[1], inCell->mIdx[2], inCell->mIdx[3],
		(int)inCell->mRemoved);
	for (int j = 0; j < 4; ++j)
		Trace("\tridge %d: nb=0x%p nb_ridge=%d", j, inCell->mNeighbour[j], inCell->mNeighbourRidge[j]);
}

void ConvexHullBuilder::DumpCells() const
{
	Trace("Dump Cells:");
	for (const Cell *c : mCells)
		if (!c->mRemoved)
			DumpCell(c);
}

void ConvexHullBuilder::ValidateCell(const Cell *inCell) const
{
	if (inCell->mRemoved)
	{
		for (int j = 0; j < 4; ++j)
			JPH_ASSERT(inCell->mNeighbour[j] == nullptr);
	}
	else
	{
		for (int j = 0; j < 4; ++j)
		{
			// All 4 vertex indices must be distinct
			for (int k = j + 1; k < 4; ++k)
				JPH_ASSERT(inCell->mIdx[j] != inCell->mIdx[k]);

			// Must have a neighbour
			const Cell *nb = inCell->mNeighbour[j];
			JPH_ASSERT(nb != nullptr);
			if (nb != nullptr)
			{
				// Neighbour must not be self
				JPH_ASSERT(nb != inCell);

				// Neighbour's reciprocal ridge must point back to us
				int nb_ridge = inCell->mNeighbourRidge[j];
				JPH_ASSERT(nb->mNeighbour[nb_ridge] == inCell);
				JPH_ASSERT(nb->mNeighbourRidge[nb_ridge] == j);

				// Shared vertices (ridge): cell's vertices minus mIdx[j], neighbour's vertices minus nb.mIdx[nb_ridge]
				// These two sets of 3 vertices must match
				int shared = 0;
				for (int a = 0; a < 4; ++a)
				{
					if (a == j) continue;
					for (int b = 0; b < 4; ++b)
					{
						if (b == nb_ridge) continue;
						if (inCell->mIdx[a] == nb->mIdx[b])
						{
							++shared;
							break;
						}
					}
				}
				JPH_ASSERT(shared == 3);
			}
		}

		// Validate no two neighbours are the same cell
		if (mCells.size() > 5)
			for (int j = 0; j < 4; ++j)
				for (int k = j + 1; k < 4; ++k)
					JPH_ASSERT(inCell->mNeighbour[j] != inCell->mNeighbour[k]);
	}
}

void ConvexHullBuilder::ValidateCells() const
{
	for (const Cell *c : mCells)
		ValidateCell(c);
}

#endif // JPH_ENABLE_ASSERTS

void ConvexHullBuilder::GetCenterOfMassAndVolume(Vec4 &outCenterOfMass, float &outVolume) const
{
	// Reference point = average of cell centroids
	Vec4 ref = Vec4::sZero();
	for (const Cell *c : mCells)
		ref += c->mCentroid;
	ref /= float(mCells.size());

	// Sum signed 4-volumes of pentachoron (cell + ref)
	outVolume = 0.0f;
	outCenterOfMass = Vec4::sZero();
	for (const Cell *c : mCells)
	{
		Vec4 e0 = mPositions[c->mIdx[0]] - ref;
		Vec4 e1 = mPositions[c->mIdx[1]] - ref;
		Vec4 e2 = mPositions[c->mIdx[2]] - ref;
		Vec4 e3 = mPositions[c->mIdx[3]] - ref;

		// Signed 4-volume = det(e0, e1, e2, e3) / 24
		// det = e3.Dot(sCross(e0, e1, e2))
		float det = e3.Dot(Vec4::sCross(e0, e1, e2));

		// Center of mass of this pentachoron = (v0 + v1 + v2 + v3 + ref) / 5
		Vec4 com5 = mPositions[c->mIdx[0]] + mPositions[c->mIdx[1]]
				  + mPositions[c->mIdx[2]] + mPositions[c->mIdx[3]] + ref;

		outVolume += det;
		outCenterOfMass += det * com5;
	}

	// Finalize
	if (abs(outVolume) > FLT_EPSILON)
		outCenterOfMass /= 5.0f * outVolume;
	else
		outCenterOfMass = ref;

	outVolume /= 24.0f;
}

void ConvexHullBuilder::DetermineMaxError(Cell *&outCellWithMaxError, float &outMaxError, int &outMaxErrorPositionIdx, float &outCoplanarDistance) const
{
	outCoplanarDistance = DetermineCoplanarDistance();

	float max_error = 0.0f;
	Cell *max_error_cell = nullptr;
	int max_error_point = -1;

	for (int i = 0; i < (int)mPositions.size(); ++i)
	{
		Vec4 v = mPositions[i];

		float min_ridge_dist_sq = FLT_MAX;
		Cell *min_ridge_dist_cell = nullptr;

		for (Cell *c : mCells)
		{
			if (c->mRemoved) continue;

			float normal_len = c->mNormal.Length();
			JPH_ASSERT(normal_len > 0.0f);
			float plane_dist = c->mNormal.Dot(v - c->mCentroid) / normal_len;

			if (plane_dist > -cCoplanarSlopFactor * outCoplanarDistance)
			{
				float ridge_dist_sq = GetDistanceToRidgeSq(v, c);
				if (ridge_dist_sq < min_ridge_dist_sq)
				{
					min_ridge_dist_sq = ridge_dist_sq;
					min_ridge_dist_cell = c;
				}

				if (ridge_dist_sq == 0.0f && plane_dist > max_error)
				{
					max_error = plane_dist;
					max_error_cell = c;
					max_error_point = i;
				}
			}
		}

		float min_ridge_dist = sqrt(min_ridge_dist_sq);
		if (min_ridge_dist_cell != nullptr && min_ridge_dist > max_error)
		{
			max_error = min_ridge_dist;
			max_error_cell = min_ridge_dist_cell;
			max_error_point = i;
		}
	}

	outCellWithMaxError = max_error_cell;
	outMaxError = max_error;
	outMaxErrorPositionIdx = max_error_point;
}

#ifdef JPH_CONVEX_BUILDER_DEBUG

void ConvexHullBuilder::DrawState([[maybe_unused]] bool inDrawConflictList) const
{
	// Minimal stub — 4D drawing not yet implemented
}

#endif // JPH_CONVEX_BUILDER_DEBUG

#ifdef JPH_CONVEX_BUILDER_DUMP_SHAPE

void ConvexHullBuilder::DumpShape() const
{
	static atomic<int> sShapeNo = 1;
	int shape_no = sShapeNo++;

	std::ofstream f;
	f.open(StringFormat("dumped_shape%d.cpp", shape_no).c_str(), std::ofstream::out | std::ofstream::trunc);
	if (!f.is_open())
		return;

	f << "{\n";
	for (const Vec4 &v : mPositions)
		f << StringFormat("\tVec4(%.9gf, %.9gf, %.9gf, %.9gf),\n",
			(double)v.GetX(), (double)v.GetY(), (double)v.GetZ(), (double)v.GetW());
	f << "},\n";
}

#endif // JPH_CONVEX_BUILDER_DUMP_SHAPE

JPH_NAMESPACE_END
