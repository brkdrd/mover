// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

// Define to validate the integrity of the hull structure
//#define JPH_EPA_CONVEX_BUILDER_VALIDATE

// Define to draw the building of the hull for debugging purposes
//#define JPH_EPA_CONVEX_BUILDER_DRAW

#include <Jolt/Core/NonCopyable.h>
#include <Jolt/Core/BinaryHeap.h>

#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
	#include <Jolt/Renderer/DebugRenderer.h>
	#include <Jolt/Core/StringTools.h>
#endif

JPH_NAMESPACE_BEGIN

/// A 4D convex hull builder specifically made for the EPA penetration depth calculation.
/// Boundary cells are tetrahedra. It trades accuracy for speed and will simply abort if
/// the hull forms defects due to numerical precision problems.
class EPAConvexHullBuilder : public NonCopyable
{
private:
#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
	/// Factor to scale convex hull when debug drawing the construction process
	static constexpr Real cDrawScale = 10;
#endif

public:
	// TODO(4D): these two limits are independent in 4D, and both had to grow a lot.
	//
	// 3D Jolt sets cMaxPoints = cMaxTriangles / 2, which is Euler's formula: a 3D simplicial polytope
	// with V vertices has exactly 2V - 4 triangles. In 4D the analogous bound is quadratic - a
	// simplicial 4-polytope with V vertices can have up to V(V-3)/2 tetrahedral cells - so the ratio
	// does not hold and the cell pool ran out long before the point limit, making AddPoint fail and
	// EPA return whatever it had (measured: it stopped at 47 points and reported a penetration depth
	// of 2.81 where the answer was 3.0).
	//
	// Approximating the boundary of a 3-sphere with a polytope also simply needs more vertices than
	// approximating a 2-sphere, so 4D EPA needs more iterations for the same accuracy on curved
	// shapes. Measured on two overlapping spheres (the worst case; boxes converge almost instantly):
	// 64 points -> 2.93, 96 -> 2.97, 128 -> 2.99, 192 -> 2.998. Cells needed scale at roughly 8 per
	// point for points spread over a 3-sphere; 768 is enough for 96 points (1024 changes nothing).
	//
	// This costs stack: the builder is a local in EPAPenetrationDepth::GetPenetrationDepth and is
	// ~104 KB at these limits. If that ever becomes a problem, lower to 512/64 (~68 KB) and accept
	// ~2.5% error on curved shapes, or move the builder off the stack.
	static constexpr int cMaxCells = 768;					///< Max cells in hull
	static constexpr int cMaxPoints = 96;					///< Max number of points in hull
	static constexpr int cMaxRidgeLength = 256;				///< Max number of horizon ridges
	static constexpr float cMinCellVolume = 1.0e-12f;		///< Minimum Gram determinant before a cell is considered degenerate
	static constexpr float cBarycentricEpsilon = 1.0e-3f;	///< Epsilon for barycentric interior check

	// Forward declare
	class Cell;

	/// Backward-compatible aliases
	using Triangle = Cell;

	/// Information about a horizon ridge collected during FindRidge
	struct HorizonRidge
	{
		Cell *		mNonFacingCell;							///< The non-facing cell on the other side
		int			mNonFacingRidge;						///< Which ridge of that cell connects back
		int			mVertex[3];								///< The 3 vertex indices of this ridge
	};

	using HorizonRidges = StaticArray<HorizonRidge, cMaxRidgeLength>;
	using NewCells = StaticArray<Cell *, cMaxRidgeLength>;
	using NewTriangles = NewCells;

	/// A tetrahedral boundary cell of the 4D EPA convex hull.
	/// Each cell has 4 vertices, 4 neighbouring cells (one per ridge).
	class Cell : public NonCopyable
	{
	public:
		/// Constructor
		inline			Cell(int inIdx0, int inIdx1, int inIdx2, int inIdx3, const Vec4 *inPositions);

		/// Check if cell is facing inPosition
		inline bool		IsFacing(Vec4Arg inPosition) const
		{
			JPH_ASSERT(!mRemoved);
			return mNormal.Dot(inPosition - mCentroid) > 0.0f;
		}

		/// Check if cell is facing the origin
		inline bool		IsFacingOrigin() const
		{
			JPH_ASSERT(!mRemoved);
			return mNormal.Dot(mCentroid) < 0.0f;
		}

		int				mIdx[4];							///< Vertex indices
		Cell *			mNeighbour[4];						///< mNeighbour[i] shares the ridge opposite vertex mIdx[i]
		int				mNeighbourRidge[4];					///< Which ridge of the neighbour connects back to this cell
		Vec4			mNormal;							///< Normal of this cell = sCross of its 3 edge vectors (not normalized), oriented by the vertex winding
		Vec4			mCentroid;							///< Center of the cell
		float			mClosestLenSq = FLT_MAX;			///< Closest distance^2 from origin to cell
		float			mLambda[3];							///< Barycentric coordinates of closest point to origin
		bool			mLambdaRelativeTo0;					///< Reference vertex: true = y0, false = y1
		bool			mClosestPointInterior = false;		///< Closest point is in the interior of the cell
		bool			mRemoved = false;					///< Cell has been removed
		bool			mInQueue = false;					///< Cell was placed in the priority queue
#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
		int				mIteration;							///< Iteration that created this cell
#endif
	};

	/// Factory that creates cells in a fixed size buffer
	class CellFactory : public NonCopyable
	{
	private:
		union alignas(Cell) Block
		{
			uint8		mData[sizeof(Cell)];
			Block *		mNextFree;
		};

		Block			mCells[cMaxCells];
		Block *			mNextFree = nullptr;
		int				mHighWatermark = 0;

	public:
		void			Clear()
		{
			mNextFree = nullptr;
			mHighWatermark = 0;
		}

		Cell *			CreateCell(int inIdx0, int inIdx1, int inIdx2, int inIdx3, const Vec4 *inPositions)
		{
			Cell *c;
			if (mNextFree != nullptr)
			{
				c = reinterpret_cast<Cell *>(&mNextFree->mData);
				mNextFree = mNextFree->mNextFree;
			}
			else
			{
				if (mHighWatermark >= cMaxCells)
					return nullptr;
				c = reinterpret_cast<Cell *>(&mCells[mHighWatermark].mData);
				++mHighWatermark;
			}
			new (c) Cell(inIdx0, inIdx1, inIdx2, inIdx3, inPositions);
			return c;
		}

		void			FreeCell(Cell *inC)
		{
			inC->~Cell();
#ifdef JPH_DEBUG
			memset(inC, 0xcd, sizeof(Cell));
#endif
			Block *bu = reinterpret_cast<Block *>(inC);
			bu->mNextFree = mNextFree;
			mNextFree = bu;
		}
	};

	// Typedefs
	using PointsBase = StaticArray<Vec4, cMaxPoints>;
	using Cells = StaticArray<Cell *, cMaxCells>;
	using Triangles = Cells;

	/// Specialized points list that allows direct access to the size
	class Points : public PointsBase
	{
	public:
		size_type &		GetSizeRef()
		{
			return mSize;
		}
	};

	/// Specialized cell list that keeps them sorted on closest distance to origin
	class CellQueue : public Cells
	{
	public:
		static bool		sCellSorter(const Cell *inC1, const Cell *inC2)
		{
			return inC1->mClosestLenSq > inC2->mClosestLenSq;
		}

		void			push_back(Cell *inC)
		{
			Cells::push_back(inC);
			inC->mInQueue = true;
			BinaryHeapPush(begin(), end(), sCellSorter);
		}

		Cell *			PeekClosest()
		{
			return front();
		}

		Cell *			PopClosest()
		{
			BinaryHeapPop(begin(), end(), sCellSorter);
			Cell *c = back();
			pop_back();
			return c;
		}
	};

	// Backward-compatible aliases
	using TriangleQueue = CellQueue;

	/// Constructor
	explicit			EPAConvexHullBuilder(const Points &inPositions) :
		mPositions(inPositions)
	{
#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
		mIteration = 0;
		mOffset = RVec3::sZero();
#endif
	}

	/// Initialize the hull with 5 points (pentachoron)
	void				Initialize(int inIdx0, int inIdx1, int inIdx2, int inIdx3, int inIdx4)
	{
		mFactory.Clear();
		mCellQueue.clear();

		int idx[5] = { inIdx0, inIdx1, inIdx2, inIdx3, inIdx4 };

		// Centroid of the seed simplex: a point inside the hull, used to orient every cell normal
		// outward. The hull only grows from here, so it stays interior.
		mInteriorPoint = (mPositions[inIdx0] + mPositions[inIdx1] + mPositions[inIdx2] + mPositions[inIdx3] + mPositions[inIdx4]) * 0.2f;

		// Create 5 cells — cell i omits vertex idx[i]
		Cell *cells[5];
		for (int i = 0; i < 5; ++i)
		{
			int v[4];
			int k = 0;
			for (int j = 0; j < 5; ++j)
				if (j != i)
					v[k++] = idx[j];
			cells[i] = CreateCell(v[0], v[1], v[2], v[3]);
		}

		// Link all pairs of cells using brute-force 3-vertex matching
		for (int i = 0; i < 5; ++i)
			for (int j = i + 1; j < 5; ++j)
				sTryLinkCells(cells[i], cells[j]);

		// Add all cells to the priority queue
		for (int i = 0; i < 5; ++i)
			mCellQueue.push_back(cells[i]);

#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
		DrawState();
		++mIteration;
#endif
	}

	/// Check if there's another cell to process
	bool				HasNextTriangle() const
	{
		return !mCellQueue.empty();
	}

	/// Access to the next closest cell (won't remove from queue)
	Cell *				PeekClosestTriangleInQueue()
	{
		return mCellQueue.PeekClosest();
	}

	/// Access to the next closest cell and remove from queue
	Cell *				PopClosestTriangleFromQueue()
	{
		return mCellQueue.PopClosest();
	}

	/// Find the cell for which inPosition is furthest to the front
	Cell *				FindFacingTriangle(Vec4Arg inPosition, float &outBestDistSq)
	{
		Cell *best = nullptr;
		float best_dist_sq = 0.0f;

		for (Cell *c : mCellQueue)
			if (!c->mRemoved)
			{
				float dot = c->mNormal.Dot(inPosition - c->mCentroid);
				if (dot > 0.0f)
				{
					float dist_sq = dot * dot / c->mNormal.LengthSq();
					if (dist_sq > best_dist_sq)
					{
						best = c;
						best_dist_sq = dist_sq;
					}
				}
			}

		outBestDistSq = best_dist_sq;
		return best;
	}

	/// Add a new point to the convex hull
	bool				AddPoint(Cell *inFacingCell, int inIdx, float inClosestDistSq, NewCells &outCells)
	{
		Vec4 pos = mPositions[inIdx];

#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
		DrawMarker(pos, Color::sYellow, 1.0f);
#endif

#ifdef JPH_EPA_CONVEX_BUILDER_VALIDATE
		ValidateCells();
#endif

		// Find horizon ridges
		HorizonRidges ridges;
		if (!FindRidge(inFacingCell, pos, ridges))
			return false;

		// Create new cells from horizon ridges + new point
		for (const HorizonRidge &r : ridges)
		{
			Cell *nc = CreateCell(r.mVertex[0], r.mVertex[1], r.mVertex[2], inIdx);
			if (nc == nullptr)
				return false;
			outCells.push_back(nc);

			// Check if we need to put this cell in the priority queue
			if ((nc->mClosestPointInterior && nc->mClosestLenSq < inClosestDistSq)
				|| nc->mClosestLenSq < 0.0f)
				mCellQueue.push_back(nc);

			// Link the new cell's ridge 3 (opposite the new point) to the non-facing cell
			sLinkCell(nc, 3, r.mNonFacingCell, r.mNonFacingRidge);
		}

		// Link new cells to each other (ridges 0-2)
		for (int i = 0; i < (int)outCells.size(); ++i)
			for (int j = i + 1; j < (int)outCells.size(); ++j)
				sTryLinkCells(outCells[i], outCells[j]);

#ifdef JPH_EPA_CONVEX_BUILDER_VALIDATE
		ValidateCells();
#endif

#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
		DrawState();
		++mIteration;
#endif

		return true;
	}

	/// Free a cell
	void				FreeTriangle(Cell *inC)
	{
#ifdef JPH_ENABLE_ASSERTS
		JPH_ASSERT(inC->mRemoved);
		for (int i = 0; i < 4; ++i)
			JPH_ASSERT(inC->mNeighbour[i] == nullptr);
#endif

#if defined(JPH_EPA_CONVEX_BUILDER_VALIDATE) || defined(JPH_EPA_CONVEX_BUILDER_DRAW)
		Cells::iterator it = std::find(mAllCells.begin(), mAllCells.end(), inC);
		JPH_ASSERT(it != mAllCells.end());
		mAllCells.erase(it);
#endif

		mFactory.FreeCell(inC);
	}

private:
	/// Create a new cell
	Cell *				CreateCell(int inIdx0, int inIdx1, int inIdx2, int inIdx3)
	{
		Cell *c = mFactory.CreateCell(inIdx0, inIdx1, inIdx2, inIdx3, mPositions.data());
		if (c == nullptr)
			return nullptr;

		// TODO(4D): orient the cell normal outward. The normal's sign comes from the vertex winding
		// (sCross of the edge vectors), and neither Initialize (which just omits one vertex at a time)
		// nor the brute-force ridge linking maintains a consistent winding across cells, so roughly
		// half of them would point inward. That breaks IsFacing()/IsFacingOrigin() and makes EPA bail
		// out immediately. mInteriorPoint is a point known to be inside the hull, so "outward" is
		// simply "away from it"; the hull only ever grows, so it stays interior.
		if (c->mNormal.Dot(c->mCentroid - mInteriorPoint) < 0.0f)
		{
			c->mNormal = -c->mNormal;
			c->mClosestLenSq = -c->mClosestLenSq;	// signed distance flips with the normal
		}

#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
		c->mIteration = mIteration;
#endif

#if defined(JPH_EPA_CONVEX_BUILDER_VALIDATE) || defined(JPH_EPA_CONVEX_BUILDER_DRAW)
		mAllCells.push_back(c);
#endif

		return c;
	}

	/// Link two cells at specified ridges
	static void			sLinkCell(Cell *inC1, int inR1, Cell *inC2, int inR2)
	{
		JPH_ASSERT(inR1 >= 0 && inR1 < 4);
		JPH_ASSERT(inR2 >= 0 && inR2 < 4);
		JPH_ASSERT(inC1->mNeighbour[inR1] == nullptr);
		JPH_ASSERT(inC2->mNeighbour[inR2] == nullptr);

		inC1->mNeighbour[inR1] = inC2;
		inC1->mNeighbourRidge[inR1] = inR2;
		inC2->mNeighbour[inR2] = inC1;
		inC2->mNeighbourRidge[inR2] = inR1;
	}

	/// Unlink a cell from its neighbours
	void				UnlinkCell(Cell *inC)
	{
		for (int i = 0; i < 4; ++i)
		{
			Cell *nb = inC->mNeighbour[i];
			if (nb != nullptr)
			{
				int nb_ridge = inC->mNeighbourRidge[i];
				JPH_ASSERT(nb->mNeighbour[nb_ridge] == inC);
				nb->mNeighbour[nb_ridge] = nullptr;
				inC->mNeighbour[i] = nullptr;
			}
		}

		// If not in priority queue, we can free it now
		if (!inC->mInQueue)
			FreeTriangle(inC);
	}

	/// Given a facing cell, find all horizon ridges via DFS.
	/// Flags facing cells for removal.
	bool				FindRidge(Cell *inFacingCell, Vec4Arg inVertex, HorizonRidges &outRidges)
	{
		JPH_ASSERT(outRidges.empty());
		JPH_ASSERT(inFacingCell->IsFacing(inVertex));

		inFacingCell->mRemoved = true;

		struct StackEntry
		{
			Cell *		mCell;
			int			mNextRidge;
		};

		StackEntry stack[cMaxRidgeLength];
		int stack_pos = 0;
		stack[0] = { inFacingCell, 0 };

		while (stack_pos >= 0)
		{
			StackEntry &entry = stack[stack_pos];

			if (entry.mNextRidge >= 4)
			{
				// Done with this cell, unlink it
				UnlinkCell(entry.mCell);
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
				JPH_ASSERT(stack_pos < cMaxRidgeLength);
				stack[stack_pos] = { neighbour, 0 };
			}
			else
			{
				// Non-facing: this is a horizon ridge
				HorizonRidge ridge;
				ridge.mNonFacingCell = neighbour;
				ridge.mNonFacingRidge = entry.mCell->mNeighbourRidge[j];

				// Ridge vertices = cell's vertices minus the one opposite this ridge
				int k = 0;
				for (int i = 0; i < 4; ++i)
					if (i != j)
						ridge.mVertex[k++] = entry.mCell->mIdx[i];

				outRidges.push_back(ridge);
			}
		}

#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
		DrawState();
#endif

		// Need at least 4 horizon ridges to form a closed surface
		return outRidges.size() >= 4;
	}

	/// Try to link two cells if they share a ridge (3 common vertices)
	static bool			sTryLinkCells(Cell *inC1, Cell *inC2)
	{
		// Only check ridges that aren't already linked
		for (int r1 = 0; r1 < 4; ++r1)
		{
			if (inC1->mNeighbour[r1] != nullptr)
				continue;

			for (int r2 = 0; r2 < 4; ++r2)
			{
				if (inC2->mNeighbour[r2] != nullptr)
					continue;

				// Count matching vertices between the two ridges
				int count = 0;
				for (int a = 0; a < 4; ++a)
				{
					if (a == r1) continue;
					for (int b = 0; b < 4; ++b)
					{
						if (b == r2) continue;
						if (inC1->mIdx[a] == inC2->mIdx[b])
						{
							++count;
							break;
						}
					}
				}

				if (count == 3)
				{
					sLinkCell(inC1, r1, inC2, r2);
					return true;
				}
			}
		}
		return false;
	}

#ifdef JPH_EPA_CONVEX_BUILDER_VALIDATE
	void				ValidateCell(const Cell *inC) const
	{
		if (inC->mRemoved)
		{
			for (int i = 0; i < 4; ++i)
				JPH_ASSERT(inC->mNeighbour[i] == nullptr);
		}
		else
		{
			for (int i = 0; i < 4; ++i)
			{
				const Cell *nb = inC->mNeighbour[i];
				JPH_ASSERT(nb != nullptr);
				if (nb != nullptr)
				{
					int nb_ridge = inC->mNeighbourRidge[i];
					JPH_ASSERT(nb->mNeighbour[nb_ridge] == inC);
					JPH_ASSERT(nb->mNeighbourRidge[nb_ridge] == i);
				}
			}
		}
	}

	void				ValidateCells() const
	{
		for (const Cell *c : mAllCells)
			ValidateCell(c);
	}
#endif

#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
public:
	void				DrawState()
	{
		// 4D hull drawing is not directly visualizable; placeholder
	}

	void				DrawLabel(const string_view &inText)
	{
	}

	void				DrawGeometry(const DebugRenderer::GeometryRef &inGeometry, ColorArg inColor)
	{
	}

	void				DrawWireTriangle(const Cell &inCell, ColorArg inColor)
	{
	}

	void				DrawMarker(Vec4Arg inPosition, ColorArg inColor, float inSize)
	{
	}

	void				DrawArrow(Vec4Arg inFrom, Vec4Arg inTo, ColorArg inColor, float inArrowSize)
	{
	}
#endif

private:
	CellFactory			mFactory;							///< Factory to create/free cells
	Vec4				mInteriorPoint { Vec4::sZero() };	///< A point known to be inside the hull, used to orient cell normals outward
	const Points &		mPositions;							///< List of positions
	CellQueue			mCellQueue;							///< Priority queue of cells sorted by closest distance to origin

#if defined(JPH_EPA_CONVEX_BUILDER_VALIDATE) || defined(JPH_EPA_CONVEX_BUILDER_DRAW)
	Cells				mAllCells;							///< All cells for debugging
#endif

#ifdef JPH_EPA_CONVEX_BUILDER_DRAW
	int					mIteration;
	RVec3				mOffset;
#endif
};

// Disable fmadd for numerical precision in the Cell constructor
JPH_PRECISE_MATH_ON

EPAConvexHullBuilder::Cell::Cell(int inIdx0, int inIdx1, int inIdx2, int inIdx3, const Vec4 *inPositions)
{
	// Fill in indices
	mIdx[0] = inIdx0;
	mIdx[1] = inIdx1;
	mIdx[2] = inIdx2;
	mIdx[3] = inIdx3;

	// Clear links
	for (int i = 0; i < 4; ++i)
	{
		mNeighbour[i] = nullptr;
		mNeighbourRidge[i] = 0;
	}

	// Get vertex positions
	Vec4 y0 = inPositions[inIdx0];
	Vec4 y1 = inPositions[inIdx1];
	Vec4 y2 = inPositions[inIdx2];
	Vec4 y3 = inPositions[inIdx3];

	// Calculate centroid
	mCentroid = (y0 + y1 + y2 + y3) * 0.25f;

	// Edge vectors from y0
	Vec4 e1 = y1 - y0;
	Vec4 e2 = y2 - y0;
	Vec4 e3 = y3 - y0;

	// TODO(4D): the normal of a tetrahedral cell is the generalized cross product of its three edge
	// vectors, i.e. the vector orthogonal to the cell's 3-flat, with its orientation coming from the
	// vertex winding exactly as 3D Jolt's y10.Cross(y20) does.
	//
	// It must NOT be "the closest point on the cell's hyperplane to the origin", which is what this
	// port used to store. That collapses to the zero vector whenever the origin lies on the cell's
	// hyperplane, leaving EPA with no direction to expand in. Worse, it silently disabled two whole
	// mechanisms: centroid . closest_point = |closest_point|^2 >= 0 always, so IsFacingOrigin() could
	// never return true and mClosestLenSq could never go negative - and those are exactly the signals
	// that drive the "encapsulate the origin" phase and the hull-defect check.
	mNormal = Vec4::sCross(e1, e2, e3);
	float normal_len_sq = mNormal.LengthSq();

	if (normal_len_sq <= cMinCellVolume)
	{
		// Degenerate cell. Leave mClosestLenSq at its FLT_MAX default so it sorts last, and make sure
		// nothing downstream can read an uninitialized normal or lambda (this was bug #7 in CLAUDE.md:
		// the old code linked the cell while leaving both uninitialized).
		mNormal = Vec4::sZero();
		mLambda[0] = mLambda[1] = mLambda[2] = 0.0f;
		mLambdaRelativeTo0 = true;
		return;
	}

	// Signed squared distance from the origin to the cell's hyperplane. The sign records which side
	// of the cell the origin is on, which is what tells EPA the origin is not enclosed yet.
	float c_dot_n = mCentroid.Dot(mNormal);
	mClosestLenSq = abs(c_dot_n) * c_dot_n / normal_len_sq;

	// Barycentric coordinates of the closest point, by solving the 3x3 Gram system. This is only used
	// to reconstruct the closest point later; the distance above is computed from the normal because
	// that is much more accurate (see the same note in 3D Jolt and in ClosestPoint::GetClosestPointOnTriangle).
	float d00 = e1.Dot(e1), d01 = e1.Dot(e2), d02 = e1.Dot(e3);
	float d11 = e2.Dot(e2), d12 = e2.Dot(e3), d22 = e3.Dot(e3);

	float C00 = d11 * d22 - d12 * d12;
	float C01 = d02 * d12 - d01 * d22;
	float C02 = d01 * d12 - d02 * d11;
	float C11 = d00 * d22 - d02 * d02;
	float C12 = d01 * d02 - d00 * d12;
	float C22 = d00 * d11 - d01 * d01;

	// The Gram determinant of 3 vectors is exactly the squared 3-volume they span, i.e. |sCross|^2
	float determinant = normal_len_sq;

	float b0 = -y0.Dot(e1);
	float b1 = -y0.Dot(e2);
	float b2 = -y0.Dot(e3);

	float l0 = (b0 * C00 + b1 * C01 + b2 * C02) / determinant;
	float l1 = (b0 * C01 + b1 * C11 + b2 * C12) / determinant;
	float l2 = (b0 * C02 + b1 * C12 + b2 * C22) / determinant;
	mLambda[0] = l0;
	mLambda[1] = l1;
	mLambda[2] = l2;
	mLambdaRelativeTo0 = true;

	// Check if the closest point is interior to the tetrahedron. As in 3D, several coplanar cells can
	// share the same hyperplane and only one of them has the origin's projection as an interior point;
	// that is the one we want to use for the contact points.
	if (l0 > -cBarycentricEpsilon && l1 > -cBarycentricEpsilon && l2 > -cBarycentricEpsilon
		&& l0 + l1 + l2 < 1.0f + cBarycentricEpsilon)
		mClosestPointInterior = true;
}

JPH_PRECISE_MATH_OFF

JPH_NAMESPACE_END
