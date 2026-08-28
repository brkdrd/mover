// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

//#define JPH_CONVEX_BUILDER_DEBUG
//#define JPH_CONVEX_BUILDER_DUMP_SHAPE

#ifdef JPH_CONVEX_BUILDER_DEBUG
	#include <Jolt/Core/Color.h>
#endif

#include <Jolt/Core/StaticArray.h>
#include <Jolt/Core/NonCopyable.h>

JPH_NAMESPACE_BEGIN

/// A convex hull builder that constructs 4D convex hulls using the quickhull algorithm.
/// Boundary cells are tetrahedra. Used for offline processing.
class JPH_EXPORT ConvexHullBuilder : public NonCopyable
{
public:
	// Forward declare
	class Cell;

	/// Conflict list: indices into mPositions
	using ConflictList = Array<int>;

	/// A tetrahedral boundary cell of the 4D convex hull.
	/// Each cell has 4 vertices, 4 neighbouring cells (one per ridge), and a hyperplane normal.
	class Cell : public NonCopyable
	{
	public:
		JPH_OVERRIDE_NEW_DELETE

		/// Initialize a cell with four vertex indices
		void					Initialize(int inIdx0, int inIdx1, int inIdx2, int inIdx3, const Vec4 *inPositions);

		/// Calculates the normal and centroid for this cell
		void					CalculateNormalAndCentroid(const Vec4 *inPositions);

		/// Check if inPosition is on the front side of this cell's hyperplane
		inline bool				IsFacing(Vec4Arg inPosition) const
		{
			JPH_ASSERT(!mRemoved);
			return mNormal.Dot(inPosition - mCentroid) > 0.0f;
		}

		int						mIdx[4];								///< Vertex indices (into mPositions)
		Cell *					mNeighbour[4] = { nullptr, nullptr, nullptr, nullptr };	///< Neighbour cell for each ridge. mNeighbour[i] shares the ridge opposite vertex mIdx[i].
		int						mNeighbourRidge[4] = { 0, 0, 0, 0 };	///< Which ridge of the neighbour connects back to this cell
		Vec4					mNormal;								///< Outward-pointing normal of the cell's hyperplane (length = 6x cell volume)
		Vec4					mCentroid;								///< Center of the cell
		ConflictList			mConflictList;							///< Positions assigned to this cell
		float					mFurthestPointDistanceSq = 0.0f;		///< Squared distance of the furthest conflict point to this cell's hyperplane
		bool					mRemoved = false;						///< Cell has been removed from the hull
#ifdef JPH_CONVEX_BUILDER_DEBUG
		int						mIteration;								///< Iteration that this cell was created
#endif
	};

	// Backward-compatible aliases
	using Face = Cell;

	// Typedefs
	using Positions = Array<Vec4>;
	using Cells = Array<Cell *>;
	using Faces = Cells;

	/// Constructor
	explicit					ConvexHullBuilder(const Positions &inPositions);

	/// Destructor
								~ConvexHullBuilder()					{ FreeCells(); }

	/// Result enum that indicates how the hull got created
	enum class EResult
	{
		Success,											///< Hull building finished successfully
		MaxVerticesReached,									///< Hull building finished, but desired accuracy not reached
		TooFewPoints,										///< Too few points to create a hull
		TooFewFaces,										///< Too few cells in the created hull
		Degenerate,											///< Degenerate hull detected
	};

	/// Takes all positions and builds a 4D convex hull.
	/// Points closer than inTolerance to the hull are discarded.
	/// @param inMaxVertices Max vertices to allow. Specify INT_MAX for no limit.
	/// @param inTolerance Max allowed distance of a point outside the hull
	/// @param outError Error message on failure
	/// @return Status code
	EResult						Initialize(int inMaxVertices, float inTolerance, const char *&outError);

	/// Returns the number of vertices currently used by the hull
	int							GetNumVerticesUsed() const;

	/// Returns true if the hull contains a cell with the given vertex indices (set equality, order-independent)
	bool						ContainsCell(const Array<int> &inIndices) const;

	/// Backward-compatible alias
	bool						ContainsFace(const Array<int> &inIndices) const { return ContainsCell(inIndices); }

	/// Calculate the center of mass and 4-volume of the current convex hull
	void						GetCenterOfMassAndVolume(Vec4 &outCenterOfMass, float &outVolume) const;

	/// Determines the point furthest outside the hull and reports the error distance
	void						DetermineMaxError(Cell *&outCellWithMaxError, float &outMaxError, int &outMaxErrorPositionIdx, float &outCoplanarDistance) const;

	/// Access to the created cells
	const Cells &				GetCells() const						{ return mCells; }

	/// Backward-compatible alias
	const Cells &				GetFaces() const						{ return mCells; }

private:
	/// Minimum Gram determinant for a non-degenerate triangle (used in initialization)
	static constexpr float		cMinTriangleGramSq = 1.0e-12f;

	/// Extra slop factor for coplanar detection
	static constexpr float		cCoplanarSlopFactor = 6.0f;

	/// A boundary ridge found during AddPoint's DFS
	struct BoundaryRidge
	{
		Cell *					mNonFacingCell;							///< The neighbouring non-facing cell
		int						mNonFacingRidgeIdx;						///< Ridge index in the non-facing cell that connects back
		int						mVertices[3];							///< The 3 vertices of this boundary ridge
	};

	using BoundaryRidges = Array<BoundaryRidge>;

	/// Determine tolerance for detecting co-hyperplanar points
	float						DetermineCoplanarDistance() const;

	/// Find the cell for which inPoint is furthest to the front
	void						GetCellForPoint(Vec4Arg inPoint, const Cells &inCells, Cell *&outCell, float &outDistSq) const;

	/// Distance from inPoint to the nearest ridge boundary of inCell.
	/// Returns 0 if point projects inside the tetrahedron.
	float						GetDistanceToRidgeSq(Vec4Arg inPoint, const Cell *inCell) const;

	/// Assign a point to the best-matching cell. Returns true if assigned to a conflict list.
	bool						AssignPointToCell(int inPositionIdx, const Cells &inCells, float inToleranceSq);

	/// Add a new point to the convex hull
	void						AddPoint(Cell *inFacingCell, int inIdx, Cells &outNewCells);

	/// Remove cells flagged as removed
	void						GarbageCollectCells();

	/// Create a new empty cell
	Cell *						CreateCell();

	/// Create a new tetrahedral cell
	Cell *						CreateTetrahedron(int inIdx0, int inIdx1, int inIdx2, int inIdx3);

	/// Delete a cell
	void						FreeCell(Cell *inCell);

	/// Release all cells
	void						FreeCells();

	/// Link two cells at the given ridges
	static void					sLinkCell(Cell *inCell1, int inRidge1, Cell *inCell2, int inRidge2);

	/// Unlink a cell from all its neighbours
	static void					sUnlinkCell(Cell *inCell);

	/// DFS to find boundary ridges between facing and non-facing cells
	void						FindBoundaryRidges(Cell *inFacingCell, Vec4Arg inVertex, BoundaryRidges &outRidges) const;

	/// Try to link two cells if they share a common ridge
	static bool					sTryLinkCells(Cell *inCell1, Cell *inCell2);

#ifdef JPH_ENABLE_ASSERTS
	void						DumpCell(const Cell *inCell) const;
	void						DumpCells() const;
	void						ValidateCell(const Cell *inCell) const;
	void						ValidateCells() const;
#endif

#ifdef JPH_CONVEX_BUILDER_DEBUG
	void						DrawState(bool inDrawConflictList = false) const;
#endif

#ifdef JPH_CONVEX_BUILDER_DUMP_SHAPE
	void						DumpShape() const;
#endif

	const Positions &			mPositions;								///< Input positions
	Cells						mCells;									///< Boundary cells of the hull
	Vec4						mInteriorPoint = Vec4::sZero();			///< A point guaranteed to be inside the hull (centroid of initial simplex)

	struct Coplanar
	{
		int						mPositionIdx;							///< Index in mPositions
		float					mDistanceSq;							///< Distance to closest cell boundary
	};
	using CoplanarList = Array<Coplanar>;
	CoplanarList				mCoplanarList;							///< Co-hyperplanar points to be added later

#ifdef JPH_CONVEX_BUILDER_DEBUG
	int							mIteration;
#endif
};

JPH_NAMESPACE_END
