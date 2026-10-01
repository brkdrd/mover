// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Geometry/Triangle.h>
#include <Jolt/Geometry/IndexedTriangle.h>
#include <Jolt/Geometry/Plane.h>
#include <Jolt/Math/Mat44.h>
#include <Jolt/Math/RMat44.h> // Transformed(RMat44Arg) below dereferences the full type

JPH_NAMESPACE_BEGIN

/// Axis aligned box (4D)
class [[nodiscard]] AABox
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructor
					AABox()												: mMin(Vec4::sReplicate(FLT_MAX)), mMax(Vec4::sReplicate(-FLT_MAX)) { }
					AABox(Vec4Arg inMin, Vec4Arg inMax)					: mMin(inMin), mMax(inMax) { }
					AABox(DVec4Arg inMin, DVec4Arg inMax)				: mMin(Vec4(float(inMin.GetX()), float(inMin.GetY()), float(inMin.GetZ()), float(inMin.GetW()))), mMax(Vec4(float(inMax.GetX()), float(inMax.GetY()), float(inMax.GetZ()), float(inMax.GetW()))) { }
					AABox(Vec4Arg inCenter, float inRadius)				: mMin(inCenter - Vec4::sReplicate(inRadius)), mMax(inCenter + Vec4::sReplicate(inRadius)) { }

	/// Create box from 2 points
	static AABox	sFromTwoPoints(Vec4Arg inP1, Vec4Arg inP2)			{ return AABox(Vec4::sMin(inP1, inP2), Vec4::sMax(inP1, inP2)); }

	/// Create box from indexed triangle
	static AABox	sFromTriangle(const VertexList &inVertices, const IndexedTriangle &inTriangle)
	{
		AABox box = sFromTwoPoints(Vec4(inVertices[inTriangle.mIdx[0]]), Vec4(inVertices[inTriangle.mIdx[1]]));
		box.Encapsulate(Vec4(inVertices[inTriangle.mIdx[2]]));
		return box;
	}

	/// Get bounding box of size FLT_MAX
	static AABox	sBiggest()
	{
		/// Max half extent of AABox is 0.5 * FLT_MAX so that GetSize() remains finite
		return AABox(Vec4::sReplicate(-0.5f * FLT_MAX), Vec4::sReplicate(0.5f * FLT_MAX));
	}

	/// Comparison operators
	bool			operator == (const AABox &inRHS) const				{ return mMin == inRHS.mMin && mMax == inRHS.mMax; }
	bool			operator != (const AABox &inRHS) const				{ return mMin != inRHS.mMin || mMax != inRHS.mMax; }

	/// Reset the bounding box to an empty bounding box
	void			SetEmpty()
	{
		mMin = Vec4::sReplicate(FLT_MAX);
		mMax = Vec4::sReplicate(-FLT_MAX);
	}

	/// Check if the bounding box is valid (max >= min)
	bool			IsValid() const
	{
		return UVec4::sAnd(Vec4::sLessOrEqual(mMin, mMax), Vec4::sLessOrEqual(mMin, mMax)).TestAllTrue();
	}

	/// Encapsulate point in bounding box
	void			Encapsulate(Vec4Arg inPos)
	{
		mMin = Vec4::sMin(mMin, inPos);
		mMax = Vec4::sMax(mMax, inPos);
	}

	/// Encapsulate bounding box in bounding box
	void			Encapsulate(const AABox &inRHS)
	{
		mMin = Vec4::sMin(mMin, inRHS.mMin);
		mMax = Vec4::sMax(mMax, inRHS.mMax);
	}

	/// Encapsulate triangle in bounding box
	void			Encapsulate(const Triangle &inRHS)
	{
		Encapsulate(Vec4(inRHS.mV[0]));
		Encapsulate(Vec4(inRHS.mV[1]));
		Encapsulate(Vec4(inRHS.mV[2]));
	}

	/// Encapsulate triangle in bounding box
	void			Encapsulate(const VertexList &inVertices, const IndexedTriangle &inTriangle)
	{
		for (uint32 idx : inTriangle.mIdx)
			Encapsulate(Vec4(inVertices[idx]));
	}

	/// Intersect this bounding box with inOther, returns the intersection
	AABox			Intersect(const AABox &inOther) const
	{
		return AABox(Vec4::sMax(mMin, inOther.mMin), Vec4::sMin(mMax, inOther.mMax));
	}

	/// Make sure that each edge of the bounding box has a minimal length
	void			EnsureMinimalEdgeLength(float inMinEdgeLength)
	{
		Vec4 min_length = Vec4::sReplicate(inMinEdgeLength);
		mMax = Vec4::sSelect(mMax, mMin + min_length, Vec4::sLess(mMax - mMin, min_length));
	}

	/// Widen the box on both sides by inVector
	void			ExpandBy(Vec4Arg inVector)
	{
		mMin -= inVector;
		mMax += inVector;
	}

	/// Get center of bounding box
	Vec4			GetCenter() const
	{
		return 0.5f * (mMin + mMax);
	}

	/// Get extent of bounding box (half of the size)
	Vec4			GetExtent() const
	{
		return 0.5f * (mMax - mMin);
	}

	/// Get size of bounding box
	Vec4			GetSize() const
	{
		return mMax - mMin;
	}

	/// Get hyper-surface area of 4D bounding box: 2*(xy*xz + xy*xw + xz*xw + yz*yw + yz*xw + yw*xz) ... actually 8 * sum of products of 3 extents
	/// In 4D, the "surface" consists of 8 cubic cells. The hyper-surface area is 2*(ex*ey*ez + ex*ey*ew + ex*ez*ew + ey*ez*ew)
	float			GetSurfaceArea() const
	{
		Vec4 extent = mMax - mMin;
		float ex = extent.GetX(), ey = extent.GetY(), ez = extent.GetZ(), ew = extent.GetW();
		return 2.0f * (ex * ey * ez + ex * ey * ew + ex * ez * ew + ey * ez * ew);
	}

	/// Get hypervolume of 4D bounding box
	float			GetVolume() const
	{
		Vec4 extent = mMax - mMin;
		return extent.GetX() * extent.GetY() * extent.GetZ() * extent.GetW();
	}

	/// Check if this box contains another box
	bool			Contains(const AABox &inOther) const
	{
		return UVec4::sAnd(Vec4::sLessOrEqual(mMin, inOther.mMin), Vec4::sGreaterOrEqual(mMax, inOther.mMax)).TestAllTrue();
	}

	/// Check if this box contains a point
	bool			Contains(Vec4Arg inOther) const
	{
		return UVec4::sAnd(Vec4::sLessOrEqual(mMin, inOther), Vec4::sGreaterOrEqual(mMax, inOther)).TestAllTrue();
	}

	/// Check if this box contains a point (double precision)
	bool			Contains(DVec4Arg inOther) const
	{
		return Contains(Vec4(float(inOther.GetX()), float(inOther.GetY()), float(inOther.GetZ()), float(inOther.GetW())));
	}

	/// Check if this box overlaps with another box
	bool			Overlaps(const AABox &inOther) const
	{
		return !UVec4::sOr(Vec4::sGreater(mMin, inOther.mMax), Vec4::sLess(mMax, inOther.mMin)).TestAnyTrue();
	}

	/// Check if this box overlaps with a plane
	bool			Overlaps(const Plane &inPlane) const
	{
		Vec4 normal = inPlane.GetNormal();
		float dist_normal = inPlane.SignedDistance(GetSupport(normal));
		float dist_min_normal = inPlane.SignedDistance(GetSupport(-normal));
		return dist_normal * dist_min_normal <= 0.0f;
	}

	/// Translate bounding box
	void			Translate(Vec4Arg inTranslation)
	{
		mMin += inTranslation;
		mMax += inTranslation;
	}

	/// Translate bounding box (double precision)
	void			Translate(DVec4Arg inTranslation)
	{
		mMin = Vec4(float(DVec4(mMin).GetX() + inTranslation.GetX()), float(DVec4(mMin).GetY() + inTranslation.GetY()), float(DVec4(mMin).GetZ() + inTranslation.GetZ()), float(DVec4(mMin).GetW() + inTranslation.GetW()));
		mMax = Vec4(float(DVec4(mMax).GetX() + inTranslation.GetX()), float(DVec4(mMax).GetY() + inTranslation.GetY()), float(DVec4(mMax).GetZ() + inTranslation.GetZ()), float(DVec4(mMax).GetW() + inTranslation.GetW()));
	}

	/// Transform bounding box by rotation matrix and translation
	AABox			Transformed(Mat44Arg inRotation, Vec4Arg inTranslation) const
	{
		// Start with translation
		Vec4 new_min, new_max;
		new_min = new_max = inTranslation;

		// Now find the extreme points by considering the product of the min and max with each column of inRotation
		for (int c = 0; c < 4; ++c)
		{
			Vec4 col = Vec4(inRotation.GetColumn4(c));

			Vec4 a = col * mMin[c];
			Vec4 b = col * mMax[c];

			new_min += Vec4::sMin(a, b);
			new_max += Vec4::sMax(a, b);
		}

		// Return the new bounding box
		return AABox(new_min, new_max);
	}

	/// Transform bounding box by rotation matrix (no translation)
	AABox			Transformed(Mat44Arg inRotation) const
	{
		return Transformed(inRotation, Vec4::sZero());
	}

	/// Transform bounding box by rigid transform (rotation + translation)
	AABox			Transformed(RMat44Arg inTransform) const
	{
		return Transformed(inTransform.GetRotation(), Vec4(inTransform.GetTranslation()));
	}

	/// Scale this bounding box, can handle non-uniform and negative scaling
	AABox			Scaled(Vec4Arg inScale) const
	{
		return AABox::sFromTwoPoints(mMin * inScale, mMax * inScale);
	}

	/// Calculate the support vector for this convex shape.
	Vec4			GetSupport(Vec4Arg inDirection) const
	{
		return Vec4::sSelect(mMax, mMin, Vec4::sLess(inDirection, Vec4::sZero()));
	}

	/// Get the vertices of the face that faces inDirection the most
	/// In 4D, the face of a hyperbox is a 3D cube (8 vertices)
	template <class VERTEX_ARRAY>
	void			GetSupportingFace(Vec4Arg inDirection, VERTEX_ARRAY &outVertices) const
	{
		outVertices.resize(8);

		int axis = inDirection.Abs().GetHighestComponentIndex();

		// The fixed coordinate on the dominant axis.
		//
		// TODO(4D): note the convention, which is NOT the same as GetSupport's. GetSupportingFace
		// returns the face that inDirection points INTO (mMax for a negative component), matching 3D
		// Jolt and ConvexHullShape::GetSupportingFace, which picks the face whose normal is most
		// anti-parallel to inDirection. Callers rely on it: ConvexShape::sCollideConvexVsConvex asks
		// shape 1 for -penetration_axis and shape 2 for +penetration_axis. An earlier pass through
		// this port "corrected" it to agree with GetSupport, which silently selected the far face of
		// each shape and made the contact manifold collapse to a single point.
		float fixed_val = (inDirection[axis] < 0.0f) ? mMax[axis] : mMin[axis];

		// Generate 2^3 = 8 combinations of min/max for the other 3 coordinates
		int other[3];
		int idx = 0;
		for (int i = 0; i < 4; ++i)
			if (i != axis)
				other[idx++] = i;

		for (int mask = 0; mask < 8; ++mask)
		{
			float coords[4];
			coords[axis] = fixed_val;
			coords[other[0]] = (mask & 1) ? mMax[other[0]] : mMin[other[0]];
			coords[other[1]] = (mask & 2) ? mMax[other[1]] : mMin[other[1]];
			coords[other[2]] = (mask & 4) ? mMax[other[2]] : mMin[other[2]];
			outVertices[mask] = Vec4(coords[0], coords[1], coords[2], coords[3]);
		}
	}

	/// Get the closest point on or in this box to inPoint
	Vec4			GetClosestPoint(Vec4Arg inPoint) const
	{
		return Vec4::sMin(Vec4::sMax(inPoint, mMin), mMax);
	}

	/// Get the squared distance between inPoint and this box (will be 0 if inPoint is inside the box)
	inline float	GetSqDistanceTo(Vec4Arg inPoint) const
	{
		return (GetClosestPoint(inPoint) - inPoint).LengthSq();
	}

	/// Bounding box min and max
	Vec4			mMin;
	Vec4			mMax;
};

JPH_NAMESPACE_END
