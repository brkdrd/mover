// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Geometry/Triangle.h>
#include <Jolt/Geometry/IndexedTriangle.h>
#include <Jolt/Geometry/AABox.h>
#include <Jolt/Math/Mat44.h>

JPH_NAMESPACE_BEGIN

class AABox;

/// Oriented box (4D)
class JPH_EXPORT_GCC_BUG_WORKAROUND [[nodiscard]] OrientedBox
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructor
					OrientedBox() = default;
					OrientedBox(Mat44Arg inOrientation, Vec4Arg inTranslation, Vec4Arg inHalfExtents) : mOrientation(inOrientation), mTranslation(inTranslation), mHalfExtents(inHalfExtents) { }

	/// Construct from axis aligned box, rotation and translation. Only works for rotation matrix (no scaling / shearing).
					OrientedBox(Mat44Arg inOrientation, Vec4Arg inTranslation, const AABox &inBox) : mOrientation(inOrientation), mTranslation(inTranslation + Vec4(inOrientation * Lane4(inBox.GetCenter()))), mHalfExtents(inBox.GetExtent()) { }

	/// Test if oriented box overlaps with axis aligned box
	bool			Overlaps(const AABox &inBox, float inEpsilon = 1.0e-6f) const;

	/// Test if two oriented boxes overlap each other
	bool			Overlaps(const OrientedBox &inBox, float inEpsilon = 1.0e-6f) const;

	Mat44			mOrientation;														///< Pure rotation matrix for the OBB
	Vec4			mTranslation;														///< Center position of the OBB
	Vec4			mHalfExtents;														///< Half extents (half the size of each edge) in local space
};

JPH_NAMESPACE_END
