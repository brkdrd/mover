// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Collision/BackFaceMode.h>
#include <Jolt/Math/RMat44.h> // TODO(4D): RayCastT is instantiated with RMat44, which needs to be complete here

JPH_NAMESPACE_BEGIN

/// Structure that holds a single ray cast in 4D space (origin + direction, no rotation along path)
template <class Vec, class Mat, class RayCastType>
struct RayCastT
{
	JPH_OVERRIDE_NEW_DELETE

	/// Constructors
								RayCastT() = default; // Allow raycast to be created uninitialized
								RayCastT(typename Vec::ArgType inOrigin, Vec4Arg inDirection) : mOrigin(inOrigin), mDirection(inDirection) { }
								RayCastT(const RayCastT<Vec, Mat, RayCastType> &) = default;

	/// Transform this ray using inTransform (rigid transform with translation)
	RayCastType					Transformed(typename Mat::ArgType inTransform) const
	{
		Vec ray_origin = inTransform * mOrigin;
		Vec4 ray_direction(inTransform * (mOrigin + mDirection) - ray_origin);
		return { ray_origin, ray_direction };
	}

	/// Translate ray using inTranslation
	RayCastType					Translated(typename Vec::ArgType inTranslation) const
	{
		return { inTranslation + mOrigin, mDirection };
	}

	/// Get point with fraction inFraction on ray (0 = start of ray, 1 = end of ray)
	inline Vec					GetPointOnRay(float inFraction) const
	{
		return mOrigin + inFraction * mDirection;
	}

	Vec							mOrigin;					///< Origin of the ray (4D position)
	Vec4						mDirection;					///< Direction and length of the ray (4D). Anything beyond this length will not be reported as a hit.
};

struct RayCast : public RayCastT<Vec4, RMat44, RayCast>
{
	using RayCastT<Vec4, RMat44, RayCast>::RayCastT;
};

struct RRayCast : public RayCastT<RVec4, RMat44, RRayCast>
{
	using RayCastT<RVec4, RMat44, RRayCast>::RayCastT;

	/// Convert from RayCast, converts single to double precision (no-op in single-precision mode)
	explicit					RRayCast(const RayCast &inRay) :
		RRayCast(RVec4(inRay.mOrigin), inRay.mDirection)
	{
	}

	/// Convert to RayCast, which implies casting from double precision to single precision
	explicit					operator RayCast() const
	{
		return RayCast(Vec4(mOrigin), mDirection);
	}
};

/// Settings to be passed with a ray cast
class RayCastSettings
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Set the backfacing mode for all shapes
	void						SetBackFaceMode(EBackFaceMode inMode) { mBackFaceModeTriangles = mBackFaceModeConvex = inMode; }

	/// How backfacing triangles should be treated (should we report back facing hits for triangle based shapes, e.g. MeshShape/HeightFieldShape?)
	EBackFaceMode				mBackFaceModeTriangles		= EBackFaceMode::IgnoreBackFaces;

	/// How backfacing convex objects should be treated (should we report back facing hits for convex shapes?)
	EBackFaceMode				mBackFaceModeConvex			= EBackFaceMode::IgnoreBackFaces;

	/// If convex shapes should be treated as solid. When true, a ray starting inside a convex shape will generate a hit at fraction 0.
	bool						mTreatConvexAsSolid			= true;
};

JPH_NAMESPACE_END
