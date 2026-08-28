// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Mat44.h>

JPH_NAMESPACE_BEGIN

/// Helper functions to get the support point for a convex object
/// Structure that transforms a convex object (rotation + translation, no scaling)
template <typename ConvexObject>
struct TransformedConvexObject
{
	/// Create transformed convex object.
	TransformedConvexObject(Mat44Arg inRotation, Vec4Arg inTranslation, const ConvexObject &inObject) :
		mRotation(inRotation),
		mTranslation(inTranslation),
		mObject(inObject)
	{
	}

	/// Calculate the support vector for this convex shape.
	Vec4					GetSupport(Vec4Arg inDirection) const
	{
		return Vec4(mRotation * Lane4(mObject.GetSupport(Vec4(mRotation.Transposed() * Lane4(inDirection))))) + mTranslation;
	}

	/// Get the vertices of the face that faces inDirection the most
	template <class VERTEX_ARRAY>
	void					GetSupportingFace(Vec4Arg inDirection, VERTEX_ARRAY &outVertices) const
	{
		mObject.GetSupportingFace(Vec4(mRotation.Transposed() * Lane4(inDirection)), outVertices);

		for (Vec4 &v : outVertices)
			v = Vec4(mRotation * Lane4(v)) + mTranslation;
	}

	Mat44					mRotation;
	Vec4					mTranslation;
	const ConvexObject &	mObject;
};

/// Structure that adds a convex radius
template <typename ConvexObject>
struct AddConvexRadius
{
	AddConvexRadius(const ConvexObject &inObject, float inRadius) :
		mObject(inObject),
		mRadius(inRadius)
	{
	}

	/// Calculate the support vector for this convex shape.
	Vec4					GetSupport(Vec4Arg inDirection) const
	{
		float length = inDirection.Length();
		return length > 0.0f ? mObject.GetSupport(inDirection) + (mRadius / length) * inDirection : mObject.GetSupport(inDirection);
	}

	const ConvexObject &	mObject;
	float					mRadius;
};

/// Structure that performs a Minkowski difference A - B
template <typename ConvexObjectA, typename ConvexObjectB>
struct MinkowskiDifference
{
	MinkowskiDifference(const ConvexObjectA &inObjectA, const ConvexObjectB &inObjectB) :
		mObjectA(inObjectA),
		mObjectB(inObjectB)
	{
	}

	/// Calculate the support vector for this convex shape.
	Vec4					GetSupport(Vec4Arg inDirection) const
	{
		return mObjectA.GetSupport(inDirection) - mObjectB.GetSupport(-inDirection);
	}

	const ConvexObjectA &	mObjectA;
	const ConvexObjectB &	mObjectB;
};

/// Class that wraps a point so that it can be used with convex collision detection
struct PointConvexSupport
{
	/// Calculate the support vector for this convex shape.
	Vec4					GetSupport([[maybe_unused]] Vec4Arg inDirection) const
	{
		return mPoint;
	}

	Vec4					mPoint;
};

/// Class that wraps a triangle so that it can used with convex collision detection
struct TriangleConvexSupport
{
	/// Constructor
							TriangleConvexSupport(Vec4Arg inV1, Vec4Arg inV2, Vec4Arg inV3) :
		mV1(inV1),
		mV2(inV2),
		mV3(inV3)
	{
	}

	/// Calculate the support vector for this convex shape.
	Vec4					GetSupport(Vec4Arg inDirection) const
	{
		// Project vertices on inDirection
		float d1 = mV1.Dot(inDirection);
		float d2 = mV2.Dot(inDirection);
		float d3 = mV3.Dot(inDirection);

		// Return vertex with biggest projection
		if (d1 > d2)
		{
			if (d1 > d3)
				return mV1;
			else
				return mV3;
		}
		else
		{
			if (d2 > d3)
				return mV2;
			else
				return mV3;
		}
	}

	/// Get the vertices of the face that faces inDirection the most
	template <class VERTEX_ARRAY>
	void					GetSupportingFace([[maybe_unused]] Vec4Arg inDirection, VERTEX_ARRAY &outVertices) const
	{
		outVertices.push_back(mV1);
		outVertices.push_back(mV2);
		outVertices.push_back(mV3);
	}

	/// The three vertices of the triangle
	Vec4					mV1;
	Vec4					mV2;
	Vec4					mV3;
};

/// Class that wraps a tetrahedron so that it can be used with convex collision detection
struct TetrahedronConvexSupport
{
	/// Constructor
							TetrahedronConvexSupport(Vec4Arg inV1, Vec4Arg inV2, Vec4Arg inV3, Vec4Arg inV4) :
		mV1(inV1),
		mV2(inV2),
		mV3(inV3),
		mV4(inV4)
	{
	}

	/// Calculate the support vector for this convex shape.
	Vec4					GetSupport(Vec4Arg inDirection) const
	{
		float d1 = mV1.Dot(inDirection);
		float d2 = mV2.Dot(inDirection);
		float d3 = mV3.Dot(inDirection);
		float d4 = mV4.Dot(inDirection);

		float max_d = d1;
		Vec4 result = mV1;
		if (d2 > max_d) { max_d = d2; result = mV2; }
		if (d3 > max_d) { max_d = d3; result = mV3; }
		if (d4 > max_d) { result = mV4; }
		return result;
	}

	Vec4					mV1;
	Vec4					mV2;
	Vec4					mV3;
	Vec4					mV4;
};

/// Class that wraps a polygon so that it can used with convex collision detection
template <class VERTEX_ARRAY>
struct PolygonConvexSupport
{
	/// Constructor
	explicit				PolygonConvexSupport(const VERTEX_ARRAY &inVertices) :
		mVertices(inVertices)
	{
	}

	/// Calculate the support vector for this convex shape.
	Vec4					GetSupport(Vec4Arg inDirection) const
	{
		Vec4 support_point = mVertices[0];
		float best_dot = mVertices[0].Dot(inDirection);

		for (typename VERTEX_ARRAY::const_iterator v = mVertices.begin() + 1; v < mVertices.end(); ++v)
		{
			float dot = v->Dot(inDirection);
			if (dot > best_dot)
			{
				best_dot = dot;
				support_point = *v;
			}
		}

		return support_point;
	}

	/// Get the vertices of the face that faces inDirection the most
	template <class VERTEX_ARRAY_ARG>
	void					GetSupportingFace([[maybe_unused]] Vec4Arg inDirection, VERTEX_ARRAY_ARG &outVertices) const
	{
		for (Vec4 v : mVertices)
			outVertices.push_back(v);
	}

	/// The vertices of the polygon
	const VERTEX_ARRAY &	mVertices;
};

JPH_NAMESPACE_END
