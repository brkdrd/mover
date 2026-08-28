// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

/// A simple triangle and its material
class Triangle
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructor
					Triangle() = default;
					Triangle(const Float4 &inV1, const Float4 &inV2, const Float4 &inV3, uint32 inMaterialIndex = 0, uint32 inUserData = 0) : mV { inV1, inV2, inV3 }, mMaterialIndex(inMaterialIndex), mUserData(inUserData) { }
					Triangle(Vec4Arg inV1, Vec4Arg inV2, Vec4Arg inV3, uint32 inMaterialIndex = 0, uint32 inUserData = 0) : mMaterialIndex(inMaterialIndex), mUserData(inUserData) { inV1.StoreFloat4(&mV[0]); inV2.StoreFloat4(&mV[1]); inV3.StoreFloat4(&mV[2]); }

	/// Get center of triangle
	Vec4			GetCentroid() const
	{
		return (Vec4(mV[0]) + Vec4(mV[1]) + Vec4(mV[2])) * (1.0f / 3.0f);
	}

	/// Vertices
	Float4			mV[3];
	uint32			mMaterialIndex = 0;
	uint32			mUserData = 0;				///< User data that can be used for anything by the application, e.g. for tracking the original index of the triangle
};

using TriangleList = Array<Triangle>;

JPH_NAMESPACE_END
