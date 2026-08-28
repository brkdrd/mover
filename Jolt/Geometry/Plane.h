// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

/// An infinite hyperplane described by the formula X . Normal + Constant = 0.
/// In 4D, a hyperplane has a 4-component normal vector and a scalar constant.
class [[nodiscard]] Plane
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructor
					Plane() = default;
					Plane(Vec4Arg inNormal, float inConstant)								: mNormal(inNormal), mConstant(inConstant) { }

	/// Create from point and normal
	static Plane	sFromPointAndNormal(Vec4Arg inPoint, Vec4Arg inNormal)					{ return Plane(inNormal, -inNormal.Dot(inPoint)); }

	/// Create from point and normal, double precision version that more accurately calculates the plane constant
	static Plane	sFromPointAndNormal(DVec4Arg inPoint, Vec4Arg inNormal)					{ return Plane(inNormal, -float(DVec4(inNormal).Dot(inPoint))); }

	/// Create from 4 counter clockwise points (4 points needed to define a hyperplane in 4D)
	static Plane	sFromPointsCCW(Vec4Arg inV1, Vec4Arg inV2, Vec4Arg inV3, Vec4Arg inV4)	{ return sFromPointAndNormal(inV1, Vec4::sCross(inV2 - inV1, inV3 - inV1, inV4 - inV1).Normalized()); }

	// Properties
	Vec4			GetNormal() const														{ return mNormal; }
	void			SetNormal(Vec4Arg inNormal)												{ mNormal = inNormal; }
	float			GetConstant() const														{ return mConstant; }
	void			SetConstant(float inConstant)											{ mConstant = inConstant; }

	/// Store normal as 4 floats + constant
	void			StoreNormal(Float4 *outNormal) const									{ mNormal.StoreFloat4(outNormal); }

	/// Offset the plane (positive value means move it in the direction of the plane normal)
	Plane			Offset(float inDistance) const											{ return Plane(mNormal, mConstant - inDistance); }

	/// Scale the plane by a (possibly non-uniform) scale applied to the space it lives in.
	/// A point x on the plane maps to inScale * x; the transformed plane normal is (mNormal / inScale), renormalized.
	Plane			Scaled(Vec4Arg inScale) const
	{
		Vec4 scaled_normal = mNormal / inScale;
		float scaled_normal_length = scaled_normal.Length();
		return Plane(scaled_normal / scaled_normal_length, mConstant / scaled_normal_length);
	}

	/// Transform the plane by a rotation matrix and translation
	inline Plane	GetTransformed(Mat44Arg inRotation, Vec4Arg inTranslation) const
	{
		Vec4 transformed_normal = Vec4(inRotation * Lane4(mNormal));
		return Plane(transformed_normal, mConstant - inTranslation.Dot(transformed_normal));
	}

	/// Distance point to plane
	float			SignedDistance(Vec4Arg inPoint) const									{ return inPoint.Dot(mNormal) + mConstant; }

	/// Project inPoint onto the plane
	Vec4			ProjectPointOnPlane(Vec4Arg inPoint) const								{ return inPoint - mNormal * SignedDistance(inPoint); }

	/// Returns intersection point of 4 hyperplanes (4 needed to get a point in 4D)
	/// Solves the 4x4 system: [n1; n2; n3; n4] * x = [-c1; -c2; -c3; -c4]
	static bool		sIntersectPlanes(const Plane &inP1, const Plane &inP2, const Plane &inP3, const Plane &inP4, Vec4 &outPoint)
	{
		Vec4 n1 = inP1.mNormal;
		Vec4 n2 = inP2.mNormal;
		Vec4 n3 = inP3.mNormal;
		Vec4 n4 = inP4.mNormal;

		// Vec4::sCross(a, b, c) is defined so that u . sCross(a, b, c) == det[u; a; b; c]. The row
		// cofactor vector v_i of M = [n1; n2; n3; n4] is the vector with u . v_i == det(M with row i
		// replaced by u), so v_i is sCross of the other three normals up to the sign of moving u from
		// row i to row 1: v1 = +sCross(n2,n3,n4), v2 = -sCross(n1,n3,n4), v3 = +sCross(n1,n2,n4),
		// v4 = -sCross(n1,n2,n3). The alternating signs are applied in the sum below.
		Vec4 v1 = Vec4::sCross(n2, n3, n4);
		Vec4 v2 = Vec4::sCross(n1, n3, n4);
		Vec4 v3 = Vec4::sCross(n1, n2, n4);
		Vec4 v4 = Vec4::sCross(n1, n2, n3);

		// det(M) = det[n1; n2; n3; n4] = n1 . v1
		float det = n1.Dot(v1);
		if (abs(det) < FLT_EPSILON)
			return false;

		// x = M^-1 rhs = adj(M) rhs / det, and column k of adj(M) is the row cofactor vector v_k,
		// so x = sum_i (-1)^(i+1) rhs_i v_i / det
		Vec4 rhs(-inP1.mConstant, -inP2.mConstant, -inP3.mConstant, -inP4.mConstant);
		outPoint = (rhs.GetX() * v1 - rhs.GetY() * v2 + rhs.GetZ() * v3 - rhs.GetW() * v4) * (1.0f / det);
		return true;
	}

private:
#ifdef JPH_OBJECT_STREAM
	friend void		CreateRTTIPlane(class RTTI &);										// For JPH_IMPLEMENT_SERIALIZABLE_OUTSIDE_CLASS
#endif

	Vec4			mNormal;														///< 4D normal vector of the hyperplane
	float			mConstant;														///< Plane constant: x . normal + constant = 0
};

JPH_NAMESPACE_END
