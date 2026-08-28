// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Math/Vec4.h>
#include <Jolt/Math/Bivec.h>
#include <Jolt/Math/Trigonometry.h>
#include <Jolt/Core/HashCombine.h>

// Create a std::hash/JPH::Hash for Rotor
JPH_MAKE_HASHABLE(JPH::Rotor, t.GetScalar(), t.GetE12(), t.GetE13(), t.GetE14(), t.GetE23(), t.GetE24(), t.GetE34(), t.GetPseudoscalar())

JPH_NAMESPACE_BEGIN

/// Spin(4) = SU(2)xSU(2). The pseudoscalar I = e1234 is central with I^2 = +1, so the projectors
/// P+- = (1 +- I)/2 split the even subalgebra into two commuting copies of the quaternions H. A
/// rotor R = [s, e12, e13, e14, e23, e24, e34, e1234] maps to two quaternions
///   qL = (s + p, e12 - e34, e13 + e24, e14 - e23)   (the +I eigenspace)
///   qR = (s - p, e12 + e34, e13 - e24, e14 + e23)   (the -I eigenspace)
/// and R is reconstructed from the pair. Normalization, exp/log and slerp are all just the
/// corresponding per-quaternion operation applied to qL and qR independently. Each quaternion is
/// stored as (w, x, y, z) in a Lane4.
namespace RotorImpl {

JPH_INLINE void sToQuats(const Vec8 &inR, Lane4 &outL, Lane4 &outR)
{
	float s = inR[0], b12 = inR[1], b13 = inR[2], b14 = inR[3], b23 = inR[4], b24 = inR[5], b34 = inR[6], p = inR[7];
	outL = Lane4(s + p, b12 - b34, b13 + b24, b14 - b23);
	outR = Lane4(s - p, b12 + b34, b13 - b24, b14 + b23);
}

JPH_INLINE Rotor sFromQuats(Lane4Arg inL, Lane4Arg inR)
{
	float l0 = inL.GetX(), l1 = inL.GetY(), l2 = inL.GetZ(), l3 = inL.GetW();
	float r0 = inR.GetX(), r1 = inR.GetY(), r2 = inR.GetZ(), r3 = inR.GetW();
	return Rotor(
		0.5f * (l0 + r0),   // s
		0.5f * (l1 + r1),   // e12
		0.5f * (l2 + r2),   // e13
		0.5f * (l3 + r3),   // e14
		0.5f * (r3 - l3),   // e23
		0.5f * (l2 - r2),   // e24
		0.5f * (r1 - l1),   // e34
		0.5f * (l0 - r0));  // e1234
}

/// Exponential of a pure-imaginary quaternion whose (x, y, z) sit in lanes (Y, Z, W):
/// exp(0, v) = (cos|v|, sin|v|/|v| * v)
JPH_INLINE Lane4 sQuatExp(Lane4Arg inImag)
{
	float x = inImag.GetY(), y = inImag.GetZ(), z = inImag.GetW();
	float angle = sqrt(x * x + y * y + z * z);
	float c = Cos(angle);
	float k = angle > 1.0e-9f? Sin(angle) / angle : 1.0f; // -> 1 as angle -> 0
	return Lane4(c, k * x, k * y, k * z);
}

/// Logarithm of a unit quaternion q = (w, x, y, z): the imaginary part (0, angle/|v| * v)
JPH_INLINE Lane4 sQuatLog(Lane4Arg inQuat)
{
	float x = inQuat.GetY(), y = inQuat.GetZ(), z = inQuat.GetW();
	float v_len = sqrt(x * x + y * y + z * z);
	float angle = atan2(v_len, inQuat.GetX());
	float k = v_len > 1.0e-9f? angle / v_len : 1.0f;
	return Lane4(0.0f, k * x, k * y, k * z);
}

} // namespace RotorImpl

Rotor Rotor::sRotation(Vec4Arg inPlaneA, Vec4Arg inPlaneB, float inAngle)
{
	JPH_ASSERT(inPlaneA.IsNormalized());
	JPH_ASSERT(inPlaneB.IsNormalized());
	JPH_ASSERT(abs(inPlaneA.Dot(inPlaneB)) < 1.0e-4f);

	float half = 0.5f * inAngle;
	float c = Cos(half);
	float s = Sin(half);

	// Wedge product components: (a ^ b)_ij = a_i*b_j - a_j*b_i
	float a1 = inPlaneA.GetX(), a2 = inPlaneA.GetY(), a3 = inPlaneA.GetZ(), a4 = inPlaneA.GetW();
	float b1 = inPlaneB.GetX(), b2 = inPlaneB.GetY(), b3 = inPlaneB.GetZ(), b4 = inPlaneB.GetW();

	return Rotor(c,
				 s * (a1 * b2 - a2 * b1),  // e12
				 s * (a1 * b3 - a3 * b1),  // e13
				 s * (a1 * b4 - a4 * b1),  // e14
				 s * (a2 * b3 - a3 * b2),  // e23
				 s * (a2 * b4 - a4 * b2),  // e24
				 s * (a3 * b4 - a4 * b3),  // e34
				 0.0f);                     // pseudoscalar = 0 for simple rotations
}

Rotor Rotor::sDoubleRotation(Vec4Arg inPlaneA1, Vec4Arg inPlaneB1, float inAngle1,
							 Vec4Arg inPlaneA2, Vec4Arg inPlaneB2, float inAngle2)
{
	// Compose two simple rotations in orthogonal planes
	Rotor r1 = sRotation(inPlaneA1, inPlaneB1, inAngle1);
	Rotor r2 = sRotation(inPlaneA2, inPlaneB2, inAngle2);
	return r1 * r2;
}

Rotor Rotor::operator * (RotorArg inRHS) const
{
	// Geometric product in the even subalgebra of Cl(4,0)
	// Layout: [s, e12, e13, e14, e23, e24, e34, e1234]
	float a0 = mValue[0], a1 = mValue[1], a2 = mValue[2], a3 = mValue[3];
	float a4 = mValue[4], a5 = mValue[5], a6 = mValue[6], a7 = mValue[7];
	float b0 = inRHS.mValue[0], b1 = inRHS.mValue[1], b2 = inRHS.mValue[2], b3 = inRHS.mValue[3];
	float b4 = inRHS.mValue[4], b5 = inRHS.mValue[5], b6 = inRHS.mValue[6], b7 = inRHS.mValue[7];

	return Rotor(
		a0*b0 - a1*b1 - a2*b2 - a3*b3 - a4*b4 - a5*b5 - a6*b6 + a7*b7,  // s
		a0*b1 + a1*b0 - a2*b4 - a3*b5 + a4*b2 + a5*b3 - a6*b7 - a7*b6,  // e12
		a0*b2 + a1*b4 + a2*b0 - a3*b6 - a4*b1 + a5*b7 + a6*b3 + a7*b5,  // e13
		a0*b3 + a1*b5 + a2*b6 + a3*b0 - a4*b7 - a5*b1 - a6*b2 - a7*b4,  // e14
		a0*b4 - a1*b2 + a2*b1 - a3*b7 + a4*b0 - a5*b6 + a6*b5 - a7*b3,  // e23
		a0*b5 - a1*b3 + a2*b7 + a3*b1 + a4*b6 + a5*b0 - a6*b4 + a7*b2,  // e24
		a0*b6 - a1*b7 - a2*b3 + a3*b2 - a4*b5 + a5*b4 + a6*b0 - a7*b1,  // e34
		a0*b7 + a1*b6 - a2*b5 + a3*b4 + a4*b3 - a5*b2 + a6*b1 + a7*b0   // e1234
	);
}

Vec4 Rotor::operator * (Vec4Arg inValue) const
{
	// Sandwich product: R * v * ~R
	// Computed in two steps through the full Clifford algebra.
	JPH_ASSERT(IsNormalized());

	float s = mValue[0], a = mValue[1], b = mValue[2], c = mValue[3];
	float d = mValue[4], e = mValue[5], f = mValue[6], p = mValue[7];
	float v1 = inValue.GetX(), v2 = inValue.GetY(), v3 = inValue.GetZ(), v4 = inValue.GetW();

	// Step 1: T = R * v (odd multivector: grade-1 and grade-3 parts)
	float t1   = s*v1 + a*v2 + b*v3 + c*v4;
	float t2   = s*v2 - a*v1 + d*v3 + e*v4;
	float t3   = s*v3 - b*v1 - d*v2 + f*v4;
	float t4   = s*v4 - c*v1 - e*v2 - f*v3;
	float t123 = d*v1 - b*v2 + a*v3 + p*v4;
	float t124 = e*v1 - c*v2 + a*v4 - p*v3;
	float t134 = f*v1 + p*v2 - c*v3 + b*v4;
	float t234 = -p*v1 + f*v2 - e*v3 + d*v4;

	// Step 2: Result = (T * ~R) extracting grade-1 components
	// ~R = [s, -a, -b, -c, -d, -e, -f, p]
	float r1 = t1*s + t2*a + t3*b + t4*c + t123*d + t124*e + t134*f + t234*p;
	float r2 = t2*s - t1*a + t3*d + t4*e - t123*b - t124*c + t234*f - t134*p;
	float r3 = t3*s - t1*b - t2*d + t4*f + t123*a - t134*c - t234*e + t124*p;
	float r4 = t4*s - t1*c - t2*e - t3*f + t124*a + t134*b + t234*d - t123*p;

	return Vec4(r1, r2, r3, r4);
}

Bivec Rotor::RotateBivec(BivecArg inBivec) const
{
	// Ad(R)(e_i ^ e_j) = (R e_i ~R) ^ (R e_j ~R). The images of the basis vectors are the columns of
	// the rotation, so the adjoint action on a bivector is linear in the wedges of those images.
	Vec4 f1 = *this * Vec4::sAxisX();
	Vec4 f2 = *this * Vec4::sAxisY();
	Vec4 f3 = *this * Vec4::sAxisZ();
	Vec4 f4 = *this * Vec4::sAxisW();

	return inBivec.GetE12() * Bivec::sWedge(f1, f2)
		 + inBivec.GetE13() * Bivec::sWedge(f1, f3)
		 + inBivec.GetE14() * Bivec::sWedge(f1, f4)
		 + inBivec.GetE23() * Bivec::sWedge(f2, f3)
		 + inBivec.GetE24() * Bivec::sWedge(f2, f4)
		 + inBivec.GetE34() * Bivec::sWedge(f3, f4);
}

Bivec Rotor::InverseRotateBivec(BivecArg inBivec) const
{
	// Ad(~R)(B): use the reversed (inverse) rotor's basis-vector images.
	Rotor r = Reversed();
	Vec4 f1 = r * Vec4::sAxisX();
	Vec4 f2 = r * Vec4::sAxisY();
	Vec4 f3 = r * Vec4::sAxisZ();
	Vec4 f4 = r * Vec4::sAxisW();

	return inBivec.GetE12() * Bivec::sWedge(f1, f2)
		 + inBivec.GetE13() * Bivec::sWedge(f1, f3)
		 + inBivec.GetE14() * Bivec::sWedge(f1, f4)
		 + inBivec.GetE23() * Bivec::sWedge(f2, f3)
		 + inBivec.GetE24() * Bivec::sWedge(f2, f4)
		 + inBivec.GetE34() * Bivec::sWedge(f3, f4);
}

Rotor Rotor::Reversed() const
{
	// Negate grade-2 (bivectors), keep grade-0 (scalar) and grade-4 (pseudoscalar)
	// ~R = [s, -e12, -e13, -e14, -e23, -e24, -e34, e1234]
	return Rotor(mValue[0],
				 -mValue[1], -mValue[2], -mValue[3],
				 -mValue[4], -mValue[5], -mValue[6],
				 mValue[7]);
}

Rotor Rotor::Inversed() const
{
	// ~R / |R|^2
	float len_sq = LengthSq();
	JPH_ASSERT(len_sq > 0.0f);
	float inv_len_sq = 1.0f / len_sq;
	Rotor rev = Reversed();
	return Rotor(rev.mValue * inv_len_sq);
}

bool Rotor::IsNormalized(float inTolerance) const
{
	Lane4 ql, qr;
	RotorImpl::sToQuats(mValue, ql, qr);
	// Both SU(2) factors must be unit length (this also implies the 8-vector has length 1)
	return abs(ql.LengthSq() - 1.0f) <= 2.0f * inTolerance
		&& abs(qr.LengthSq() - 1.0f) <= 2.0f * inTolerance;
}

Rotor Rotor::Normalized() const
{
	Lane4 ql, qr;
	RotorImpl::sToQuats(mValue, ql, qr);
	float ll = ql.Length(), lr = qr.Length();
	JPH_ASSERT(ll > 0.0f && lr > 0.0f);
	return RotorImpl::sFromQuats(ql * (1.0f / ll), qr * (1.0f / lr));
}

Rotor Rotor::sExp(BivecArg inBivector)
{
	// Split the bivector into its two SU(2) parts and exponentiate each like a quaternion.
	// Build the bivector in Rotor lane layout [s, e12, e13, e14, e23, e24, e34, e1234] with s=p=0.
	Vec8 b(0.0f, inBivector.GetE12(), inBivector.GetE13(), inBivector.GetE14(),
		   inBivector.GetE23(), inBivector.GetE24(), inBivector.GetE34(), 0.0f);
	Lane4 vl, vr;
	RotorImpl::sToQuats(b, vl, vr); // pure-imaginary: the real (X) lane is 0 for a bivector
	return RotorImpl::sFromQuats(RotorImpl::sQuatExp(vl), RotorImpl::sQuatExp(vr));
}

Bivec Rotor::Log() const
{
	JPH_ASSERT(IsNormalized());
	Lane4 ql, qr;
	RotorImpl::sToQuats(mValue, ql, qr);
	// Log of each factor is pure imaginary; reconstruct gives a rotor with s = p = 0, i.e. a bivector
	Rotor b = RotorImpl::sFromQuats(RotorImpl::sQuatLog(ql), RotorImpl::sQuatLog(qr));
	return Bivec(b.GetE12(), b.GetE13(), b.GetE14(), b.GetE23(), b.GetE24(), b.GetE34());
}

Rotor Rotor::LERP(RotorArg inDestination, float inFraction) const
{
	// Blend then project back onto Spin(4) so the result is a valid rotor
	float scale0 = 1.0f - inFraction;
	return Rotor(mValue * scale0 + inDestination.mValue * inFraction).Normalized();
}

Rotor Rotor::SLERP(RotorArg inDestination, float inFraction) const
{
	// Difference at which to LERP instead of SLERP
	const float delta = 0.0001f;

	// Calc cosine
	float sign_scale1 = 1.0f;
	float cos_omega = Dot(inDestination);

	// Adjust signs — R and -R represent the same rotation
	if (cos_omega < 0.0f)
	{
		cos_omega = -cos_omega;
		sign_scale1 = -1.0f;
	}

	// Calculate coefficients
	float scale0, scale1;
	if (1.0f - cos_omega > delta)
	{
		// Standard case (slerp)
		float omega = ACos(cos_omega);
		float sin_omega = Sin(omega);
		scale0 = Sin((1.0f - inFraction) * omega) / sin_omega;
		scale1 = sign_scale1 * Sin(inFraction * omega) / sin_omega;
	}
	else
	{
		// Rotors are very close so we can do a linear interpolation
		scale0 = 1.0f - inFraction;
		scale1 = sign_scale1 * inFraction;
	}

	return Rotor(mValue * scale0 + inDestination.mValue * scale1).Normalized();
}

JPH_NAMESPACE_END
