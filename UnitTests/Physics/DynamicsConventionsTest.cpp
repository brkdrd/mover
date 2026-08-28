// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Math/Bivec.h>
#include <Jolt/Math/Rotor.h>

JPH_SUPPRESS_WARNINGS_STD_BEGIN
#include <random>
JPH_SUPPRESS_WARNINGS_STD_END

// TODO(4D): Pins the sign conventions that tie together the three halves of 4D rigid body
// dynamics: the rotor integrator (Body::AddRotationStep), the point velocity operator
// (MotionProperties::GetPointVelocityCOM -> Bivec::Contract) and the torque / angular impulse
// operator (Body::AddForce/AddImpulse -> Bivec::sWedge). In 3D all three use the cross product
// and the conventions agree trivially; in 4D they do NOT, because the left contraction is the
// *negative* adjoint of the wedge:
//
//     <Omega _| r, F>_R4 = -<Omega, r ^ F>_bivec = +<Omega, F ^ r>_bivec
//
// so the correct 4D analog of the 3D expression "a x b" appearing in a torque or Jacobian role
// is "b ^ a", with the operands swapped. These tests fix that choice against first principles
// (finite differencing the actual integrator, and conservation of energy) rather than against
// an analogy with the 3D code.
//
// Deliberately depends only on Jolt/Math headers so it can run in the standalone Math test
// runner while the physics library does not link yet (CLAUDE.md, Tests section).
TEST_SUITE("DynamicsConventionsTests")
{
	static Vec4 sRandomVec(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		return Vec4(dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom));
	}

	static Bivec sRandomBivec(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> dist(-1.0f, 1.0f);
		return Bivec(dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom));
	}

	/// A generic (non simple) rotor: a product of simple rotations in three different planes
	static Rotor sRandomRotor(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> angle(-JPH_PI, JPH_PI);
		return Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), angle(ioRandom))
			 * Rotor::sRotation(Vec4::sAxisZ(), Vec4::sAxisW(), angle(ioRandom))
			 * Rotor::sRotation(Vec4::sAxisY(), Vec4::sAxisW(), angle(ioRandom));
	}

	/// The velocity of a body fixed point must be the time derivative of that point under the
	/// rotation integrator that Body::AddRotationStep actually applies:
	///     R(t) = exp(Omega t / 2) R(0),   r(t) = R(t) r_body ~R(t)
	/// Finite difference r(t) and compare against Omega.Contract(r). This is what forces the sign
	/// of every torque in the engine: once the point velocity operator is pinned, the generalized
	/// force conjugate to Omega is pinned with it.
	TEST_CASE("TestPointVelocityIsDerivativeOfIntegratedRotation")
	{
		UnitTestRandom random(0x4D1234);
		for (int n = 0; n < 25; ++n)
		{
			Rotor r0 = sRandomRotor(random);
			Bivec omega = sRandomBivec(random);
			Vec4 r_body = sRandomVec(random);

			// Position of the body fixed point at t = 0, in world space relative to the COM
			Vec4 r = r0 * r_body;

			// Central difference of the integrated rotation
			const float h = 1.0e-3f;
			Vec4 r_plus = (Rotor::sExp(0.5f * (h * omega)) * r0) * r_body;
			Vec4 r_minus = (Rotor::sExp(0.5f * (-h * omega)) * r0) * r_body;
			Vec4 numerical = (r_plus - r_minus) / (2.0f * h);

			CHECK((numerical - omega.Contract(r)).Length() < 1.0e-3f);
		}
	}

	/// Power / virtual work identity. A force F applied at world space offset r from the center of
	/// mass delivers power F . v_point. Splitting v_point = v + Omega _| r shows that the bivector
	/// conjugate to Omega (i.e. the torque) is F ^ r, NOT r ^ F.
	TEST_CASE("TestTorqueIsForceWedgeOffset")
	{
		UnitTestRandom random(0x4D5678);
		for (int n = 0; n < 50; ++n)
		{
			Vec4 r = sRandomVec(random);
			Vec4 f = sRandomVec(random);
			Bivec omega = sRandomBivec(random);

			// Power delivered into the rotational degrees of freedom
			float power = f.Dot(omega.Contract(r));

			// It equals <torque, Omega> for torque = F ^ r
			CHECK_APPROX_EQUAL(power, Bivec::sWedge(f, r).Dot(omega), 1.0e-5f);

			// and is the negative of the naive 3D transliteration r ^ F
			CHECK_APPROX_EQUAL(power, -Bivec::sWedge(r, f).Dot(omega), 1.0e-5f);
		}
	}

	/// Conservation of energy for an impulse applied at a point of a body at rest.
	///
	/// For a rigid body starting at rest, the work done by an impulse P applied at offset r is
	/// exactly half the dot product of P with the resulting velocity of that material point (the
	/// average of the velocity before and after). That must equal the resulting kinetic energy:
	///
	///     1/2 m |v|^2 + 1/2 <Omega, I Omega>  ==  1/2 P . v_point_after
	///
	/// This is the decisive test for the torque sign: the rotational term enters the right hand
	/// side linearly through Contract but the left hand side quadratically, so flipping the sign of
	/// the angular impulse breaks the balance by exactly twice the rotational energy.
	TEST_CASE("TestImpulseAtPointConservesEnergy")
	{
		UnitTestRandom random(0x4D9ABC);
		uniform_real_distribution<float> positive(0.5f, 4.0f);

		for (int n = 0; n < 50; ++n)
		{
			float mass = positive(random);

			// Diagonal (principal frame) bivector inertia, body aligned with the world for this test
			Bivec inertia(positive(random), positive(random), positive(random), positive(random), positive(random), positive(random));
			Bivec inv_inertia(1.0f / inertia.GetE12(), 1.0f / inertia.GetE13(), 1.0f / inertia.GetE14(),
							  1.0f / inertia.GetE23(), 1.0f / inertia.GetE24(), 1.0f / inertia.GetE34());

			Vec4 r = sRandomVec(random);
			Vec4 p = sRandomVec(random);

			// Apply the impulse: linear P / m, angular I^-1 (P ^ r)
			Vec4 v = p / mass;
			Bivec omega = inv_inertia * Bivec::sWedge(p, r);

			// Resulting kinetic energy
			float energy = 0.5f * mass * v.LengthSq() + 0.5f * omega.Dot(inertia * omega);

			// Work done by the impulse = P . (average velocity of the material point)
			Vec4 v_point = v + omega.Contract(r);
			float work = 0.5f * p.Dot(v_point);

			CHECK_APPROX_EQUAL(energy, work, 1.0e-4f * max(1.0f, abs(energy)));

			// An impulse always does non negative work on a body at rest
			CHECK(p.Dot(v_point) >= 0.0f);
		}
	}

	/// The contact tangent space in 4D is 3 dimensional, so friction needs 3 orthonormal tangents
	/// rather than 2. Vec4::GetNormalizedPerpendicularBasis produces them; check that they really are
	/// an orthonormal basis of the hyperplane perpendicular to the normal.
	TEST_CASE("TestNormalizedPerpendicularBasisIsOrthonormal")
	{
		UnitTestRandom random(0x4DBA51);
		for (int n = 0; n < 50; ++n)
		{
			Vec4 normal = sRandomVec(random);
			if (normal.LengthSq() < 1.0e-6f)
				continue;
			normal = normal.Normalized();

			Vec4 t1, t2, t3;
			normal.GetNormalizedPerpendicularBasis(t1, t2, t3);

			// Unit length
			CHECK_APPROX_EQUAL(t1.Length(), 1.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL(t2.Length(), 1.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL(t3.Length(), 1.0f, 1.0e-4f);

			// Mutually perpendicular, and all perpendicular to the normal
			CHECK_APPROX_EQUAL(normal.Dot(t1), 0.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL(normal.Dot(t2), 0.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL(normal.Dot(t3), 0.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL(t1.Dot(t2), 0.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL(t1.Dot(t3), 0.0f, 1.0e-4f);
			CHECK_APPROX_EQUAL(t2.Dot(t3), 0.0f, 1.0e-4f);

			// {normal, t1, t2, t3} spans R^4: any vector is recovered from its 4 components
			Vec4 v = sRandomVec(random);
			Vec4 recomposed = normal * v.Dot(normal) + t1 * v.Dot(t1) + t2 * v.Dot(t2) + t3 * v.Dot(t3);
			CHECK((recomposed - v).Length() < 1.0e-4f);
		}
	}

	/// The angular row of a constraint Jacobian. For a constraint that measures the relative
	/// velocity of a point along a world space axis n, the angular Jacobian bivector J must satisfy
	/// <J, Omega> == n . (velocity of the point due to Omega), which makes J = n ^ r. The impulse
	/// conjugate to it (the angular impulse from applying scalar impulse lambda along n at r) is
	/// lambda * (n ^ r), i.e. the same bivector -- exactly as in 3D, where both are r x n.
	TEST_CASE("TestConstraintAngularJacobianRow")
	{
		UnitTestRandom random(0x4DDEF0);
		for (int n = 0; n < 50; ++n)
		{
			Vec4 r = sRandomVec(random);
			Vec4 axis = sRandomVec(random).Normalized();
			Bivec omega = sRandomBivec(random);

			Bivec jacobian = Bivec::sWedge(axis, r);

			// J . Omega reproduces the constrained component of the point velocity
			CHECK_APPROX_EQUAL(jacobian.Dot(omega), axis.Dot(omega.Contract(r)), 1.0e-5f);

			// and J^T lambda is the angular impulse of the linear impulse lambda * axis at r
			float lambda = 0.75f;
			CHECK((Bivec::sWedge(lambda * axis, r) - lambda * jacobian).IsNearZero(1.0e-10f));
		}
	}
}
