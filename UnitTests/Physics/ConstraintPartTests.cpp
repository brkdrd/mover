// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Math/Bivec.h>
#include <Jolt/Math/Rotor.h>
#include <Jolt/Math/Mat44.h>

JPH_SUPPRESS_WARNINGS_STD_BEGIN
#include <random>
JPH_SUPPRESS_WARNINGS_STD_END

// TODO(4D): Conservation tests for the 4D constraint part algebra.
//
// CLAUDE.md requires every hand derived Clifford formula to be validated against first principles
// before anything is built on it. The point constraint in particular is a genuine 4D redesign (4
// linear DOF, a 4x4 effective mass, a Vec4 multiplier and a bivector angular impulse), so it is
// validated here by the properties a correct constraint solve must have, independent of how it was
// derived: after one solve the relative velocity of the constrained points is zero, and the impulse
// pair conserves linear momentum, conserves angular momentum about the world origin, and never adds
// kinetic energy.
//
// These mirror the formulas in Physics/Constraints/ConstraintPart/PointConstraintPart.h and
// AxisConstraintPart.h. They are reproduced here rather than driving the real classes because those
// take Body references and the physics library does not link yet (CLAUDE.md, Tests section); once it
// does, drive the real parts through PhysicsTestContext and delete the duplication.
TEST_SUITE("ConstraintPartTests")
{
	/// A rigid body reduced to what a constraint part actually touches
	struct TestBody
	{
		float				mMass;
		float				mInvMass;
		Rotor				mRotation;
		Bivec				mInertiaDiagonal;			///< Principal frame (body space) bivector inertia
		Bivec				mInvInertiaDiagonal;
		Vec4				mPosition;					///< Center of mass
		Vec4				mLinearVelocity;
		Bivec				mAngularVelocity;			///< World space
		Vec4				mR;							///< World space offset from the COM to the constraint point

		/// World space inverse inertia operator: Ad(R) D^-1 Ad(R^-1), as
		/// MotionProperties::MultiplyWorldSpaceInverseInertiaByVector computes it
		Bivec				MultiplyInvInertia(BivecArg inV) const
		{
			return mRotation.RotateBivec(mInvInertiaDiagonal * mRotation.InverseRotateBivec(inV));
		}

		/// World space inertia operator (for momentum / energy bookkeeping)
		Bivec				MultiplyInertia(BivecArg inV) const
		{
			return mRotation.RotateBivec(mInertiaDiagonal * mRotation.InverseRotateBivec(inV));
		}

		/// Velocity of the constrained material point
		Vec4				GetPointVelocity() const
		{
			return mLinearVelocity + mAngularVelocity.Contract(mR);
		}

		/// Angular momentum bivector about the world origin. The wedge operands follow the torque
		/// convention (see DynamicsConventionsTest): d/dt (v ^ x) = F ^ x is the torque about the origin.
		Bivec				GetAngularMomentum() const
		{
			return MultiplyInertia(mAngularVelocity) + mMass * Bivec::sWedge(mLinearVelocity, mPosition);
		}

		Vec4				GetLinearMomentum() const
		{
			return mMass * mLinearVelocity;
		}

		float				GetKineticEnergy() const
		{
			return 0.5f * mMass * mLinearVelocity.LengthSq() + 0.5f * mAngularVelocity.Dot(MultiplyInertia(mAngularVelocity));
		}

		/// Apply an impulse at the constrained point (inSign = +1 to add, -1 to subtract)
		void				ApplyImpulse(Vec4Arg inLambda, float inSign)
		{
			mLinearVelocity += (inSign * mInvMass) * inLambda;
			mAngularVelocity += inSign * MultiplyInvInertia(Bivec::sWedge(inLambda, mR));
		}
	};

	static Rotor sRandomRotor(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> angle(-JPH_PI, JPH_PI);
		return Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), angle(ioRandom))
			 * Rotor::sRotation(Vec4::sAxisZ(), Vec4::sAxisW(), angle(ioRandom))
			 * Rotor::sRotation(Vec4::sAxisY(), Vec4::sAxisW(), angle(ioRandom));
	}

	static Vec4 sRandomVec(UnitTestRandom &ioRandom, float inRange = 1.0f)
	{
		uniform_real_distribution<float> dist(-inRange, inRange);
		return Vec4(dist(ioRandom), dist(ioRandom), dist(ioRandom), dist(ioRandom));
	}

	/// Build a dynamic test body with a random pose, velocity and (anisotropic) inertia
	static TestBody sRandomBody(UnitTestRandom &ioRandom)
	{
		uniform_real_distribution<float> positive(0.5f, 3.0f);

		TestBody b;
		b.mMass = positive(ioRandom);
		b.mInvMass = 1.0f / b.mMass;
		b.mRotation = sRandomRotor(ioRandom);
		b.mInertiaDiagonal = Bivec(positive(ioRandom), positive(ioRandom), positive(ioRandom),
								   positive(ioRandom), positive(ioRandom), positive(ioRandom));
		b.mInvInertiaDiagonal = Bivec(1.0f / b.mInertiaDiagonal.GetE12(), 1.0f / b.mInertiaDiagonal.GetE13(),
									  1.0f / b.mInertiaDiagonal.GetE14(), 1.0f / b.mInertiaDiagonal.GetE23(),
									  1.0f / b.mInertiaDiagonal.GetE24(), 1.0f / b.mInertiaDiagonal.GetE34());
		b.mPosition = sRandomVec(ioRandom, 3.0f);
		b.mLinearVelocity = sRandomVec(ioRandom, 2.0f);
		b.mAngularVelocity = Bivec(sRandomVec(ioRandom, 2.0f).GetX(), sRandomVec(ioRandom, 2.0f).GetY(),
								   sRandomVec(ioRandom, 2.0f).GetZ(), sRandomVec(ioRandom, 2.0f).GetW(),
								   sRandomVec(ioRandom, 2.0f).GetX(), sRandomVec(ioRandom, 2.0f).GetY());
		b.mR = sRandomVec(ioRandom, 1.5f);
		return b;
	}

	/// Turn a body into a static one: infinite mass, no inertia, no velocity
	static void sMakeStatic(TestBody &ioBody)
	{
		ioBody.mMass = 0.0f;					// Only used for momentum bookkeeping, which excludes static bodies
		ioBody.mInvMass = 0.0f;
		ioBody.mInertiaDiagonal = Bivec::sZero();
		ioBody.mInvInertiaDiagonal = Bivec::sZero();
		ioBody.mLinearVelocity = Vec4::sZero();
		ioBody.mAngularVelocity = Bivec::sZero();
	}

	/// The 4x4 effective mass of the point constraint, assembled exactly as
	/// PointConstraintPart::CalculateConstraintProperties does it
	static Mat44 sPointInverseEffectiveMass(const TestBody &inBody1, const TestBody &inBody2)
	{
		const Vec4 basis[4] = { Vec4::sAxisX(), Vec4::sAxisY(), Vec4::sAxisZ(), Vec4::sAxisW() };
		float summed_inv_mass = inBody1.mInvMass + inBody2.mInvMass;

		Mat44 k;
		for (uint c = 0; c < 4; ++c)
		{
			Bivec col1 = inBody1.MultiplyInvInertia(Bivec::sWedge(basis[c], inBody1.mR));
			Bivec col2 = inBody2.MultiplyInvInertia(Bivec::sWedge(basis[c], inBody2.mR));
			k.SetColumn4(c, summed_inv_mass * basis[c] + col1.Contract(inBody1.mR) + col2.Contract(inBody2.mR));
		}
		return k;
	}

	/// Place body 2 so that the two constraint points coincide, which is the state a satisfied point
	/// constraint is in (and what makes angular momentum conservation an exact statement)
	static void sAlignConstraintPoints(const TestBody &inBody1, TestBody &ioBody2)
	{
		ioBody2.mPosition = inBody1.mPosition + inBody1.mR - ioBody2.mR;
	}

	/// The inverse effective mass matrix must be symmetric positive definite: it is
	/// (invM1 + invM2) E + sum of C_r I^-1 C_r^T, and each of those is symmetric PSD.
	TEST_CASE("TestPointConstraintEffectiveMassIsSymmetricPositiveDefinite")
	{
		UnitTestRandom random(0xC0FFEE);
		for (int n = 0; n < 25; ++n)
		{
			TestBody b1 = sRandomBody(random), b2 = sRandomBody(random);
			Mat44 k = sPointInverseEffectiveMass(b1, b2);

			// Symmetry
			for (uint i = 0; i < 4; ++i)
				for (uint j = 0; j < 4; ++j)
					CHECK(abs(k(i, j) - k(j, i)) < 1.0e-4f * max(1.0f, abs(k(i, j))));

			// Positive definiteness, sampled along random directions
			for (int s = 0; s < 10; ++s)
			{
				Vec4 u = sRandomVec(random).Normalized();
				CHECK(u.Dot(Vec4(k * u)) > 0.0f);
			}

			// Invertible, so the constraint is never spuriously deactivated for two dynamic bodies
			CHECK(abs(k.GetDeterminant()) > 1.0e-12f);
		}
	}

	/// The defining property of the velocity solve: one exact solve removes the relative velocity of
	/// the two constrained points along all 4 spatial directions. Every sign in the derivation (the
	/// Contract in the Jacobian, the wedge in the impulse, and the assembly of K from both) has to be
	/// consistent for this to hold -- flipping any one of them leaves a non zero residual.
	TEST_CASE("TestPointConstraintSolveRemovesRelativePointVelocity")
	{
		UnitTestRandom random(0xC0FFE1);
		for (int n = 0; n < 50; ++n)
		{
			TestBody b1 = sRandomBody(random), b2 = sRandomBody(random);
			sAlignConstraintPoints(b1, b2);

			Mat44 effective_mass = sPointInverseEffectiveMass(b1, b2).Inversed();
			Vec4 jv = b1.GetPointVelocity() - b2.GetPointVelocity();
			Vec4 lambda = effective_mass * jv;

			b1.ApplyImpulse(lambda, -1.0f);
			b2.ApplyImpulse(lambda, 1.0f);

			Vec4 residual = b1.GetPointVelocity() - b2.GetPointVelocity();
			CHECK(residual.Length() < 1.0e-3f * max(1.0f, jv.Length()));
		}
	}

	/// Same, with body 2 static (the common case: a body pinned to the world)
	TEST_CASE("TestPointConstraintSolveAgainstStaticBody")
	{
		UnitTestRandom random(0xC0FFE2);
		for (int n = 0; n < 50; ++n)
		{
			TestBody b1 = sRandomBody(random), b2 = sRandomBody(random);
			sMakeStatic(b2);
			sAlignConstraintPoints(b1, b2);

			Mat44 k = sPointInverseEffectiveMass(b1, b2);
			CHECK(abs(k.GetDeterminant()) > 1.0e-12f);

			Vec4 jv = b1.GetPointVelocity() - b2.GetPointVelocity();
			Vec4 lambda = k.Inversed() * jv;

			b1.ApplyImpulse(lambda, -1.0f);

			// The point of the dynamic body is brought to rest (the static body's point never moves)
			CHECK(b1.GetPointVelocity().Length() < 1.0e-3f * max(1.0f, jv.Length()));
		}
	}

	/// The constraint impulse is internal: it must not change the total linear or angular momentum of
	/// the pair, and it must not create kinetic energy. Angular momentum is the sharp one -- the two
	/// impulses are equal and opposite but act through the bodies' own offsets r1 and r2, so the
	/// bookkeeping only closes if the angular impulse operator is the exact adjoint of the point
	/// velocity operator.
	TEST_CASE("TestPointConstraintConservesMomentumAndEnergy")
	{
		UnitTestRandom random(0xC0FFE3);
		for (int n = 0; n < 50; ++n)
		{
			TestBody b1 = sRandomBody(random), b2 = sRandomBody(random);
			sAlignConstraintPoints(b1, b2);

			Vec4 linear_before = b1.GetLinearMomentum() + b2.GetLinearMomentum();
			Bivec angular_before = b1.GetAngularMomentum() + b2.GetAngularMomentum();
			float energy_before = b1.GetKineticEnergy() + b2.GetKineticEnergy();

			Vec4 jv = b1.GetPointVelocity() - b2.GetPointVelocity();
			Vec4 lambda = sPointInverseEffectiveMass(b1, b2).Inversed() * jv;
			b1.ApplyImpulse(lambda, -1.0f);
			b2.ApplyImpulse(lambda, 1.0f);

			Vec4 linear_after = b1.GetLinearMomentum() + b2.GetLinearMomentum();
			Bivec angular_after = b1.GetAngularMomentum() + b2.GetAngularMomentum();
			float energy_after = b1.GetKineticEnergy() + b2.GetKineticEnergy();

			CHECK((linear_after - linear_before).Length() < 1.0e-3f * max(1.0f, linear_before.Length()));
			CHECK((angular_after - angular_before).Length() < 1.0e-2f * max(1.0f, angular_before.Length()));

			// A constraint that only removes relative velocity can only dissipate energy
			CHECK(energy_after <= energy_before + 1.0e-3f * max(1.0f, energy_before));
		}
	}

	/// The single row (AxisConstraintPart) case: the same machinery restricted to one direction. The
	/// effective mass is the scalar invM1 + invM2 + (n^r1).(I1^-1 (n^r1)) + (n^r2).(I2^-1 (n^r2)) and
	/// the solve must remove the relative point velocity *along that axis only*.
	TEST_CASE("TestAxisConstraintSolveRemovesRelativeVelocityAlongAxis")
	{
		UnitTestRandom random(0xC0FFE4);
		for (int n = 0; n < 50; ++n)
		{
			TestBody b1 = sRandomBody(random), b2 = sRandomBody(random);
			sAlignConstraintPoints(b1, b2);
			Vec4 axis = sRandomVec(random).Normalized();

			// Angular Jacobian rows: n ^ r, not r ^ n
			Bivec j1 = Bivec::sWedge(axis, b1.mR);
			Bivec j2 = Bivec::sWedge(axis, b2.mR);
			Bivec inv_i1_j1 = b1.MultiplyInvInertia(j1);
			Bivec inv_i2_j2 = b2.MultiplyInvInertia(j2);

			float inv_effective_mass = b1.mInvMass + b2.mInvMass + j1.Dot(inv_i1_j1) + j2.Dot(inv_i2_j2);
			CHECK(inv_effective_mass > 0.0f);

			float jv = axis.Dot(b1.GetPointVelocity() - b2.GetPointVelocity());
			float lambda = jv / inv_effective_mass;

			Vec4 linear_before = b1.GetLinearMomentum() + b2.GetLinearMomentum();
			Bivec angular_before = b1.GetAngularMomentum() + b2.GetAngularMomentum();

			b1.ApplyImpulse(lambda * axis, -1.0f);
			b2.ApplyImpulse(lambda * axis, 1.0f);

			// Relative velocity along the axis is gone
			CHECK(abs(axis.Dot(b1.GetPointVelocity() - b2.GetPointVelocity())) < 1.0e-3f * max(1.0f, abs(jv)));

			// and the impulse pair was internal
			CHECK((b1.GetLinearMomentum() + b2.GetLinearMomentum() - linear_before).Length() < 1.0e-3f * max(1.0f, linear_before.Length()));
			Bivec angular_after = b1.GetAngularMomentum() + b2.GetAngularMomentum();
			CHECK((angular_after - angular_before).Length() < 1.0e-2f * max(1.0f, angular_before.Length()));
		}
	}

	/// The angular Jacobian row and the angular impulse must be the *same* bivector, which is what
	/// makes J M^-1 J^T (the effective mass) the correct scaling for the solve. This is the identity
	/// that fails if r ^ n is used instead of n ^ r: the effective mass is unchanged (it is quadratic)
	/// but the solve then converges to the wrong velocity state.
	TEST_CASE("TestAxisConstraintEffectiveMassPredictsVelocityChange")
	{
		UnitTestRandom random(0xC0FFE5);
		for (int n = 0; n < 50; ++n)
		{
			TestBody b1 = sRandomBody(random), b2 = sRandomBody(random);
			sAlignConstraintPoints(b1, b2);
			Vec4 axis = sRandomVec(random).Normalized();

			Bivec j1 = Bivec::sWedge(axis, b1.mR);
			Bivec j2 = Bivec::sWedge(axis, b2.mR);
			float inv_effective_mass = b1.mInvMass + b2.mInvMass
									 + j1.Dot(b1.MultiplyInvInertia(j1)) + j2.Dot(b2.MultiplyInvInertia(j2));

			// Applying a unit impulse pair along the axis must change the relative velocity along the
			// axis by exactly the inverse effective mass
			float before = axis.Dot(b1.GetPointVelocity() - b2.GetPointVelocity());
			b1.ApplyImpulse(axis, -1.0f);
			b2.ApplyImpulse(axis, 1.0f);
			float after = axis.Dot(b1.GetPointVelocity() - b2.GetPointVelocity());

			CHECK_APPROX_EQUAL(before - after, inv_effective_mass, 1.0e-3f * max(1.0f, inv_effective_mass));
		}
	}
}
