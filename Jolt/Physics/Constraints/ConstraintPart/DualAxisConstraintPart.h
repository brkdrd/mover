// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Body/Body.h>
#include <Jolt/Physics/StateRecorder.h>
#include <Jolt/Math/Vector.h>
#include <Jolt/Math/Matrix.h>

JPH_NAMESPACE_BEGIN

/**
	Constrains movement on 2 axis (4D version)

	@see "Constraints Derivation for Rigid Body Simulation in 3D" - Daniel Chappuis, section 2.3.1

	Constraint equation (eq 51): C = [ (p2 - p1) . n1 ; (p2 - p1) . n2 ]

	In 4D the angular Jacobian rows use the wedge product n ^ r (a bivector) in place of r x n -- the
	operands are swapped relative to the 3D cross product, see AxisConstraintPart. The angular
	velocities w1/w2 are bivectors, and the world-space inverse inertia is applied via
	MotionProperties::MultiplyWorldSpaceInverseInertiaByVector. The effective mass stays a 2x2 matrix
	(one row per constrained axis, independent of the spatial dimension).

	n1, n2 = constraint axes (normalized), p1/p2 = constraint points, r1 = p1 - x1, r2 = p2 - x2,
	u = p2 - p1, x1/x2 = center of mass. v = [v1, w1, v2, w2].
**/
class DualAxisConstraintPart
{
public:
	using Vec2 = Vector<2>;
	using Mat22 = Matrix<2, 2>;

private:
	/// Internal helper function to update velocities of bodies after Lagrange multiplier is calculated
	JPH_INLINE bool				ApplyVelocityStep(Body &ioBody1, Body &ioBody2, Vec4Arg inN1, Vec4Arg inN2, const Vec2 &inLambda) const
	{
		// Apply impulse if delta is not zero
		if (!inLambda.IsZero())
		{
			// Impulse P = J^T lambda, Euler velocity integration v' = v + M^-1 P
			Vec4 impulse = inN1 * inLambda[0] + inN2 * inLambda[1];
			if (ioBody1.IsDynamic())
			{
				MotionProperties *mp1 = ioBody1.GetMotionProperties();
				mp1->SubLinearVelocityStep(mp1->GetInverseMass() * impulse);
				mp1->SubAngularVelocityStep(mInvI1_R1PlusUxN1 * inLambda[0] + mInvI1_R1PlusUxN2 * inLambda[1]);
			}
			if (ioBody2.IsDynamic())
			{
				MotionProperties *mp2 = ioBody2.GetMotionProperties();
				mp2->AddLinearVelocityStep(mp2->GetInverseMass() * impulse);
				mp2->AddAngularVelocityStep(mInvI2_R2xN1 * inLambda[0] + mInvI2_R2xN2 * inLambda[1]);
			}
			return true;
		}

		return false;
	}

	/// Internal helper function to calculate the lagrange multiplier
	inline void					CalculateLagrangeMultiplier(const Body &inBody1, const Body &inBody2, Vec4Arg inN1, Vec4Arg inN2, Vec2 &outLambda) const
	{
		// lambda = -K^-1 (J v + b)
		Vec4 delta_lin = inBody1.GetLinearVelocity() - inBody2.GetLinearVelocity();
		Vec2 jv;
		jv[0] = inN1.Dot(delta_lin) + mR1PlusUxN1.Dot(inBody1.GetAngularVelocity()) - mR2xN1.Dot(inBody2.GetAngularVelocity());
		jv[1] = inN2.Dot(delta_lin) + mR1PlusUxN2.Dot(inBody1.GetAngularVelocity()) - mR2xN2.Dot(inBody2.GetAngularVelocity());
		outLambda = mEffectiveMass * jv;
	}

public:
	/// Calculate properties used during the functions below
	/// All input vectors are in world space; inRotation1/2 are the bodies' rotors (for the inverse inertia).
	inline void					CalculateConstraintProperties(const Body &inBody1, RotorArg inRotation1, Vec4Arg inR1PlusU, const Body &inBody2, RotorArg inRotation2, Vec4Arg inR2, Vec4Arg inN1, Vec4Arg inN2)
	{
		JPH_ASSERT(inN1.IsNormalized(1.0e-5f));
		JPH_ASSERT(inN2.IsNormalized(1.0e-5f));

		// Angular Jacobian bivectors n ^ r (operands swapped w.r.t. the 3D r x n, see
		// AxisConstraintPart / UnitTests/Physics/DynamicsConventionsTest)
		mR1PlusUxN1 = Bivec::sWedge(inN1, inR1PlusU);
		mR1PlusUxN2 = Bivec::sWedge(inN2, inR1PlusU);
		mR2xN1 = Bivec::sWedge(inN1, inR2);
		mR2xN2 = Bivec::sWedge(inN2, inR2);

		// Calculate effective mass: K^-1 = (J M^-1 J^T)^-1, eq 59
		Mat22 inv_effective_mass;
		if (inBody1.IsDynamic())
		{
			const MotionProperties *mp1 = inBody1.GetMotionProperties();
			mInvI1_R1PlusUxN1 = mp1->MultiplyWorldSpaceInverseInertiaByVector(inRotation1, mR1PlusUxN1);
			mInvI1_R1PlusUxN2 = mp1->MultiplyWorldSpaceInverseInertiaByVector(inRotation1, mR1PlusUxN2);

			inv_effective_mass(0, 0) = mp1->GetInverseMass() + mR1PlusUxN1.Dot(mInvI1_R1PlusUxN1);
			inv_effective_mass(0, 1) = mR1PlusUxN1.Dot(mInvI1_R1PlusUxN2);
			inv_effective_mass(1, 0) = mR1PlusUxN2.Dot(mInvI1_R1PlusUxN1);
			inv_effective_mass(1, 1) = mp1->GetInverseMass() + mR1PlusUxN2.Dot(mInvI1_R1PlusUxN2);
		}
		else
		{
			inv_effective_mass = Mat22::sZero();
		}

		if (inBody2.IsDynamic())
		{
			const MotionProperties *mp2 = inBody2.GetMotionProperties();
			mInvI2_R2xN1 = mp2->MultiplyWorldSpaceInverseInertiaByVector(inRotation2, mR2xN1);
			mInvI2_R2xN2 = mp2->MultiplyWorldSpaceInverseInertiaByVector(inRotation2, mR2xN2);

			inv_effective_mass(0, 0) += mp2->GetInverseMass() + mR2xN1.Dot(mInvI2_R2xN1);
			inv_effective_mass(0, 1) += mR2xN1.Dot(mInvI2_R2xN2);
			inv_effective_mass(1, 0) += mR2xN2.Dot(mInvI2_R2xN1);
			inv_effective_mass(1, 1) += mp2->GetInverseMass() + mR2xN2.Dot(mInvI2_R2xN2);
		}

		if (!mEffectiveMass.SetInversed(inv_effective_mass))
			Deactivate();
	}

	/// Deactivate this constraint
	inline void					Deactivate()
	{
		mEffectiveMass.SetZero();
		mTotalLambda.SetZero();
	}

	/// Check if constraint is active
	inline bool					IsActive() const
	{
		return !mEffectiveMass.IsZero();
	}

	/// Must be called from the WarmStartVelocityConstraint call to apply the previous frame's impulses
	inline void					WarmStart(Body &ioBody1, Body &ioBody2, Vec4Arg inN1, Vec4Arg inN2, float inWarmStartImpulseRatio)
	{
		mTotalLambda *= inWarmStartImpulseRatio;
		ApplyVelocityStep(ioBody1, ioBody2, inN1, inN2, mTotalLambda);
	}

	/// Iteratively update the velocity constraint. Makes sure d/dt C(...) = 0, where C is the constraint equation.
	inline bool					SolveVelocityConstraint(Body &ioBody1, Body &ioBody2, Vec4Arg inN1, Vec4Arg inN2)
	{
		Vec2 lambda;
		CalculateLagrangeMultiplier(ioBody1, ioBody2, inN1, inN2, lambda);

		// Store accumulated lambda
		mTotalLambda += lambda;

		return ApplyVelocityStep(ioBody1, ioBody2, inN1, inN2, lambda);
	}

	/// Iteratively update the position constraint. Makes sure C(...) = 0.
	inline bool					SolvePositionConstraint(Body &ioBody1, Body &ioBody2, Vec4Arg inU, Vec4Arg inN1, Vec4Arg inN2, float inBaumgarte) const
	{
		Vec2 c;
		c[0] = inU.Dot(inN1);
		c[1] = inU.Dot(inN2);
		if (!c.IsZero())
		{
			// lambda = -K^-1 * beta / dt * C (the dt cancels with the position integrate below)
			Vec2 lambda = -inBaumgarte * (mEffectiveMass * c);

			// Directly integrate the position change for one time step, discarding the velocity change (Catto GDC 2007)
			Vec4 impulse = inN1 * lambda[0] + inN2 * lambda[1];
			if (ioBody1.IsDynamic())
			{
				ioBody1.SubPositionStep(ioBody1.GetMotionProperties()->GetInverseMass() * impulse);
				ioBody1.SubRotationStep(mInvI1_R1PlusUxN1 * lambda[0] + mInvI1_R1PlusUxN2 * lambda[1]);
			}
			if (ioBody2.IsDynamic())
			{
				ioBody2.AddPositionStep(ioBody2.GetMotionProperties()->GetInverseMass() * impulse);
				ioBody2.AddRotationStep(mInvI2_R2xN1 * lambda[0] + mInvI2_R2xN2 * lambda[1]);
			}
			return true;
		}

		return false;
	}

	/// Override total lagrange multiplier, can be used to set the initial value for warm starting
	inline void					SetTotalLambda(const Vec2 &inLambda)
	{
		mTotalLambda = inLambda;
	}

	/// Return lagrange multiplier
	inline const Vec2 &			GetTotalLambda() const
	{
		return mTotalLambda;
	}

	/// Save state of this constraint part
	void						SaveState(StateRecorder &inStream) const
	{
		inStream.Write(mTotalLambda);
	}

	/// Restore state of this constraint part
	void						RestoreState(StateRecorder &inStream)
	{
		inStream.Read(mTotalLambda);
	}

private:
	Bivec						mR1PlusUxN1;
	Bivec						mR1PlusUxN2;
	Bivec						mR2xN1;
	Bivec						mR2xN2;
	Bivec						mInvI1_R1PlusUxN1;
	Bivec						mInvI1_R1PlusUxN2;
	Bivec						mInvI2_R2xN1;
	Bivec						mInvI2_R2xN2;
	Mat22						mEffectiveMass;
	Vec2						mTotalLambda { Vec2::sZero() };
};

JPH_NAMESPACE_END
