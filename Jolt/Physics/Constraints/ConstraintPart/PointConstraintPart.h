// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Body/Body.h>
#include <Jolt/Physics/StateRecorder.h>

JPH_NAMESPACE_BEGIN

/// Constrains movement along all 4 spatial axes (4D version)
///
/// @see "Constraints Derivation for Rigid Body Simulation in 3D" - Daniel Chappuis, section 2.2.1
///
/// Constraint equation (eq 45):
///
/// \f[C = p_2 - p_1\f]
///
/// TODO(4D): this is a redesign of the 3D part rather than a type swap. A point-to-point constraint
/// has 4 linear degrees of freedom in 4D, so C has 4 rows, the effective mass is a 4x4 matrix and the
/// Lagrange multiplier is a Vec4 (it was 3 / 3x3 / Vec3). The two operators that replace the 3D cross
/// product are adjoint to each other up to a sign:
///
///   - the velocity of a body fixed point due to the angular velocity bivector w is w _| r
///     (Bivec::Contract), and
///   - the angular impulse of a linear impulse P applied at offset r is P ^ r (Bivec::sWedge),
///
/// with \f$\left<w \angle r, P\right> = -\left<w, r \wedge P\right> = \left<w, P \wedge r\right>\f$.
/// Contract is the *negative* adjoint of the wedge, which is why the operands of the wedge are
/// swapped relative to the 3D "r x P". Getting this backwards silently injects energy; it is pinned
/// by UnitTests/Physics/DynamicsConventionsTest.cpp.
///
/// Writing \f$C_r\f$ for the linear map \f$w \mapsto w \angle r\f$ (4x6) its transpose is
/// \f$C_r^T: u \mapsto u \wedge r\f$ (6x4), so the effective mass is
///
/// \f[K = \left(m_1^{-1} + m_2^{-1}\right) E + C_{r_1} I_1^{-1} C_{r_1}^T + C_{r_2} I_2^{-1} C_{r_2}^T\f]
///
/// which we build one column at a time by feeding it the 4 basis vectors. The per body 6x4 operator
/// \f$I^{-1} C_r^T\f$ is cached as its 4 columns (a bivector each), and is exactly what turns a
/// Lagrange multiplier into an angular velocity change.
///
/// Used terms (everything in world space unless noted):\n
/// p1, p2 = constraint points.\n
/// r1 = p1 - x1, r2 = p2 - x2.\n
/// x1, x2 = center of mass for the bodies.\n
/// v = [v1, w1, v2, w2], with w1, w2 bivector angular velocities.\n
/// M = mass matrix, block diagonal [m1, I1, m2, I2] with I the 6x6 bivector inertia.\n
/// \f$K^{-1} = \left( J M^{-1} J^T \right)^{-1}\f$ = effective mass.\n
/// \f$\beta\f$ = baumgarte constant.\n
/// E = identity matrix.
class PointConstraintPart
{
	/// Apply the cached 6x4 operator I^-1 C_r^T to a Lagrange multiplier: returns I^-1 (lambda ^ r)
	static JPH_INLINE Bivec	sApplyInvInertiaWedge(const Bivec (&inColumns)[4], Vec4Arg inLambda)
	{
		return inLambda.GetX() * inColumns[0] + inLambda.GetY() * inColumns[1]
			 + inLambda.GetZ() * inColumns[2] + inLambda.GetW() * inColumns[3];
	}

	/// Cache the 4 columns of I^-1 C_r^T, i.e. column j = I^-1 (e_j ^ r)
	static JPH_INLINE void	sCalculateInvInertiaWedge(const MotionProperties *inMP, RotorArg inRotation, Vec4Arg inR, Bivec (&outColumns)[4])
	{
		outColumns[0] = inMP->MultiplyWorldSpaceInverseInertiaByVector(inRotation, Bivec::sWedge(Vec4::sAxisX(), inR));
		outColumns[1] = inMP->MultiplyWorldSpaceInverseInertiaByVector(inRotation, Bivec::sWedge(Vec4::sAxisY(), inR));
		outColumns[2] = inMP->MultiplyWorldSpaceInverseInertiaByVector(inRotation, Bivec::sWedge(Vec4::sAxisZ(), inR));
		outColumns[3] = inMP->MultiplyWorldSpaceInverseInertiaByVector(inRotation, Bivec::sWedge(Vec4::sAxisW(), inR));
	}

	static JPH_INLINE void	sZeroColumns(Bivec (&outColumns)[4])
	{
		for (int i = 0; i < 4; ++i)
			outColumns[i] = Bivec::sZero();
	}

	JPH_INLINE bool			ApplyVelocityStep(Body &ioBody1, Body &ioBody2, Vec4Arg inLambda) const
	{
		// Apply impulse if delta is not zero
		if (inLambda != Vec4::sZero())
		{
			// Calculate velocity change due to constraint
			//
			// Impulse:
			// P = J^T lambda
			//
			// Euler velocity integration:
			// v' = v + M^-1 P
			if (ioBody1.IsDynamic())
			{
				MotionProperties *mp1 = ioBody1.GetMotionProperties();
				mp1->SubLinearVelocityStep(mp1->GetInverseMass() * inLambda);
				mp1->SubAngularVelocityStep(sApplyInvInertiaWedge(mInvI1_R1X, inLambda));
			}
			if (ioBody2.IsDynamic())
			{
				MotionProperties *mp2 = ioBody2.GetMotionProperties();
				mp2->AddLinearVelocityStep(mp2->GetInverseMass() * inLambda);
				mp2->AddAngularVelocityStep(sApplyInvInertiaWedge(mInvI2_R2X, inLambda));
			}
			return true;
		}

		return false;
	}

public:
	/// Calculate properties used during the functions below
	/// @param inBody1 The first body that this constraint is attached to
	/// @param inBody2 The second body that this constraint is attached to
	/// @param inRotation1 The rotor for body 1 (rotates inR1 to world space and orients the inertia)
	/// @param inRotation2 The rotor for body 2 (rotates inR2 to world space and orients the inertia)
	/// @param inR1 Local space vector from center of mass to constraint point for body 1
	/// @param inR2 Local space vector from center of mass to constraint point for body 2
	inline void				CalculateConstraintProperties(const Body &inBody1, RotorArg inRotation1, Vec4Arg inR1, const Body &inBody2, RotorArg inRotation2, Vec4Arg inR2)
	{
		// Positions where the point constraint acts on, relative to the centers of mass, in world space
		mR1 = inRotation1 * inR1;
		mR2 = inRotation2 * inR2;

		// Calculate effective mass: K^-1 = (J M^-1 J^T)^-1, one column at a time
		float summed_inv_mass = 0.0f;
		if (inBody1.IsDynamic())
		{
			const MotionProperties *mp1 = inBody1.GetMotionProperties();
			summed_inv_mass = mp1->GetInverseMass();
			sCalculateInvInertiaWedge(mp1, inRotation1, mR1, mInvI1_R1X);
		}
		else
		{
			sZeroColumns(mInvI1_R1X);
		}

		if (inBody2.IsDynamic())
		{
			const MotionProperties *mp2 = inBody2.GetMotionProperties();
			summed_inv_mass += mp2->GetInverseMass();
			sCalculateInvInertiaWedge(mp2, inRotation2, mR2, mInvI2_R2X);
		}
		else
		{
			sZeroColumns(mInvI2_R2X);
		}

		// Column j is K e_j = (invM1 + invM2) e_j + (I1^-1 (e_j ^ r1)) _| r1 + (I2^-1 (e_j ^ r2)) _| r2.
		// The two angular terms are C_r (I^-1 C_r^T), which is symmetric positive semi definite, so K
		// is symmetric and (whenever a body is dynamic, giving a positive invM on the diagonal)
		// positive definite.
		const Vec4 basis[4] = { Vec4::sAxisX(), Vec4::sAxisY(), Vec4::sAxisZ(), Vec4::sAxisW() };
		Mat44 inv_effective_mass;
		for (uint c = 0; c < 4; ++c)
		{
			Vec4 column = summed_inv_mass * basis[c] + mInvI1_R1X[c].Contract(mR1) + mInvI2_R2X[c].Contract(mR2);
			inv_effective_mass.SetColumn4(c, column);
		}

		// A singular K means the constraint has no effect along at least one direction (e.g. both
		// bodies non dynamic), in which case we disable it rather than divide by zero
		if (abs(inv_effective_mass.GetDeterminant()) < 1.0e-24f)
			Deactivate();
		else
			mEffectiveMass = inv_effective_mass.Inversed();
	}

	/// Deactivate this constraint
	inline void				Deactivate()
	{
		mEffectiveMass = Mat44::sZero();
		mTotalLambda = Vec4::sZero();
	}

	/// Check if constraint is active
	inline bool				IsActive() const
	{
		return mEffectiveMass(3, 3) != 0.0f;
	}

	/// Must be called from the WarmStartVelocityConstraint call to apply the previous frame's impulses
	/// @param ioBody1 The first body that this constraint is attached to
	/// @param ioBody2 The second body that this constraint is attached to
	/// @param inWarmStartImpulseRatio Ratio of new step to old time step (dt_new / dt_old) for scaling the lagrange multiplier of the previous frame
	inline void				WarmStart(Body &ioBody1, Body &ioBody2, float inWarmStartImpulseRatio)
	{
		mTotalLambda *= inWarmStartImpulseRatio;
		ApplyVelocityStep(ioBody1, ioBody2, mTotalLambda);
	}

	/// Iteratively update the velocity constraint. Makes sure d/dt C(...) = 0, where C is the constraint equation.
	/// @param ioBody1 The first body that this constraint is attached to
	/// @param ioBody2 The second body that this constraint is attached to
	inline bool				SolveVelocityConstraint(Body &ioBody1, Body &ioBody2)
	{
		// Calculate lagrange multiplier:
		//
		// lambda = -K^-1 (J v + b)
		//
		// J v is the difference of the two constraint point velocities, v + w _| r per body
		Vec4 jv = ioBody1.GetLinearVelocity() + ioBody1.GetAngularVelocity().Contract(mR1)
				- ioBody2.GetLinearVelocity() - ioBody2.GetAngularVelocity().Contract(mR2);
		Vec4 lambda = mEffectiveMass * jv;
		mTotalLambda += lambda; // Store accumulated lambda
		return ApplyVelocityStep(ioBody1, ioBody2, lambda);
	}

	/// Iteratively update the position constraint. Makes sure C(...) = 0.
	/// @param ioBody1 The first body that this constraint is attached to
	/// @param ioBody2 The second body that this constraint is attached to
	/// @param inBaumgarte Baumgarte constant (fraction of the error to correct)
	inline bool				SolvePositionConstraint(Body &ioBody1, Body &ioBody2, float inBaumgarte) const
	{
		Vec4 separation = Vec4(ioBody2.GetCenterOfMassPosition() - ioBody1.GetCenterOfMassPosition()) + mR2 - mR1;
		if (separation != Vec4::sZero())
		{
			// Calculate lagrange multiplier (lambda) for Baumgarte stabilization:
			//
			// lambda = -K^-1 * beta / dt * C
			//
			// We should divide by inDeltaTime, but we should multiply by inDeltaTime in the Euler step below so they're cancelled out
			Vec4 lambda = mEffectiveMass * (-inBaumgarte * separation);

			// Directly integrate velocity change for one time step
			//
			// Euler velocity integration:
			// dv = M^-1 P
			//
			// Impulse:
			// P = J^T lambda
			//
			// Euler position integration:
			// x' = x + dv * dt
			//
			// Note we don't accumulate velocities for the stabilization. This is using the approach described in 'Modeling and
			// Solving Constraints' by Erin Catto presented at GDC 2007. On slide 78 it is suggested to split up the Baumgarte
			// stabilization for positional drift so that it does not actually add to the momentum. We combine an Euler velocity
			// integrate + a position integrate and then discard the velocity change.
			if (ioBody1.IsDynamic())
			{
				ioBody1.SubPositionStep(ioBody1.GetMotionProperties()->GetInverseMass() * lambda);
				ioBody1.SubRotationStep(sApplyInvInertiaWedge(mInvI1_R1X, lambda));
			}
			if (ioBody2.IsDynamic())
			{
				ioBody2.AddPositionStep(ioBody2.GetMotionProperties()->GetInverseMass() * lambda);
				ioBody2.AddRotationStep(sApplyInvInertiaWedge(mInvI2_R2X, lambda));
			}

			return true;
		}

		return false;
	}

	/// Return lagrange multiplier
	Vec4					GetTotalLambda() const
	{
		return mTotalLambda;
	}

	/// Save state of this constraint part
	void					SaveState(StateRecorder &inStream) const
	{
		inStream.Write(mTotalLambda);
	}

	/// Restore state of this constraint part
	void					RestoreState(StateRecorder &inStream)
	{
		inStream.Read(mTotalLambda);
	}

private:
	Vec4					mR1;
	Vec4					mR2;
	Bivec					mInvI1_R1X[4];			///< Columns of I1^-1 C_r1^T: column j = I1^-1 (e_j ^ r1)
	Bivec					mInvI2_R2X[4];			///< Columns of I2^-1 C_r2^T: column j = I2^-1 (e_j ^ r2)
	Mat44					mEffectiveMass;
	Vec4					mTotalLambda { Vec4::sZero() };
};

JPH_NAMESPACE_END
