// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

void MotionProperties::MoveKinematic(Vec4Arg inDeltaPosition, RotorArg inDeltaRotation, float inDeltaTime)
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite));
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::Read));
	JPH_ASSERT(mCachedBodyType == EBodyType::RigidBody);
	JPH_ASSERT(mCachedMotionType != EMotionType::Static);

	// Calculate required linear velocity
	mLinearVelocity = LockTranslation(inDeltaPosition / inDeltaTime);

	// Calculate required angular velocity from rotor delta
	// For small rotations: R ≈ 1 + dt/2 * Omega, so Omega ≈ 2 * bivec(R) / dt
	float inv_dt = 1.0f / inDeltaTime;
	mAngularVelocity = LockAngular(Bivec(
		2.0f * inv_dt * inDeltaRotation.GetE12(),
		2.0f * inv_dt * inDeltaRotation.GetE13(),
		2.0f * inv_dt * inDeltaRotation.GetE14(),
		2.0f * inv_dt * inDeltaRotation.GetE23(),
		2.0f * inv_dt * inDeltaRotation.GetE24(),
		2.0f * inv_dt * inDeltaRotation.GetE34()));
}

void MotionProperties::ClampLinearVelocity()
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite));

	float len_sq = mLinearVelocity.LengthSq();
	JPH_ASSERT(isfinite(len_sq));
	if (len_sq > Square(mMaxLinearVelocity))
		mLinearVelocity *= mMaxLinearVelocity / sqrt(len_sq);
}

void MotionProperties::ClampAngularVelocity()
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite));

	float len_sq = mAngularVelocity.LengthSq();
	JPH_ASSERT(isfinite(len_sq));
	if (len_sq > Square(mMaxAngularVelocity))
		mAngularVelocity *= mMaxAngularVelocity / sqrt(len_sq);
}

inline void MotionProperties::ScaleToMass(float inMass)
{
	JPH_ASSERT(mInvMass > 0.0f, "Body must have finite mass");
	JPH_ASSERT(inMass > 0.0f, "New mass cannot be zero");

	float new_inv_mass = 1.0f / inMass;
	mInvInertiaDiagonal *= new_inv_mass / mInvMass;
	mInvMass = new_inv_mass;
}

Bivec MotionProperties::MultiplyWorldSpaceInverseInertiaByVector(RotorArg inBodyRotation, BivecArg inV) const
{
	JPH_ASSERT(mCachedMotionType == EMotionType::Dynamic);

	// Mask out components of DOFs that are not allowed
	Bivec angular_dofs_mask = GetAngularDOFsBivecMask();
	Bivec v = Bivec::sAnd(inV, angular_dofs_mask);

	// The world-space inverse inertia operator is I_world^{-1} = Ad(R) D^{-1} Ad(R)^{-1}, where D^{-1}
	// is the diagonal inverse inertia in the principal frame and R = body_rotation * inertia_rotation
	// maps the principal frame into world space. Applied to a world-space bivector:
	//   1. bring it into the principal frame:  v_local  = Ad(R^{-1}) v_world
	//   2. apply the diagonal inverse inertia:  a_local  = D^{-1} v_local
	//   3. bring the result back to world space: a_world = Ad(R) a_local
	Rotor combined = inBodyRotation * mInertiaRotation;

	Bivec v_local = combined.InverseRotateBivec(v);
	Bivec a_local = mInvInertiaDiagonal * v_local;
	Bivec result = combined.RotateBivec(a_local);

	// Mask out rows of DOFs that are not allowed
	return Bivec::sAnd(result, angular_dofs_mask);
}

InverseInertiaOperator MotionProperties::GetInverseInertiaForRotation(RotorArg inBodyRotation, float inInvInertiaScale) const
{
	JPH_ASSERT(mCachedMotionType == EMotionType::Dynamic);

	// I_world^-1 = Ad(R) D^-1 Ad(R)^T with R = body_rotation * inertia_rotation. Ad(R) is orthogonal,
	// so Ad(R)^-1 = Ad(R)^T and the operator is symmetric:
	//   col_k = Ad(R) D^-1 Ad(R)^T e_k = sum_j d_j * (A_j . e_k) * A_j,   A_j = Ad(R) e_j
	// Ad(R)(e_i ^ e_j) is (R e_i ~R) ^ (R e_j ~R), so all six A_j come from the same four sandwich
	// products - the whole cost of RotateBivec, paid once for the entire operator.
	Rotor combined = inBodyRotation * mInertiaRotation;
	Vec4 f1 = combined * Vec4::sAxisX();
	Vec4 f2 = combined * Vec4::sAxisY();
	Vec4 f3 = combined * Vec4::sAxisZ();
	Vec4 f4 = combined * Vec4::sAxisW();

	Bivec a[6] = { Bivec::sWedge(f1, f2), Bivec::sWedge(f1, f3), Bivec::sWedge(f1, f4),
				   Bivec::sWedge(f2, f3), Bivec::sWedge(f2, f4), Bivec::sWedge(f3, f4) };

	// DOF locks are in world space, whereas the diagonal is in the principal frame.
	// Mask columns and rows of the assembled world operator, not the local diagonal.
	Bivec mask = GetAngularDOFsBivecMask();
	Bivec inv_d = mInvInertiaDiagonal * inInvInertiaScale;

	InverseInertiaOperator op;
	for (int k = 0; k < 6; ++k)
	{
		if ((uint16(mAllowedDOFs) & (1u << (k + 4))) == 0)
		{
			op.mColumn[k] = Bivec::sZero();
			continue;
		}
		Vec8 col = Vec8::sZero();
		for (int j = 0; j < 6; ++j)
			col = col + a[j].GetVec8() * Vec8::sReplicate(inv_d[j] * a[j][k]);

		// Rows of locked DOFs have to come out zero too; the operator is symmetric, so masking the
		// column is what masking the row would do
		op.mColumn[k] = Bivec::sAnd(Bivec(col), mask);
	}

	return op;
}

namespace MotionPropertiesImpl
{
	// Step-local double adjoint. A float rotor is only approximately unit: four
	// successive float adjoints need not be inverses and can amplify momentum
	// even for dt = 0. Evaluate the normalized SO(4) map in double, and round
	// only the final velocity. This normalizes the ROTATION, never L or energy.
	class GyroscopicAdjoint
	{
	public:
		explicit GyroscopicAdjoint(RotorArg inRotation)
		{
			double s = inRotation.GetScalar(), a = inRotation.GetE12(), b = inRotation.GetE13(), c = inRotation.GetE14();
			double d = inRotation.GetE23(), e = inRotation.GetE24(), f = inRotation.GetE34(), p = inRotation.GetPseudoscalar();
			double norm = s*s + a*a + b*b + c*c + d*d + e*e + f*f + p*p;
			double defect = 2 * (s*p - a*f + b*e - c*d);
			double inverse_scale = 1.0 / sqrt(norm*norm - defect*defect);
			double q[4][4];
			for (int j = 0; j < 4; ++j)
			{
				double v1 = j == 0? 1 : 0, v2 = j == 1? 1 : 0, v3 = j == 2? 1 : 0, v4 = j == 3? 1 : 0;
				double t1 = s*v1 + a*v2 + b*v3 + c*v4;
				double t2 = s*v2 - a*v1 + d*v3 + e*v4;
				double t3 = s*v3 - b*v1 - d*v2 + f*v4;
				double t4 = s*v4 - c*v1 - e*v2 - f*v3;
				double t123 = d*v1 - b*v2 + a*v3 + p*v4;
				double t124 = e*v1 - c*v2 + a*v4 - p*v3;
				double t134 = f*v1 + p*v2 - c*v3 + b*v4;
				double t234 = -p*v1 + f*v2 - e*v3 + d*v4;
				q[0][j] = (t1*s + t2*a + t3*b + t4*c + t123*d + t124*e + t134*f + t234*p) * inverse_scale;
				q[1][j] = (t2*s - t1*a + t3*d + t4*e - t123*b - t124*c + t234*f - t134*p) * inverse_scale;
				q[2][j] = (t3*s - t1*b - t2*d + t4*f + t123*a - t134*c - t234*e + t124*p) * inverse_scale;
				q[3][j] = (t4*s - t1*c - t2*e - t3*f + t124*a + t134*b + t234*d - t123*p) * inverse_scale;
			}
			constexpr int i[6] = { 0, 0, 0, 1, 1, 2 }, j[6] = { 1, 2, 3, 2, 3, 3 };
			for (int row = 0; row < 6; ++row)
				for (int col = 0; col < 6; ++col)
					mA[row][col] = q[i[row]][i[col]] * q[j[row]][j[col]] - q[i[row]][j[col]] * q[j[row]][i[col]];
		}

		void ToWorld(const double *inLocal, double *outWorld) const
		{
			for (int i = 0; i < 6; ++i)
			{
				outWorld[i] = 0;
				for (int j = 0; j < 6; ++j) outWorld[i] += mA[i][j] * inLocal[j];
			}
		}
		void ToLocal(const double *inWorld, double *outLocal) const
		{
			for (int i = 0; i < 6; ++i)
			{
				outLocal[i] = 0;
				for (int j = 0; j < 6; ++j) outLocal[i] += mA[j][i] * inWorld[j];
			}
		}

	private:
		double mA[6][6];
	};
}

Rotor MotionProperties::IntegrateGyroscopicRotationInternal(RotorArg inBodyRotation, float inDeltaTime)
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite));
	JPH_ASSERT(mCachedBodyType == EBodyType::RigidBody);
	JPH_ASSERT(mCachedMotionType == EMotionType::Dynamic);

	// Free rotation requires all six rotational DOFs and a nonsingular inertia.
	// World-space locks model external constraint torques, not an isolated top.
	if ((mAllowedDOFs & EAllowedDOFs::AllRotation) != EAllowedDOFs::AllRotation)
		return (Rotor::sExp(0.5f * inDeltaTime * mAngularVelocity) * inBodyRotation).Normalized();
	for (int k = 0; k < 6; ++k)
		if (mInvInertiaDiagonal[k] <= 0.0f)
			return (Rotor::sExp(0.5f * inDeltaTime * mAngularVelocity) * inBodyRotation).Normalized();

	MotionPropertiesImpl::GyroscopicAdjoint initial(inBodyRotation * mInertiaRotation);
	double velocity[6], momentum[6], world_momentum[6];
	for (int k = 0; k < 6; ++k) velocity[k] = mAngularVelocity[k];
	initial.ToLocal(velocity, momentum);
	for (int k = 0; k < 6; ++k) momentum[k] /= mInvInertiaDiagonal[k];
	initial.ToWorld(momentum, world_momentum);
	Bivec local_momentum { float(momentum[0]), float(momentum[1]), float(momentum[2]),
		float(momentum[3]), float(momentum[4]), float(momentum[5]) };

	// H(L) = sum_k L_k^2 / (2 I_k). Each single-plane Hamiltonian has an
	// EXACT flow: Q <- Q exp(e_k L_k dt / (2 I_k)),
	// L <- Ad(exp(-e_k L_k dt / (2 I_k))) L. L_k is constant in that flow.
	// Compose these six flows symmetrically (Strang, 11 plane rotations).
	// This is a second-order, time-reversible Lie-Poisson splitting: world
	// momentum and the coadjoint orbit are preserved, energy has bounded
	// splitting error, NOT an artificial momentum/energy length projection.
	// The old forward-Euler commutator had the wrong sign as well as secular
	// error: in this convention local Ldot = -[Omega_matrix, L_matrix].
	Rotor delta = Rotor::sIdentity();
	for (int stage = 0; stage < 11; ++stage)
	{
		int k = stage <= 5? stage : 10 - stage;
		float h = (stage == 5? inDeltaTime : 0.5f * inDeltaTime);
		float angle = h * mInvInertiaDiagonal[k] * local_momentum[k];
		float plane[6] = { 0, 0, 0, 0, 0, 0 };
		plane[k] = 0.5f * angle;
		Rotor step = Rotor::sExp(Bivec(plane[0], plane[1], plane[2], plane[3], plane[4], plane[5]));
		delta = delta * step;
		local_momentum = step.InverseRotateBivec(local_momentum);
	}

	// Derive velocity from THE SAME normalized endpoint orientation used by
	// the body. Momentum is a step-local derived value, never separate state.
	// Reconstructing from world L avoids accumulating eleven adjoint rounding
	// errors per step. No torque/impulse is applied here: those have one owner
	// in the existing force kick and contact/CCD solvers.
	Rotor rotation = (inBodyRotation * mInertiaRotation * delta * mInertiaRotation.Reversed()).Normalized();
	MotionPropertiesImpl::GyroscopicAdjoint final(rotation * mInertiaRotation);
	final.ToLocal(world_momentum, momentum);
	for (int k = 0; k < 6; ++k) momentum[k] *= mInvInertiaDiagonal[k];
	final.ToWorld(momentum, velocity);
	mAngularVelocity = Bivec(float(velocity[0]), float(velocity[1]), float(velocity[2]),
		float(velocity[3]), float(velocity[4]), float(velocity[5]));
	return rotation;
}

void MotionProperties::ApplyForceTorqueAndDragInternal(RotorArg inBodyRotation, Vec4Arg inGravity, float inDeltaTime)
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite));
	JPH_ASSERT(mCachedBodyType == EBodyType::RigidBody);
	JPH_ASSERT(mCachedMotionType == EMotionType::Dynamic);

	// Update linear velocity
	mLinearVelocity = LockTranslation(mLinearVelocity + inDeltaTime * (mGravityFactor * inGravity + mInvMass * GetAccumulatedForce()));

	// Update angular velocity
	mAngularVelocity += inDeltaTime * MultiplyWorldSpaceInverseInertiaByVector(inBodyRotation, GetAccumulatedTorque());

	// Linear damping: dv/dt = -c * v
	// Solution: v(t) = v(0) * e^(-c * t) or v2 = v1 * e^(-c * dt)
	// Taylor expansion of e^(-c * dt) = 1 - c * dt + ...
	mLinearVelocity *= max(0.0f, 1.0f - mLinearDamping * inDeltaTime);
	mAngularVelocity *= max(0.0f, 1.0f - mAngularDamping * inDeltaTime);

	// Clamp velocities
	ClampLinearVelocity();
	ClampAngularVelocity();
}

void MotionProperties::ResetSleepTestSpheres(const RVec4 *inPoints)
{
#ifdef JPH_DOUBLE_PRECISION
	// Make spheres relative to the first point and initialize them to zero radius
	DVec4 offset = inPoints[0];
	offset.StoreDouble4(&mSleepTestOffset);
	mSleepTestSpheres[0] = Sphere(Vec4::sZero(), 0.0f);
	for (int i = 1; i < 3; ++i)
	{
		DVec4 diff = inPoints[i] - offset;
		mSleepTestSpheres[i] = Sphere(Vec4(float(diff.GetX()), float(diff.GetY()), float(diff.GetZ()), float(diff.GetW())), 0.0f);
	}
#else
	// Initialize the spheres to zero radius around the supplied points
	for (int i = 0; i < 3; ++i)
		mSleepTestSpheres[i] = Sphere(inPoints[i], 0.0f);
#endif

	mSleepTestTimer = 0.0f;
}

ECanSleep MotionProperties::AccumulateSleepTime(float inDeltaTime, float inTimeBeforeSleep)
{
	mSleepTestTimer += inDeltaTime;
	return mSleepTestTimer >= inTimeBeforeSleep? ECanSleep::CanSleep : ECanSleep::CannotSleep;
}

JPH_NAMESPACE_END
