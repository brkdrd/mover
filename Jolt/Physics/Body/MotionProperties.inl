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

void MotionProperties::ApplyGyroscopicForceInternal(RotorArg inBodyRotation, float inDeltaTime)
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite));
	JPH_ASSERT(mCachedBodyType == EBodyType::RigidBody);
	JPH_ASSERT(mCachedMotionType == EMotionType::Dynamic);

	// Calculate local space inertia (reciprocal of inverse, avoiding division by zero)
	Bivec inv_diag = mInvInertiaDiagonal;
	Bivec local_inertia(
		inv_diag[0] != 0.0f ? 1.0f / inv_diag[0] : 0.0f,
		inv_diag[1] != 0.0f ? 1.0f / inv_diag[1] : 0.0f,
		inv_diag[2] != 0.0f ? 1.0f / inv_diag[2] : 0.0f,
		inv_diag[3] != 0.0f ? 1.0f / inv_diag[3] : 0.0f,
		inv_diag[4] != 0.0f ? 1.0f / inv_diag[4] : 0.0f,
		inv_diag[5] != 0.0f ? 1.0f / inv_diag[5] : 0.0f);

	// Calculate local space angular velocity: bring the world-space angular velocity into the body's
	// inertia principal frame via the inverse adjoint (Ad(R^{-1}), R = body_rotation * inertia_rotation).
	Rotor combined = inBodyRotation * mInertiaRotation;
	Bivec local_angular_velocity = combined.InverseRotateBivec(mAngularVelocity);
	Bivec local_momentum = local_inertia * local_angular_velocity;

	// The gyroscopic torque in 4D: T = [Omega, L] (commutator of bivectors)
	// [Omega, L]_ij = sum_k (Omega_ik * L_kj - L_ik * Omega_kj)
	// This is the Lie bracket of the bivectors viewed as elements of so(4).
	//
	// For now, use the simplified energy-conserving approach:
	// new_momentum = momentum - dt * [omega, momentum], then rescale to preserve |momentum|
	//
	// The commutator of two bivectors A, B in so(4):
	float w12 = local_angular_velocity.GetE12(), w13 = local_angular_velocity.GetE13(), w14 = local_angular_velocity.GetE14();
	float w23 = local_angular_velocity.GetE23(), w24 = local_angular_velocity.GetE24(), w34 = local_angular_velocity.GetE34();
	float L12 = local_momentum.GetE12(), L13 = local_momentum.GetE13(), L14 = local_momentum.GetE14();
	float L23 = local_momentum.GetE23(), L24 = local_momentum.GetE24(), L34 = local_momentum.GetE34();

	// [w, L]_ij = sum_k (w_ik * L_kj - L_ik * w_kj)
	// Antisymmetric matrix representation: w_ij = -w_ji
	Bivec commutator(
		w13*L23 + w14*L24 - w23*L13 - w24*L14,                    // [w,L]_12
		-w12*L23 + w14*L34 + w23*L12 - w34*L14,                   // [w,L]_13
		-w12*L24 - w13*L34 + w24*L12 + w34*L13,                   // [w,L]_14
		w12*L13 - w13*L12 + w24*L34 - w34*L24,                    // [w,L]_23
		w12*L14 - w14*L12 - w23*L34 + w34*L23,                    // [w,L]_24
		w13*L14 - w14*L13 + w23*L24 - w24*L23                     // [w,L]_34
	);

	Bivec new_local_momentum = local_momentum - inDeltaTime * commutator;
	float new_len_sq = new_local_momentum.LengthSq();
	new_local_momentum = new_len_sq > 0.0f ? new_local_momentum * sqrt(local_momentum.LengthSq() / new_len_sq) : Bivec::sZero();

	// Convert body-frame momentum back to a body-frame angular velocity, then rotate it to world space
	Bivec new_local_angular_velocity = mInvInertiaDiagonal * new_local_momentum;
	mAngularVelocity = combined.RotateBivec(new_local_angular_velocity);
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
	offset.StoreDouble3(&mSleepTestOffset);
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
