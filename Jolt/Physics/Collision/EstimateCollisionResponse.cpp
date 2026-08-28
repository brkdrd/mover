// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Collision/EstimateCollisionResponse.h>
#include <Jolt/Physics/Body/Body.h>

JPH_NAMESPACE_BEGIN

void EstimateCollisionResponse(const Body &inBody1, const Body &inBody2, const ContactManifold &inManifold, CollisionEstimationResult &outResult, float inCombinedFriction, float inCombinedRestitution, float inMinVelocityForRestitution, uint inNumIterations)
{
	// Note this code is based on AxisConstraintPart, see that class for more comments on the math

	ContactPoints::size_type num_points = inManifold.mRelativeContactPointsOn1.size();
	JPH_ASSERT(num_points == inManifold.mRelativeContactPointsOn2.size());

	// Start with zero impulses
	outResult.mImpulses.resize(num_points);
	memset(outResult.mImpulses.data(), 0, num_points * sizeof(CollisionEstimationResult::Impulse));

	// Calculate friction directions. TODO(4D): 3 of them, the contact tangent space is a hyperplane.
	inManifold.mWorldSpaceNormal.GetNormalizedPerpendicularBasis(outResult.mTangent1, outResult.mTangent2, outResult.mTangent3);

	// Get body velocities
	EMotionType motion_type1 = inBody1.GetMotionType();
	const MotionProperties *motion_properties1 = inBody1.GetMotionPropertiesUnchecked();
	if (motion_type1 != EMotionType::Static)
	{
		outResult.mLinearVelocity1 = motion_properties1->GetLinearVelocity();
		outResult.mAngularVelocity1 = motion_properties1->GetAngularVelocity();
	}
	else
	{
		outResult.mLinearVelocity1 = Vec4::sZero();
		outResult.mAngularVelocity1 = Bivec::sZero();
	}

	EMotionType motion_type2 = inBody2.GetMotionType();
	const MotionProperties *motion_properties2 = inBody2.GetMotionPropertiesUnchecked();
	if (motion_type2 != EMotionType::Static)
	{
		outResult.mLinearVelocity2 = motion_properties2->GetLinearVelocity();
		outResult.mAngularVelocity2 = motion_properties2->GetAngularVelocity();
	}
	else
	{
		outResult.mLinearVelocity2 = Vec4::sZero();
		outResult.mAngularVelocity2 = Bivec::sZero();
	}

	// Get inverse mass. TODO(4D): the world space inverse inertia is a 6x6 operator on bivectors, not a
	// Mat44, so it is applied through the body's rotor rather than precomputed as a matrix.
	float inv_m1 = motion_type1 == EMotionType::Dynamic? motion_properties1->GetInverseMass() : 0.0f;
	float inv_m2 = motion_type2 == EMotionType::Dynamic? motion_properties2->GetInverseMass() : 0.0f;
	const MotionProperties *inv_i1 = motion_type1 == EMotionType::Dynamic? motion_properties1 : nullptr;
	const MotionProperties *inv_i2 = motion_type2 == EMotionType::Dynamic? motion_properties2 : nullptr;
	Rotor rot1 = inBody1.GetRotation();
	Rotor rot2 = inBody2.GetRotation();

	// Get center of masses relative to the base offset
	Vec4 com1 = Vec4(inBody1.GetCenterOfMassPosition() - inManifold.mBaseOffset);
	Vec4 com2 = Vec4(inBody2.GetCenterOfMassPosition() - inManifold.mBaseOffset);

	struct AxisConstraint
	{
		inline void		Initialize(Vec4Arg inR1, Vec4Arg inR2, Vec4Arg inWorldSpaceNormal, float inInvM1, float inInvM2, RotorArg inRotation1, const MotionProperties *inInvI1, RotorArg inRotation2, const MotionProperties *inInvI2)
		{
			// Calculate effective mass: K^-1 = (J M^-1 J^T)^-1
			// TODO(4D): the angular Jacobian is the bivector n ^ r -- operands swapped w.r.t. the 3D
			// r x n, because Contract is the negative adjoint of the wedge (see AxisConstraintPart).
			mR1PlusUxAxis = Bivec::sWedge(inWorldSpaceNormal, inR1);
			mR2xAxis = Bivec::sWedge(inWorldSpaceNormal, inR2);
			mInvI1_R1PlusUxAxis = inInvI1 != nullptr? inInvI1->MultiplyWorldSpaceInverseInertiaByVector(inRotation1, mR1PlusUxAxis) : Bivec::sZero();
			mInvI2_R2xAxis = inInvI2 != nullptr? inInvI2->MultiplyWorldSpaceInverseInertiaByVector(inRotation2, mR2xAxis) : Bivec::sZero();
			mEffectiveMass = 1.0f / (inInvM1 + mInvI1_R1PlusUxAxis.Dot(mR1PlusUxAxis) + inInvM2 + mInvI2_R2xAxis.Dot(mR2xAxis));
			mBias = 0.0f;
		}

		inline float	SolveGetLambda(Vec4Arg inWorldSpaceNormal, const CollisionEstimationResult &inResult) const
		{
			// Calculate jacobian multiplied by linear/angular velocity
			float jv = inWorldSpaceNormal.Dot(inResult.mLinearVelocity1 - inResult.mLinearVelocity2) + mR1PlusUxAxis.Dot(inResult.mAngularVelocity1) - mR2xAxis.Dot(inResult.mAngularVelocity2);

			// Lagrange multiplier is:
			//
			// lambda = -K^-1 (J v + b)
			return mEffectiveMass * (jv - mBias);
		}

		inline void		SolveApplyLambda(Vec4Arg inWorldSpaceNormal, float inInvM1, float inInvM2, float inLambda, CollisionEstimationResult &ioResult) const
		{
			// Apply impulse to body velocities
			ioResult.mLinearVelocity1 -= (inLambda * inInvM1) * inWorldSpaceNormal;
			ioResult.mAngularVelocity1 -= inLambda * mInvI1_R1PlusUxAxis;
			ioResult.mLinearVelocity2 += (inLambda * inInvM2) * inWorldSpaceNormal;
			ioResult.mAngularVelocity2 += inLambda * mInvI2_R2xAxis;
		}

		inline void		Solve(Vec4Arg inWorldSpaceNormal, float inInvM1, float inInvM2, float inMinLambda, float inMaxLambda, float &ioTotalLambda, CollisionEstimationResult &ioResult) const
		{
			// Calculate new total lambda
			float total_lambda = ioTotalLambda + SolveGetLambda(inWorldSpaceNormal, ioResult);

			// Clamp impulse
			total_lambda = Clamp(total_lambda, inMinLambda, inMaxLambda);

			SolveApplyLambda(inWorldSpaceNormal, inInvM1, inInvM2, total_lambda - ioTotalLambda, ioResult);

			ioTotalLambda = total_lambda;
		}

		Bivec			mR1PlusUxAxis;
		Bivec			mR2xAxis;
		Bivec			mInvI1_R1PlusUxAxis;
		Bivec			mInvI2_R2xAxis;
		float			mEffectiveMass;
		float			mBias;
	};

	struct Constraint
	{
		AxisConstraint	mContact;
		AxisConstraint	mFriction1;
		AxisConstraint	mFriction2;
		AxisConstraint	mFriction3;
	};

	// Initialize the constraint properties
	Constraint constraints[ContactPoints::Capacity];
	for (uint c = 0; c < num_points; ++c)
	{
		Constraint &constraint = constraints[c];

		// Calculate contact points relative to body 1 and 2
		Vec4 p = 0.5f * (inManifold.mRelativeContactPointsOn1[c] + inManifold.mRelativeContactPointsOn2[c]);
		Vec4 r1 = p - com1;
		Vec4 r2 = p - com2;

		// Initialize contact constraint
		constraint.mContact.Initialize(r1, r2, inManifold.mWorldSpaceNormal, inv_m1, inv_m2, rot1, inv_i1, rot2, inv_i2);

		// Handle elastic collisions
		if (inCombinedRestitution > 0.0f)
		{
			// Calculate velocity of contact point
			Vec4 relative_velocity = outResult.mLinearVelocity2 + outResult.mAngularVelocity2.Contract(r2) - outResult.mLinearVelocity1 - outResult.mAngularVelocity1.Contract(r1);
			float normal_velocity = relative_velocity.Dot(inManifold.mWorldSpaceNormal);

			// If it is big enough, apply restitution
			if (normal_velocity < -inMinVelocityForRestitution)
				constraint.mContact.mBias = inCombinedRestitution * normal_velocity;
		}

		if (inCombinedFriction > 0.0f)
		{
			// Initialize friction constraints
			constraint.mFriction1.Initialize(r1, r2, outResult.mTangent1, inv_m1, inv_m2, rot1, inv_i1, rot2, inv_i2);
			constraint.mFriction2.Initialize(r1, r2, outResult.mTangent2, inv_m1, inv_m2, rot1, inv_i1, rot2, inv_i2);
			constraint.mFriction3.Initialize(r1, r2, outResult.mTangent3, inv_m1, inv_m2, rot1, inv_i1, rot2, inv_i2);
		}
	}

	// If there's only 1 contact point, we only need 1 iteration
	int num_iterations = inCombinedFriction <= 0.0f && num_points == 1? 1 : inNumIterations;

	// Solve iteratively
	for (int iteration = 0; iteration < num_iterations; ++iteration)
	{
		// Solve friction constraints first
		if (inCombinedFriction > 0.0f && iteration > 0) // For first iteration the contact impulse is zero so there's no point in applying friction
			for (uint c = 0; c < num_points; ++c)
			{
				const Constraint &constraint = constraints[c];
				CollisionEstimationResult::Impulse &impulse = outResult.mImpulses[c];

				float lambda1 = impulse.mFrictionImpulse1 + constraint.mFriction1.SolveGetLambda(outResult.mTangent1, outResult);
				float lambda2 = impulse.mFrictionImpulse2 + constraint.mFriction2.SolveGetLambda(outResult.mTangent2, outResult);
				float lambda3 = impulse.mFrictionImpulse3 + constraint.mFriction3.SolveGetLambda(outResult.mTangent3, outResult);

				// Calculate max impulse based on contact impulse
				float max_impulse = inCombinedFriction * impulse.mContactImpulse;

				// If the total lambda that we will apply is too large, scale it back (the friction cone
				// clamps a 3 dimensional tangential impulse in 4D)
				float total_lambda_sq = Square(lambda1) + Square(lambda2) + Square(lambda3);
				if (total_lambda_sq > Square(max_impulse))
				{
					float scale = max_impulse / sqrt(total_lambda_sq);
					lambda1 *= scale;
					lambda2 *= scale;
					lambda3 *= scale;
				}

				constraint.mFriction1.SolveApplyLambda(outResult.mTangent1, inv_m1, inv_m2, lambda1 - impulse.mFrictionImpulse1, outResult);
				constraint.mFriction2.SolveApplyLambda(outResult.mTangent2, inv_m1, inv_m2, lambda2 - impulse.mFrictionImpulse2, outResult);
				constraint.mFriction3.SolveApplyLambda(outResult.mTangent3, inv_m1, inv_m2, lambda3 - impulse.mFrictionImpulse3, outResult);

				impulse.mFrictionImpulse1 = lambda1;
				impulse.mFrictionImpulse2 = lambda2;
				impulse.mFrictionImpulse3 = lambda3;
			}

		// Solve contact constraints last
		for (uint c = 0; c < num_points; ++c)
			constraints[c].mContact.Solve(inManifold.mWorldSpaceNormal, inv_m1, inv_m2, 0.0f, FLT_MAX, outResult.mImpulses[c].mContactImpulse, outResult);
	}
}

JPH_NAMESPACE_END
