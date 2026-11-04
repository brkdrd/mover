// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Physics/Body/Body.h>
#include <Jolt/Physics/Constraints/ConstraintPart/SpringPart.h>
#include <Jolt/Physics/Constraints/SpringSettings.h>
#include <Jolt/Physics/StateRecorder.h>
#include <Jolt/Physics/DeterminismLog.h>

JPH_NAMESPACE_BEGIN

/// Constraint that constrains motion along 1 axis (4D version).
///
/// @see "Constraints Derivation for Rigid Body Simulation in 3D" - Daniel Chappuis, section 2.1.1
///
/// Constraint equation (eq 25): C = (p2 - p1) . n
///
/// In 4D the angular part of the Jacobian is a bivector: where 3D used the cross product r x n, we use
/// the wedge product n ^ r -- note the swapped operands. The point velocity operator Bivec::Contract is
/// the *negative* adjoint of the wedge, so n ^ r (not r ^ n) is the bivector satisfying
/// <n ^ r, w> == n . (w _| r), which is what a Jacobian row has to compute. The angular velocity w and
/// the angular Jacobian are bivectors, and the world-space inverse inertia is the 6x6 operator applied
/// via MotionProperties::MultiplyWorldSpaceInverseInertiaByVector. For a single axis the effective mass
/// is still a scalar: K = J M^-1 J^T = invM1 + invM2 + (n^r1).(I1^-1 (n^r1)) + (n^r2).(I2^-1 (n^r2)).
///
/// n = constraint axis (normalized), p1/p2 = constraint points, r1 = p1 - x1, r2 = p2 - x2,
/// u = p2 - p1, x1/x2 = center of mass. v = [v1, w1, v2, w2] with w1/w2 bivector angular velocities.
template <bool StoreAngularJacobians>
class AxisConstraintPartT
{
	/// Solver rows are streamed on every iteration. Store only the six physical components,
	/// reconstructing the two zero SIMD lanes on load. Arithmetic still uses Bivec/Vec8.
	class StoredBivec
	{
	public:
		JPH_INLINE StoredBivec &operator = (const Bivec &inValue)
		{
			for (uint i = 0; i < 6; ++i)
				mComponents[i] = inValue[i];
			return *this;
		}

		JPH_INLINE Bivec Load() const
		{
			return Bivec(mComponents[0], mComponents[1], mComponents[2],
				mComponents[3], mComponents[4], mComponents[5]);
		}

	private:
		float mComponents[6];
	};
	static_assert(sizeof(StoredBivec) == 6 * sizeof(float));

	/// Build the world-space inverse inertia operator for a body, or a zero operator if it has none.
	/// TODO(4D): the operator is built once per body here and reused for every row; see
	/// InverseInertiaOperator for why that matters far more in 4D than it does in 3D.
	static JPH_INLINE InverseInertiaOperator sGetInvInertia(const Body &inBody, float inInvInertiaScale)
	{
		return inBody.IsDynamic()? inBody.GetMotionPropertiesUnchecked()->GetInverseInertiaForRotation(inBody.GetRotation(), inInvInertiaScale) : InverseInertiaOperator::sZero();
	}

	/// Internal helper function to update velocities of bodies after Lagrange multiplier is calculated
	template <EMotionType Type1, EMotionType Type2, class VelocityState>
	JPH_INLINE bool				ApplyVelocityStep(VelocityState *ioMotionProperties1, float inInvMass1, VelocityState *ioMotionProperties2, float inInvMass2, Vec4Arg inWorldSpaceAxis, float inLambda) const
	{
		// Apply impulse if delta is not zero
		if (inLambda != 0.0f)
		{
			// Impulse P = J^T lambda, Euler velocity integration v' = v + M^-1 P
			if constexpr (Type1 == EMotionType::Dynamic)
			{
				ioMotionProperties1->SubLinearVelocityStep((inLambda * inInvMass1) * inWorldSpaceAxis);
				ioMotionProperties1->SubAngularVelocityStep(inLambda * mInvI1_R1PlusUxAxis.Load());
			}
			if constexpr (Type2 == EMotionType::Dynamic)
			{
				ioMotionProperties2->AddLinearVelocityStep((inLambda * inInvMass2) * inWorldSpaceAxis);
				ioMotionProperties2->AddAngularVelocityStep(inLambda * mInvI2_R2xAxis.Load());
			}
			return true;
		}

		return false;
	}

	/// Internal helper function to calculate the inverse effective mass
	template <EMotionType Type1, EMotionType Type2>
	JPH_INLINE float			TemplatedCalculateInverseEffectiveMass(float inInvMass1, const InverseInertiaOperator &inInvI1, Vec4Arg inR1PlusU, float inInvMass2, const InverseInertiaOperator &inInvI2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis)
	{
		JPH_ASSERT(inWorldSpaceAxis.IsNormalized(1.0e-5f));

		// Angular Jacobian bivectors n ^ r (the 4D analog of the 3D r x n has its operands swapped,
		// because Bivec::Contract is the negative adjoint of the wedge -- this is what makes
		// <n ^ r, w> equal n . (w _| r), the constrained component of the point velocity)
		Bivec r1_axis, r2_axis;
		if constexpr (Type1 != EMotionType::Static)
		{
			r1_axis = Bivec::sWedge(inWorldSpaceAxis, inR1PlusU);
			if constexpr (StoreAngularJacobians) mR1PlusUxAxis = r1_axis;
		}
		if constexpr (Type2 != EMotionType::Static)
		{
			r2_axis = Bivec::sWedge(inWorldSpaceAxis, inR2);
			if constexpr (StoreAngularJacobians) mR2xAxis = r2_axis;
		}

		// Calculate inverse effective mass: K = J M^-1 J^T
		float inv_effective_mass;

		if constexpr (Type1 == EMotionType::Dynamic)
		{
			mInvI1_R1PlusUxAxis = inInvI1.Multiply(r1_axis);
			inv_effective_mass = inInvMass1 + mInvI1_R1PlusUxAxis.Load().Dot(r1_axis);
		}
		else
		{
			inv_effective_mass = 0.0f;
		}

		if constexpr (Type2 == EMotionType::Dynamic)
		{
			mInvI2_R2xAxis = inInvI2.Multiply(r2_axis);
			inv_effective_mass += inInvMass2 + mInvI2_R2xAxis.Load().Dot(r2_axis);
		}

		return inv_effective_mass;
	}

	/// Internal helper function to calculate the inverse effective mass
	JPH_INLINE float			CalculateInverseEffectiveMass(const Body &inBody1, Vec4Arg inR1PlusU, const Body &inBody2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis)
	{
		return CalculateInverseEffectiveMassWithMassOverride(inBody1, inBody1.IsDynamic()? inBody1.GetMotionPropertiesUnchecked()->GetInverseMass() : 0.0f, 1.0f, inR1PlusU, inBody2, inBody2.IsDynamic()? inBody2.GetMotionPropertiesUnchecked()->GetInverseMass() : 0.0f, 1.0f, inR2, inWorldSpaceAxis);
	}

	/// Internal helper function to calculate the inverse effective mass, version that supports mass scaling
	JPH_INLINE float			CalculateInverseEffectiveMassWithMassOverride(const Body &inBody1, float inInvMass1, float inInvInertiaScale1, Vec4Arg inR1PlusU, const Body &inBody2, float inInvMass2, float inInvInertiaScale2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis)
	{
		InverseInertiaOperator inv_i1 = sGetInvInertia(inBody1, inInvInertiaScale1);
		InverseInertiaOperator inv_i2 = sGetInvInertia(inBody2, inInvInertiaScale2);
		InverseInertiaOperator zero = InverseInertiaOperator::sZero();

		// Dispatch to the correct templated form
		switch (inBody1.GetMotionType())
		{
		case EMotionType::Dynamic:
			switch (inBody2.GetMotionType())
			{
			case EMotionType::Dynamic:
				return TemplatedCalculateInverseEffectiveMass<EMotionType::Dynamic, EMotionType::Dynamic>(inInvMass1, inv_i1, inR1PlusU, inInvMass2, inv_i2, inR2, inWorldSpaceAxis);

			case EMotionType::Kinematic:
				return TemplatedCalculateInverseEffectiveMass<EMotionType::Dynamic, EMotionType::Kinematic>(inInvMass1, inv_i1, inR1PlusU, 0.0f, zero, inR2, inWorldSpaceAxis);

			case EMotionType::Static:
				return TemplatedCalculateInverseEffectiveMass<EMotionType::Dynamic, EMotionType::Static>(inInvMass1, inv_i1, inR1PlusU, 0.0f, zero, inR2, inWorldSpaceAxis);

			default:
				break;
			}
			break;

		case EMotionType::Kinematic:
			JPH_ASSERT(inBody2.IsDynamic());
			return TemplatedCalculateInverseEffectiveMass<EMotionType::Kinematic, EMotionType::Dynamic>(0.0f, zero, inR1PlusU, inInvMass2, inv_i2, inR2, inWorldSpaceAxis);

		case EMotionType::Static:
			JPH_ASSERT(inBody2.IsDynamic());
			return TemplatedCalculateInverseEffectiveMass<EMotionType::Static, EMotionType::Dynamic>(0.0f, zero, inR1PlusU, inInvMass2, inv_i2, inR2, inWorldSpaceAxis);

		default:
			break;
		}

		JPH_ASSERT(false);
		return 0.0f;
	}

public:
	/// Templated form of CalculateConstraintProperties with the motion types baked in
	template <EMotionType Type1, EMotionType Type2>
	JPH_INLINE void				TemplatedCalculateConstraintProperties(float inInvMass1, const InverseInertiaOperator &inInvI1, Vec4Arg inR1PlusU, float inInvMass2, const InverseInertiaOperator &inInvI2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis, float inBias = 0.0f)
	{
		float inv_effective_mass = TemplatedCalculateInverseEffectiveMass<Type1, Type2>(inInvMass1, inInvI1, inR1PlusU, inInvMass2, inInvI2, inR2, inWorldSpaceAxis);

		if (inv_effective_mass == 0.0f)
			Deactivate();
		else
		{
			mEffectiveMass = 1.0f / inv_effective_mass;
			mSpringPart.CalculateSpringPropertiesWithBias(inBias);
		}

		JPH_DET_LOG("TemplatedCalculateConstraintProperties: invM1: " << inInvMass1 << " invM2: " << inInvMass2 << " bias: " << inBias << " r1PlusUxAxis: " << Bivec::sWedge(inWorldSpaceAxis, inR1PlusU) << " r2xAxis: " << Bivec::sWedge(inWorldSpaceAxis, inR2) << " invI1_R1PlusUxAxis: " << mInvI1_R1PlusUxAxis.Load() << " invI2_R2xAxis: " << mInvI2_R2xAxis.Load() << " effectiveMass: " << mEffectiveMass << " totalLambda: " << mTotalLambda);
	}

	/// Calculate properties used during the functions below
	/// @param inBody1 The first body that this constraint is attached to
	/// @param inBody2 The second body that this constraint is attached to
	/// @param inR1PlusU See equations above (r1 + u)
	/// @param inR2 See equations above (r2)
	/// @param inWorldSpaceAxis Axis along which the constraint acts (normalized, pointing from body 1 to 2)
	/// @param inBias Bias term (b) for the constraint impulse: lambda = J v + b
	inline void					CalculateConstraintProperties(const Body &inBody1, Vec4Arg inR1PlusU, const Body &inBody2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis, float inBias = 0.0f)
	{
		float inv_effective_mass = CalculateInverseEffectiveMass(inBody1, inR1PlusU, inBody2, inR2, inWorldSpaceAxis);

		if (inv_effective_mass == 0.0f)
			Deactivate();
		else
		{
			mEffectiveMass = 1.0f / inv_effective_mass;
			mSpringPart.CalculateSpringPropertiesWithBias(inBias);
		}
	}

	/// Calculate properties used during the functions below, version that supports mass scaling
	inline void					CalculateConstraintPropertiesWithMassOverride(const Body &inBody1, float inInvMass1, float inInvInertiaScale1, Vec4Arg inR1PlusU, const Body &inBody2, float inInvMass2, float inInvInertiaScale2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis, float inBias = 0.0f)
	{
		float inv_effective_mass = CalculateInverseEffectiveMassWithMassOverride(inBody1, inInvMass1, inInvInertiaScale1, inR1PlusU, inBody2, inInvMass2, inInvInertiaScale2, inR2, inWorldSpaceAxis);

		if (inv_effective_mass == 0.0f)
			Deactivate();
		else
		{
			mEffectiveMass = 1.0f / inv_effective_mass;
			mSpringPart.CalculateSpringPropertiesWithBias(inBias);
		}
	}

	/// Calculate properties used during the functions below
	inline void					CalculateConstraintPropertiesWithFrequencyAndDamping(float inDeltaTime, const Body &inBody1, Vec4Arg inR1PlusU, const Body &inBody2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis, float inBias, float inC, float inFrequency, float inDamping)
	{
		float inv_effective_mass = CalculateInverseEffectiveMass(inBody1, inR1PlusU, inBody2, inR2, inWorldSpaceAxis);

		if (inv_effective_mass == 0.0f)
			Deactivate();
		else
			mSpringPart.CalculateSpringPropertiesWithFrequencyAndDamping(inDeltaTime, inv_effective_mass, inBias, inC, inFrequency, inDamping, mEffectiveMass);
	}

	/// Calculate properties used during the functions below
	inline void					CalculateConstraintPropertiesWithStiffnessAndDamping(float inDeltaTime, const Body &inBody1, Vec4Arg inR1PlusU, const Body &inBody2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis, float inBias, float inC, float inStiffness, float inDamping)
	{
		float inv_effective_mass = CalculateInverseEffectiveMass(inBody1, inR1PlusU, inBody2, inR2, inWorldSpaceAxis);

		if (inv_effective_mass == 0.0f)
			Deactivate();
		else
			mSpringPart.CalculateSpringPropertiesWithStiffnessAndDamping(inDeltaTime, inv_effective_mass, inBias, inC, inStiffness, inDamping, mEffectiveMass);
	}

	/// Selects one of the above functions based on the spring settings
	inline void					CalculateConstraintPropertiesWithSettings(float inDeltaTime, const Body &inBody1, Vec4Arg inR1PlusU, const Body &inBody2, Vec4Arg inR2, Vec4Arg inWorldSpaceAxis, float inBias, float inC, const SpringSettings &inSpringSettings)
	{
		float inv_effective_mass = CalculateInverseEffectiveMass(inBody1, inR1PlusU, inBody2, inR2, inWorldSpaceAxis);

		if (inv_effective_mass == 0.0f)
			Deactivate();
		else if (inSpringSettings.mMode == ESpringMode::FrequencyAndDamping)
			mSpringPart.CalculateSpringPropertiesWithFrequencyAndDamping(inDeltaTime, inv_effective_mass, inBias, inC, inSpringSettings.mFrequency, inSpringSettings.mDamping, mEffectiveMass);
		else
			mSpringPart.CalculateSpringPropertiesWithStiffnessAndDamping(inDeltaTime, inv_effective_mass, inBias, inC, inSpringSettings.mStiffness, inSpringSettings.mDamping, mEffectiveMass);
	}

	/// Deactivate this constraint
	inline void					Deactivate()
	{
		mEffectiveMass = 0.0f;
		mTotalLambda = 0.0f;
	}

	/// Check if constraint is active
	inline bool					IsActive() const
	{
		return mEffectiveMass != 0.0f;
	}

	/// Templated form of WarmStart with the motion types baked in
	template <EMotionType Type1, EMotionType Type2>
	inline void					TemplatedWarmStart(MotionProperties *ioMotionProperties1, float inInvMass1, MotionProperties *ioMotionProperties2, float inInvMass2, Vec4Arg inWorldSpaceAxis, float inWarmStartImpulseRatio)
	{
		mTotalLambda *= inWarmStartImpulseRatio;

		ApplyVelocityStep<Type1, Type2>(ioMotionProperties1, inInvMass1, ioMotionProperties2, inInvMass2, inWorldSpaceAxis, mTotalLambda);
	}

	/// Must be called from the WarmStartVelocityConstraint call to apply the previous frame's impulses
	inline void					WarmStart(Body &ioBody1, Body &ioBody2, Vec4Arg inWorldSpaceAxis, float inWarmStartImpulseRatio)
	{
		EMotionType motion_type1 = ioBody1.GetMotionType();
		MotionProperties *motion_properties1 = ioBody1.GetMotionPropertiesUnchecked();

		EMotionType motion_type2 = ioBody2.GetMotionType();
		MotionProperties *motion_properties2 = ioBody2.GetMotionPropertiesUnchecked();

		// Note: Warm starting doesn't differentiate between kinematic/static bodies so we handle both as static bodies
		if (motion_type1 == EMotionType::Dynamic)
		{
			if (motion_type2 == EMotionType::Dynamic)
				TemplatedWarmStart<EMotionType::Dynamic, EMotionType::Dynamic>(motion_properties1, motion_properties1->GetInverseMass(), motion_properties2, motion_properties2->GetInverseMass(), inWorldSpaceAxis, inWarmStartImpulseRatio);
			else
				TemplatedWarmStart<EMotionType::Dynamic, EMotionType::Static>(motion_properties1, motion_properties1->GetInverseMass(), motion_properties2, 0.0f /* Unused */, inWorldSpaceAxis, inWarmStartImpulseRatio);
		}
		else
		{
			JPH_ASSERT(motion_type2 == EMotionType::Dynamic);
			TemplatedWarmStart<EMotionType::Static, EMotionType::Dynamic>(motion_properties1, 0.0f /* Unused */, motion_properties2, motion_properties2->GetInverseMass(), inWorldSpaceAxis, inWarmStartImpulseRatio);
		}
	}

	/// Templated form of SolveVelocityConstraint with the motion types baked in, part 1: get the total lambda
	template <EMotionType Type1, EMotionType Type2>
	JPH_INLINE float			TemplatedSolveVelocityConstraintGetTotalLambda(const MotionProperties *ioMotionProperties1, const MotionProperties *ioMotionProperties2, Vec4Arg inWorldSpaceAxis) const
	{
		static_assert(StoreAngularJacobians, "Contact rows use the shared point velocity instead");
		// Calculate jacobian multiplied by linear velocity
		float jv;
		if constexpr (Type1 != EMotionType::Static && Type2 != EMotionType::Static)
			jv = inWorldSpaceAxis.Dot(ioMotionProperties1->GetLinearVelocity() - ioMotionProperties2->GetLinearVelocity());
		else if constexpr (Type1 != EMotionType::Static)
			jv = inWorldSpaceAxis.Dot(ioMotionProperties1->GetLinearVelocity());
		else if constexpr (Type2 != EMotionType::Static)
			jv = inWorldSpaceAxis.Dot(-ioMotionProperties2->GetLinearVelocity());
		else
			JPH_ASSERT(false); // Static vs static is nonsensical!

		// Calculate jacobian multiplied by angular velocity (bivector inner products)
		if constexpr (Type1 != EMotionType::Static)
			jv += mR1PlusUxAxis.Load().Dot(ioMotionProperties1->GetAngularVelocity());
		if constexpr (Type2 != EMotionType::Static)
			jv -= mR2xAxis.Load().Dot(ioMotionProperties2->GetAngularVelocity());

		return GetTotalLambdaFromVelocity(jv);
	}

	/// Contact rows share a point velocity across all tangents. By the wedge/contraction
	/// adjoint identity, axis . relative_point_velocity is the same Jacobian product.
	JPH_INLINE float GetTotalLambdaFromVelocity(float inJacobianVelocity) const
	{
		return mTotalLambda + mEffectiveMass * (inJacobianVelocity - mSpringPart.GetBias(mTotalLambda));
	}

	/// Templated form of SolveVelocityConstraint with the motion types baked in, part 2: apply new lambda
	template <EMotionType Type1, EMotionType Type2, class VelocityState>
	JPH_INLINE bool				TemplatedSolveVelocityConstraintApplyLambda(VelocityState *ioMotionProperties1, float inInvMass1, VelocityState *ioMotionProperties2, float inInvMass2, Vec4Arg inWorldSpaceAxis, float inTotalLambda)
	{
		float delta_lambda = inTotalLambda - mTotalLambda; // Calculate change in lambda
		mTotalLambda = inTotalLambda; // Store accumulated impulse

		return ApplyVelocityStep<Type1, Type2>(ioMotionProperties1, inInvMass1, ioMotionProperties2, inInvMass2, inWorldSpaceAxis, delta_lambda);
	}

	/// Templated form of SolveVelocityConstraint with the motion types baked in
	template <EMotionType Type1, EMotionType Type2>
	inline bool					TemplatedSolveVelocityConstraint(MotionProperties *ioMotionProperties1, float inInvMass1, MotionProperties *ioMotionProperties2, float inInvMass2, Vec4Arg inWorldSpaceAxis, float inMinLambda, float inMaxLambda)
	{
		float total_lambda = TemplatedSolveVelocityConstraintGetTotalLambda<Type1, Type2>(ioMotionProperties1, ioMotionProperties2, inWorldSpaceAxis);

		// Clamp impulse to specified range
		total_lambda = Clamp(total_lambda, inMinLambda, inMaxLambda);

		return TemplatedSolveVelocityConstraintApplyLambda<Type1, Type2>(ioMotionProperties1, inInvMass1, ioMotionProperties2, inInvMass2, inWorldSpaceAxis, total_lambda);
	}

	/// Iteratively update the velocity constraint. Makes sure d/dt C(...) = 0, where C is the constraint equation.
	inline bool					SolveVelocityConstraint(Body &ioBody1, Body &ioBody2, Vec4Arg inWorldSpaceAxis, float inMinLambda, float inMaxLambda)
	{
		EMotionType motion_type1 = ioBody1.GetMotionType();
		MotionProperties *motion_properties1 = ioBody1.GetMotionPropertiesUnchecked();

		EMotionType motion_type2 = ioBody2.GetMotionType();
		MotionProperties *motion_properties2 = ioBody2.GetMotionPropertiesUnchecked();

		// Dispatch to the correct templated form
		switch (motion_type1)
		{
		case EMotionType::Dynamic:
			switch (motion_type2)
			{
			case EMotionType::Dynamic:
				return TemplatedSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Dynamic>(motion_properties1, motion_properties1->GetInverseMass(), motion_properties2, motion_properties2->GetInverseMass(), inWorldSpaceAxis, inMinLambda, inMaxLambda);

			case EMotionType::Kinematic:
				return TemplatedSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Kinematic>(motion_properties1, motion_properties1->GetInverseMass(), motion_properties2, 0.0f /* Unused */, inWorldSpaceAxis, inMinLambda, inMaxLambda);

			case EMotionType::Static:
				return TemplatedSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Static>(motion_properties1, motion_properties1->GetInverseMass(), motion_properties2, 0.0f /* Unused */, inWorldSpaceAxis, inMinLambda, inMaxLambda);

			default:
				JPH_ASSERT(false);
				break;
			}
			break;

		case EMotionType::Kinematic:
			JPH_ASSERT(motion_type2 == EMotionType::Dynamic);
			return TemplatedSolveVelocityConstraint<EMotionType::Kinematic, EMotionType::Dynamic>(motion_properties1, 0.0f /* Unused */, motion_properties2, motion_properties2->GetInverseMass(), inWorldSpaceAxis, inMinLambda, inMaxLambda);

		case EMotionType::Static:
			JPH_ASSERT(motion_type2 == EMotionType::Dynamic);
			return TemplatedSolveVelocityConstraint<EMotionType::Static, EMotionType::Dynamic>(motion_properties1, 0.0f /* Unused */, motion_properties2, motion_properties2->GetInverseMass(), inWorldSpaceAxis, inMinLambda, inMaxLambda);

		default:
			JPH_ASSERT(false);
			break;
		}

		return false;
	}

	/// Iteratively update the velocity constraint, version that supports mass scaling.
	inline bool					SolveVelocityConstraintWithMassOverride(Body &ioBody1, float inInvMass1, Body &ioBody2, float inInvMass2, Vec4Arg inWorldSpaceAxis, float inMinLambda, float inMaxLambda)
	{
		EMotionType motion_type1 = ioBody1.GetMotionType();
		MotionProperties *motion_properties1 = ioBody1.GetMotionPropertiesUnchecked();

		EMotionType motion_type2 = ioBody2.GetMotionType();
		MotionProperties *motion_properties2 = ioBody2.GetMotionPropertiesUnchecked();

		// Dispatch to the correct templated form
		switch (motion_type1)
		{
		case EMotionType::Dynamic:
			switch (motion_type2)
			{
			case EMotionType::Dynamic:
				return TemplatedSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Dynamic>(motion_properties1, inInvMass1, motion_properties2, inInvMass2, inWorldSpaceAxis, inMinLambda, inMaxLambda);

			case EMotionType::Kinematic:
				return TemplatedSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Kinematic>(motion_properties1, inInvMass1, motion_properties2, 0.0f /* Unused */, inWorldSpaceAxis, inMinLambda, inMaxLambda);

			case EMotionType::Static:
				return TemplatedSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Static>(motion_properties1, inInvMass1, motion_properties2, 0.0f /* Unused */, inWorldSpaceAxis, inMinLambda, inMaxLambda);

			default:
				JPH_ASSERT(false);
				break;
			}
			break;

		case EMotionType::Kinematic:
			JPH_ASSERT(motion_type2 == EMotionType::Dynamic);
			return TemplatedSolveVelocityConstraint<EMotionType::Kinematic, EMotionType::Dynamic>(motion_properties1, 0.0f /* Unused */, motion_properties2, inInvMass2, inWorldSpaceAxis, inMinLambda, inMaxLambda);

		case EMotionType::Static:
			JPH_ASSERT(motion_type2 == EMotionType::Dynamic);
			return TemplatedSolveVelocityConstraint<EMotionType::Static, EMotionType::Dynamic>(motion_properties1, 0.0f /* Unused */, motion_properties2, inInvMass2, inWorldSpaceAxis, inMinLambda, inMaxLambda);

		default:
			JPH_ASSERT(false);
			break;
		}

		return false;
	}

	/// Iteratively update the position constraint. Makes sure C(...) = 0.
	inline bool					SolvePositionConstraint(Body &ioBody1, Body &ioBody2, Vec4Arg inWorldSpaceAxis, float inC, float inBaumgarte) const
	{
		// Only apply position constraint when the constraint is hard, otherwise the velocity bias will fix the constraint
		if (inC != 0.0f && !mSpringPart.IsActive())
		{
			// lambda = -K^-1 * beta / dt * C (the dt cancels with the position integrate below)
			float lambda = -mEffectiveMass * inBaumgarte * inC;

			// Directly integrate the position change for one time step, discarding the velocity change (Catto GDC 2007)
			if (ioBody1.IsDynamic())
			{
				ioBody1.SubPositionStep((lambda * ioBody1.GetMotionProperties()->GetInverseMass()) * inWorldSpaceAxis);
				ioBody1.SubRotationStep(lambda * mInvI1_R1PlusUxAxis.Load());
			}
			if (ioBody2.IsDynamic())
			{
				ioBody2.AddPositionStep((lambda * ioBody2.GetMotionProperties()->GetInverseMass()) * inWorldSpaceAxis);
				ioBody2.AddRotationStep(lambda * mInvI2_R2xAxis.Load());
			}
			return true;
		}

		return false;
	}

	/// Iteratively update the position constraint, version that supports mass scaling.
	inline bool					SolvePositionConstraintWithMassOverride(Body &ioBody1, float inInvMass1, Body &ioBody2, float inInvMass2, Vec4Arg inWorldSpaceAxis, float inC, float inBaumgarte) const
	{
		if (inC != 0.0f && !mSpringPart.IsActive())
		{
			float lambda = -mEffectiveMass * inBaumgarte * inC;

			if (ioBody1.IsDynamic())
			{
				ioBody1.SubPositionStep((lambda * inInvMass1) * inWorldSpaceAxis);
				ioBody1.SubRotationStep(lambda * mInvI1_R1PlusUxAxis.Load());
			}
			if (ioBody2.IsDynamic())
			{
				ioBody2.AddPositionStep((lambda * inInvMass2) * inWorldSpaceAxis);
				ioBody2.AddRotationStep(lambda * mInvI2_R2xAxis.Load());
			}
			return true;
		}

		return false;
	}

	/// Override total lagrange multiplier, can be used to set the initial value for warm starting
	inline void					SetTotalLambda(float inLambda)
	{
		mTotalLambda = inLambda;
	}

	/// Return lagrange multiplier
	inline float				GetTotalLambda() const
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
	struct NoJacobian { };
	// Contact rows have restitution / surface-velocity bias, but no spring softness.
	// Removing the unused coefficient also makes a compact contact row one cache line.
	class ContactBias
	{
	public:
		void CalculateSpringPropertiesWithBias(float inBias) { mBias = inBias; }
		float GetBias(float) const { return mBias; }
		bool IsActive() const { return false; }
	private:
		float mBias = 0.0f;
	};
	using AngularJacobian = std::conditional_t<StoreAngularJacobians, StoredBivec, NoJacobian>;
	AngularJacobian				mR1PlusUxAxis;			///< Angular Jacobian bivector (r1 + u) ^ n
	AngularJacobian				mR2xAxis;				///< Angular Jacobian bivector r2 ^ n
	StoredBivec					mInvI1_R1PlusUxAxis;	///< I1^-1 applied to the body 1 angular Jacobian
	StoredBivec					mInvI2_R2xAxis;			///< I2^-1 applied to the body 2 angular Jacobian
	float						mEffectiveMass = 0.0f;
	std::conditional_t<StoreAngularJacobians, SpringPart, ContactBias> mSpringPart;
	float						mTotalLambda = 0.0f;
};

using AxisConstraintPart = AxisConstraintPartT<true>;
using ContactConstraintPart = AxisConstraintPartT<false>;
static_assert(sizeof(ContactConstraintPart) < sizeof(AxisConstraintPart));
static_assert(sizeof(ContactConstraintPart) == 64);

JPH_NAMESPACE_END
