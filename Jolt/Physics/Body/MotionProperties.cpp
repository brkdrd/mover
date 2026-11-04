// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Body/MotionProperties.h>
#include <Jolt/Physics/StateRecorder.h>

JPH_NAMESPACE_BEGIN

void MotionProperties::SetMassProperties(EAllowedDOFs inAllowedDOFs, const MassProperties &inMassProperties)
{
	// Store allowed DOFs
	mAllowedDOFs = inAllowedDOFs;

	// Decompose DOFs: 4 translation bits (0-3), 6 rotation bits (4-9)
	uint allowed_translation_axis = uint(inAllowedDOFs) & 0b1111;
	uint allowed_rotation_axis = (uint(inAllowedDOFs) >> 4) & 0b111111;

	// Set inverse mass
	if (allowed_translation_axis == 0)
	{
		// No translation possible
		mInvMass = 0.0f;
	}
	else
	{
		JPH_ASSERT(inMassProperties.mMass > 0.0f, "Invalid mass. "
			"Some shapes like MeshShape or TriangleShape cannot calculate mass automatically, "
			"in this case you need to provide it by setting BodyCreationSettings::mOverrideMassProperties and mMassPropertiesOverride.");
		mInvMass = 1.0f / inMassProperties.mMass;
	}

	if (allowed_rotation_axis == 0)
	{
		// No rotation possible
		mInvInertiaDiagonal = Bivec::sZero();
		mInertiaRotation = Rotor::sIdentity();
	}
	else
	{
		// In 4D the inertia tensor is a 6x6 map bivector->bivector.
		// We store it as: principal bivector diagonal I_ij = M_ii + M_jj (6 values)
		// plus a Rotor aligning the body frame to the principal axes of the 4x4 second moment tensor M.
		Mat44 principal_rotation;
		Vec4 diagonal;
		// A diagonal body-local tensor already supplies a principal frame.
		// Keeping it avoids a permutation -> rotor -> adjoint round trip,
		// which needlessly perturbs even coordinate-aligned null modes.
		bool is_diagonal = true;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				if (i != j && inMassProperties.mInertia(i, j) != 0.0f)
					is_diagonal = false;
		bool decomposed;
		if (is_diagonal)
		{
			principal_rotation = Mat44::sIdentity();
			diagonal = Vec4(inMassProperties.mInertia.GetDiagonal4());
			decomposed = true;
		}
		else
			decomposed = inMassProperties.DecomposePrincipalMomentsOfInertia(principal_rotation, diagonal);
		// Zero is a physical particle moment, not a decomposition failure.
		// An absolute near-zero test would also replace valid small positive
		// moments by invented unit-sphere inertia.
		if (decomposed)
		{
			// Derive the 6 bivector principal moments from the 4 second-moment eigenvalues
			Bivec bivec_inertia = MassProperties::sBivecInertiaFromSecondMomentDiagonal(diagonal);

			// Invert component-wise (guard against zero components)
			auto safe_inv = [](float v) { return v > 0.0f ? 1.0f / v : 0.0f; };
			mInvInertiaDiagonal = Bivec(
				safe_inv(bivec_inertia.GetE12()),
				safe_inv(bivec_inertia.GetE13()),
				safe_inv(bivec_inertia.GetE14()),
				safe_inv(bivec_inertia.GetE23()),
				safe_inv(bivec_inertia.GetE24()),
				safe_inv(bivec_inertia.GetE34()));

			// Convert 4x4 orthogonal rotation to Rotor
			mInertiaRotation = principal_rotation.GetRotor();
		}
		else
		{
			// Failed to decompose. Fall back to a uniform 4D hypersphere of radius 1.
			// Solid 4D unit hypersphere second moment tensor is diagonal M_ii = m / 6,
			// so bivector inertia is I_ij = m / 3.
			float fallback_inv = 3.0f * mInvMass;
			mInvInertiaDiagonal = Bivec(fallback_inv, fallback_inv, fallback_inv, fallback_inv, fallback_inv, fallback_inv);
			mInertiaRotation = Rotor::sIdentity();
		}
	}

	JPH_ASSERT(mInvMass != 0.0f || mInvInertiaDiagonal != Bivec::sZero(), "Can't lock all axes, use a static body for this.");
}

void MotionProperties::SaveState(StateRecorder &inStream) const
{
	// Only write properties that can change at runtime
	inStream.Write(mLinearVelocity);
	inStream.Write(mAngularVelocity);
	inStream.Write(mForce);
	inStream.Write(mTorque);
#ifdef JPH_DOUBLE_PRECISION
	inStream.Write(mSleepTestOffset);
#endif // JPH_DOUBLE_PRECISION
	inStream.Write(mSleepTestSpheres);
	inStream.Write(mSleepTestTimer);
	inStream.Write(mAllowSleeping);
}

void MotionProperties::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mLinearVelocity);
	inStream.Read(mAngularVelocity);
	inStream.Read(mForce);
	inStream.Read(mTorque);
#ifdef JPH_DOUBLE_PRECISION
	inStream.Read(mSleepTestOffset);
#endif // JPH_DOUBLE_PRECISION
	inStream.Read(mSleepTestSpheres);
	inStream.Read(mSleepTestTimer);
	inStream.Read(mAllowSleeping);
}

JPH_NAMESPACE_END
