// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Body/Body.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
// TODO(4D): SoftBody is a documented v1 scope cut (CLAUDE.md "Scope cuts for v1"). The SoftBody
// classes are still 3D and their sources are excluded from Jolt.cmake, so every reference to them
// from the rigid body core is compiled out behind JPH_ENABLE_SOFT_BODY, which is never defined for
// v1. The code is left in place verbatim so it can be turned back on once SoftBody is ported.
#ifdef JPH_ENABLE_SOFT_BODY
	#include <Jolt/Physics/SoftBody/SoftBodyCreationSettings.h>
	#include <Jolt/Physics/SoftBody/SoftBodyMotionProperties.h>
#endif // JPH_ENABLE_SOFT_BODY
#include <Jolt/Physics/PhysicsSettings.h>
#include <Jolt/Physics/StateRecorder.h>
#include <Jolt/Physics/Collision/Shape/EmptyShape.h>
#include <Jolt/Core/StringTools.h>
#include <Jolt/Core/Profiler.h>
#include <Jolt/Core/InsertionSort.h>
#include <Jolt/Math/Matrix.h>
#include <Jolt/Math/Vector.h>
#include <Jolt/Math/EigenValueSymmetric.h>
#include <cmath>
#include <cstdlib>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_DEBUG_RENDERER

JPH_NAMESPACE_BEGIN

static const EmptyShape sFixedToWorldShape;
Body Body::sFixedToWorld(false);

Body::Body(bool) :
	mPosition(RVec4::sZero()),
	mRotation(Rotor::sIdentity()),
	mShape(&sFixedToWorldShape), // Dummy shape
	mFriction(0.0f),
	mRestitution(0.0f),
	mObjectLayer(cObjectLayerInvalid),
	mMotionType(EMotionType::Static)
{
	sFixedToWorldShape.SetEmbedded();
}

void Body::SetMotionType(EMotionType inMotionType)
{
	if (mMotionType == inMotionType)
		return;

	JPH_ASSERT(inMotionType == EMotionType::Static || mMotionProperties != nullptr, "Body needs to be created with mAllowDynamicOrKinematic set to true");
	JPH_ASSERT(inMotionType != EMotionType::Static || !IsActive(), "Deactivate body first");
	JPH_ASSERT(inMotionType == EMotionType::Dynamic || !IsSoftBody(), "Soft bodies can only be dynamic, you can make individual vertices kinematic by setting their inverse mass to 0");

	// Store new motion type
	mMotionType = inMotionType;

	if (mMotionProperties != nullptr)
	{
		// Update cache
		JPH_IF_ENABLE_ASSERTS(mMotionProperties->mCachedMotionType = inMotionType;)

		switch (inMotionType)
		{
		case EMotionType::Static:
			// Stop the object
			mMotionProperties->mLinearVelocity = Vec4::sZero();
			mMotionProperties->mAngularVelocity = Bivec::sZero();
			[[fallthrough]];

		case EMotionType::Kinematic:
			// Cancel forces
			mMotionProperties->ResetForce();
			mMotionProperties->ResetTorque();
			break;

		case EMotionType::Dynamic:
			break;
		}
	}
}

void Body::SetAllowSleeping(bool inAllow)
{
	mMotionProperties->mAllowSleeping = inAllow;
	if (inAllow)
		ResetSleepTimer();
}

void Body::MoveKinematic(RVec4Arg inTargetPosition, RotorArg inTargetRotation, float inDeltaTime)
{
	JPH_ASSERT(IsRigidBody()); // Only valid for rigid bodies
	JPH_ASSERT(!IsStatic());
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::Read));

	// Calculate center of mass at end situation
	RVec4 new_com = inTargetPosition + RVec4(Mat44::sRotation(inTargetRotation) * Lane4(mShape->GetCenterOfMass()));

	// Calculate delta position and rotation
	Vec4 delta_pos = Vec4(new_com - mPosition);
	Rotor delta_rotation = inTargetRotation * mRotation.Reversed();

	mMotionProperties->MoveKinematic(delta_pos, delta_rotation, inDeltaTime);
}

void Body::CalculateWorldSpaceBoundsInternal()
{
	mBounds = mShape->GetWorldSpaceBounds(GetCenterOfMassTransform(), Vec4::sOne());
}

void Body::SetPositionAndRotationInternal(RVec4Arg inPosition, RotorArg inRotation, bool inResetSleepTimer)
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::ReadWrite));

	mPosition = inPosition + RVec4(Mat44::sRotation(inRotation) * Lane4(mShape->GetCenterOfMass()));
	mRotation = inRotation;

	// Initialize bounding box
	CalculateWorldSpaceBoundsInternal();

	// Reset sleeping test
	if (inResetSleepTimer && mMotionProperties != nullptr)
		ResetSleepTimer();
}

void Body::UpdateCenterOfMassInternal(Vec4Arg inPreviousCenterOfMass, bool inUpdateMassProperties)
{
	// Update center of mass position so the world position for this body stays the same
	mPosition += RVec4(Mat44::sRotation(mRotation) * Lane4(mShape->GetCenterOfMass() - inPreviousCenterOfMass));

	// Recalculate mass and inertia if requested
	if (inUpdateMassProperties && mMotionProperties != nullptr)
		mMotionProperties->SetMassProperties(mMotionProperties->GetAllowedDOFs(), mShape->GetMassProperties());
}

void Body::SetShapeInternal(const Shape *inShape, bool inUpdateMassProperties)
{
	JPH_ASSERT(IsRigidBody()); // Only valid for rigid bodies
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::ReadWrite));

	// Get the old center of mass
	Vec4 old_com = mShape->GetCenterOfMass();

	// Update the shape
	mShape = inShape;

	// Update center of mass
	UpdateCenterOfMassInternal(old_com, inUpdateMassProperties);

	// Recalculate bounding box
	CalculateWorldSpaceBoundsInternal();
}

ECanSleep Body::UpdateSleepStateInternal(float inDeltaTime, float inMaxMovement, float inTimeBeforeSleep)
{
	// Check override & sensors will never go to sleep (they would stop detecting collisions with sleeping bodies)
	if (!mMotionProperties->mAllowSleeping || IsSensor())
		return ECanSleep::CannotSleep;

	// Get the points to test
	RVec4 points[3];
	GetSleepTestPoints(points);

#ifdef JPH_DOUBLE_PRECISION
	// Get base offset for spheres
	DVec4 offset = mMotionProperties->GetSleepTestOffset();
#endif // JPH_DOUBLE_PRECISION

	for (int i = 0; i < 3; ++i)
	{
		Sphere &sphere = mMotionProperties->mSleepTestSpheres[i];

		// Make point relative to base offset
#ifdef JPH_DOUBLE_PRECISION
		DVec4 diff = points[i] - offset;
		Vec4 p(float(diff.GetX()), float(diff.GetY()), float(diff.GetZ()), float(diff.GetW()));
#else
		Vec4 p = points[i];
#endif // JPH_DOUBLE_PRECISION

		// Encapsulate the point in a sphere
		sphere.EncapsulatePoint(p);

		// Test if it exceeded the max movement
		if (sphere.GetRadius() > inMaxMovement)
		{
			// Body is not sleeping, reset test
			mMotionProperties->ResetSleepTestSpheres(points);
			return ECanSleep::CannotSleep;
		}
	}

	return mMotionProperties->AccumulateSleepTime(inDeltaTime, inTimeBeforeSleep);
}

void Body::GetSubmergedVolume(RVec4Arg inSurfacePosition, Vec4Arg inSurfaceNormal, float &outTotalVolume, float &outSubmergedVolume, Vec4 &outRelativeCenterOfBuoyancy) const
{
	// For GetSubmergedVolume we transform the surface relative to the body position for increased precision
	Mat44 rotation = Mat44::sRotation(mRotation);
	Plane surface_relative_to_body = Plane::sFromPointAndNormal(Vec4(inSurfacePosition - mPosition), inSurfaceNormal);

	// Calculate amount of volume that is submerged and what the center of buoyancy is
	// TODO(4D): inBaseOffset is unconditional in the 4D Shape interface (it was debug-renderer only
	// upstream), so it is always passed here
	mShape->GetSubmergedVolume(RMat44(rotation, mPosition), Vec4::sOne(), surface_relative_to_body, outTotalVolume, outSubmergedVolume, outRelativeCenterOfBuoyancy, mPosition);
}

bool Body::ApplyBuoyancyImpulse(float inTotalVolume, float inSubmergedVolume, Vec4Arg inRelativeCenterOfBuoyancy, float inBuoyancy, float inLinearDrag, float inAngularDrag, Vec4Arg inFluidVelocity, Vec4Arg inGravity, float inDeltaTime)
{
	JPH_ASSERT(IsRigidBody()); // Only implemented for rigid bodies currently

	// We follow the approach from 'Game Programming Gems 6' 2.5 Exact Buoyancy for Polyhedra
	// All quantities below are in world space

	// If we're not submerged, there's no point in doing the rest of the calculations
	if (inSubmergedVolume > 0.0f)
	{
	#ifdef JPH_DEBUG_RENDERER
		// Draw submerged volume properties
		if (Shape::sDrawSubmergedVolumes)
		{
			RVec4 center_of_buoyancy = mPosition + RVec4(Lane4(inRelativeCenterOfBuoyancy));
			DebugRenderer::sInstance->DrawMarker(center_of_buoyancy, Color::sWhite, 2.0f);
			DebugRenderer::sInstance->DrawText3D(center_of_buoyancy, StringFormat("%.3f / %.3f", (double)inSubmergedVolume, (double)inTotalVolume));
		}
	#endif // JPH_DEBUG_RENDERER

		// When buoyancy is 1 we want neutral buoyancy, this means that the density of the liquid is the same as the density of the body at that point.
		// Buoyancy > 1 should make the object float, < 1 should make it sink.
		float inverse_mass = mMotionProperties->GetInverseMass();
		float fluid_density = inBuoyancy / (inTotalVolume * inverse_mass);

		// Buoyancy force = Density of Fluid * Submerged volume * Magnitude of gravity * Up direction (eq 2.5.1)
		// Impulse = Force * Delta time
		// We should apply this at the center of buoyancy (= center of mass of submerged volume)
		Vec4 buoyancy_impulse = -fluid_density * inSubmergedVolume * mMotionProperties->GetGravityFactor() * inGravity * inDeltaTime;

		// Calculate the velocity of the center of buoyancy relative to the fluid
		Vec4 linear_velocity = mMotionProperties->GetLinearVelocity();
		Bivec angular_velocity = mMotionProperties->GetAngularVelocity();
		Vec4 center_of_buoyancy_velocity = linear_velocity + angular_velocity.Contract(inRelativeCenterOfBuoyancy);
		Vec4 relative_center_of_buoyancy_velocity = inFluidVelocity - center_of_buoyancy_velocity;

		// Here we deviate from the article, instead of eq 2.5.14 we use a quadratic drag formula: https://en.wikipedia.org/wiki/Drag_%28physics%29
		// Drag force = 0.5 * Fluid Density * (Velocity of fluid - Velocity of center of buoyancy)^2 * Linear Drag * Area Facing the Relative Fluid Velocity

		// Get size of local bounding box
		Vec4 size = mShape->GetLocalBounds().GetSize();

		// Determine projected area of the local space bounding box in the direction of the relative velocity
		// In 4D, the "area" facing a direction is a 3D cross-section (hyperarea)
		float area = 0.0f;
		float relative_center_of_buoyancy_velocity_len_sq = relative_center_of_buoyancy_velocity.LengthSq();
		if (relative_center_of_buoyancy_velocity_len_sq > 1.0e-12f)
		{
			// Approximate cross-sectional hyperarea: project box onto velocity direction
			// For a 4D box, the cross section facing direction d has hyperarea = product of 3 edge lengths
			// perpendicular to d. Approximate by: sum(|d_i| * product of other 3 sizes) / |d|
			Vec4 local_vel = Vec4(Mat44::sRotation(GetRotation()).Transposed() * Lane4(relative_center_of_buoyancy_velocity));
			Vec4 abs_vel = local_vel.Abs();
			float inv_len = 1.0f / sqrt(relative_center_of_buoyancy_velocity_len_sq);
			// Cross-sectional area perpendicular to each axis
			float sx = size.GetX(), sy = size.GetY(), sz = size.GetZ(), sw = size.GetW();
			area = inv_len * (abs_vel.GetX() * sy * sz * sw
			                + abs_vel.GetY() * sx * sz * sw
			                + abs_vel.GetZ() * sx * sy * sw
			                + abs_vel.GetW() * sx * sy * sz);
		}

		// Calculate the impulse
		Vec4 drag_impulse = (0.5f * fluid_density * inLinearDrag * area * inDeltaTime) * relative_center_of_buoyancy_velocity * relative_center_of_buoyancy_velocity.Length();

		// Clamp magnitude against current linear velocity to prevent overshoot
		float linear_velocity_len_sq = linear_velocity.LengthSq();
		float drag_delta_linear_velocity_len_sq = (drag_impulse * inverse_mass).LengthSq();
		if (drag_delta_linear_velocity_len_sq > linear_velocity_len_sq)
			drag_impulse *= sqrt(linear_velocity_len_sq / drag_delta_linear_velocity_len_sq);

		// Calculate the resulting delta linear velocity due to buoyancy and drag
		Vec4 delta_linear_velocity = (drag_impulse + buoyancy_impulse) * inverse_mass;
		mMotionProperties->AddLinearVelocityStep(delta_linear_velocity);

		// Determine average width of the body (across the four axes)
		float l = (size.GetX() + size.GetY() + size.GetZ() + size.GetW()) / 4.0f;

		// Drag torque (bivector) = -Angular Drag * Mass * Submerged volume / Total volume * (Average width)^2 * Angular velocity
		Bivec drag_angular_impulse = (-inAngularDrag * inSubmergedVolume / inTotalVolume * inDeltaTime * Square(l) / inverse_mass) * angular_velocity;
		Bivec drag_delta_angular_velocity = mMotionProperties->MultiplyWorldSpaceInverseInertiaByVector(mRotation, drag_angular_impulse);

		// Clamp magnitude against the current angular velocity to prevent overshoot
		float angular_velocity_len_sq = angular_velocity.LengthSq();
		float drag_delta_angular_velocity_len_sq = drag_delta_angular_velocity.LengthSq();
		if (drag_delta_angular_velocity_len_sq > angular_velocity_len_sq)
			drag_delta_angular_velocity *= sqrt(angular_velocity_len_sq / drag_delta_angular_velocity_len_sq);

		// Calculate total delta angular velocity due to drag and buoyancy
		// Torque from buoyancy+drag forces at center of buoyancy: (buoyancy_impulse + drag_impulse) ^ r
		// (impulse first -- see the operand order note in Body::AddForce)
		Bivec buoyancy_angular_impulse = Bivec::sWedge(buoyancy_impulse + drag_impulse, inRelativeCenterOfBuoyancy);
		Bivec delta_angular_velocity = drag_delta_angular_velocity + mMotionProperties->MultiplyWorldSpaceInverseInertiaByVector(mRotation, buoyancy_angular_impulse);
		mMotionProperties->AddAngularVelocityStep(delta_angular_velocity);
		return true;
	}

	return false;
}

bool Body::ApplyBuoyancyImpulse(RVec4Arg inSurfacePosition, Vec4Arg inSurfaceNormal, float inBuoyancy, float inLinearDrag, float inAngularDrag, Vec4Arg inFluidVelocity, Vec4Arg inGravity, float inDeltaTime)
{
	JPH_PROFILE_FUNCTION();

	float total_volume, submerged_volume;
	Vec4 relative_center_of_buoyancy;
	GetSubmergedVolume(inSurfacePosition, inSurfaceNormal, total_volume, submerged_volume, relative_center_of_buoyancy);

	return ApplyBuoyancyImpulse(total_volume, submerged_volume, relative_center_of_buoyancy, inBuoyancy, inLinearDrag, inAngularDrag, inFluidVelocity, inGravity, inDeltaTime);
}

void Body::SaveState(StateRecorder &inStream) const
{
	// Only write properties that can change at runtime
	inStream.Write(mPosition);
	inStream.Write(mRotation);

	if (mMotionProperties != nullptr)
	{
#ifdef JPH_ENABLE_SOFT_BODY
		if (IsSoftBody())
			static_cast<const SoftBodyMotionProperties *>(mMotionProperties)->SaveState(inStream);
		else
#endif // JPH_ENABLE_SOFT_BODY
			mMotionProperties->SaveState(inStream);
	}
}

void Body::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mPosition);
	inStream.Read(mRotation);

	if (mMotionProperties != nullptr)
	{
#ifdef JPH_ENABLE_SOFT_BODY
		if (IsSoftBody())
			static_cast<SoftBodyMotionProperties *>(mMotionProperties)->RestoreState(inStream);
		else
#endif // JPH_ENABLE_SOFT_BODY
			mMotionProperties->RestoreState(inStream);

		JPH_IF_ENABLE_ASSERTS(mMotionProperties->mCachedMotionType = mMotionType);
	}

	// Initialize bounding box
	CalculateWorldSpaceBoundsInternal();
}

BodyCreationSettings Body::GetBodyCreationSettings() const
{
	BodyCreationSettings result;
	if (!TryGetBodyCreationSettings(result))
	{
		// This legacy API has no error channel. Never return plausible settings
		// which silently turn a null rotational mode into a huge inverse.
		Trace("Body::GetBodyCreationSettings: unrepresentable mass state; use TryGetBodyCreationSettings");
		std::abort();
	}
	return result;
}

bool Body::TryGetBodyCreationSettings(BodyCreationSettings &outSettings) const
{
	JPH_ASSERT(IsRigidBody());

	BodyCreationSettings result;

	result.mPosition = GetPosition();
	result.mRotation = GetRotation();
	result.mLinearVelocity = mMotionProperties != nullptr? mMotionProperties->GetLinearVelocity() : Vec4::sZero();
	result.mAngularVelocity = mMotionProperties != nullptr? mMotionProperties->GetAngularVelocity() : Bivec::sZero();
	result.mObjectLayer = GetObjectLayer();
	result.mUserData = mUserData;
	result.mCollisionGroup = GetCollisionGroup();
	result.mMotionType = GetMotionType();
	result.mAllowedDOFs = mMotionProperties != nullptr? mMotionProperties->GetAllowedDOFs() : EAllowedDOFs::All;
	result.mAllowDynamicOrKinematic = mMotionProperties != nullptr;
	result.mIsSensor = IsSensor();
	result.mCollideKinematicVsNonDynamic = GetCollideKinematicVsNonDynamic();
	result.mUseManifoldReduction = GetUseManifoldReduction();
	result.mApplyGyroscopicForce = GetApplyGyroscopicForce();
	result.mMotionQuality = mMotionProperties != nullptr? mMotionProperties->GetMotionQuality() : EMotionQuality::Discrete;
	result.mEnhancedInternalEdgeRemoval = GetEnhancedInternalEdgeRemoval();
	result.mAllowSleeping = mMotionProperties != nullptr? GetAllowSleeping() : true;
	result.mFriction = GetFriction();
	result.mRestitution = GetRestitution();
	result.mLinearDamping = mMotionProperties != nullptr? mMotionProperties->GetLinearDamping() : 0.0f;
	result.mAngularDamping = mMotionProperties != nullptr? mMotionProperties->GetAngularDamping() : 0.0f;
	result.mMaxLinearVelocity = mMotionProperties != nullptr? mMotionProperties->GetMaxLinearVelocity() : 0.0f;
	result.mMaxAngularVelocity = mMotionProperties != nullptr? mMotionProperties->GetMaxAngularVelocity() : 0.0f;
	result.mGravityFactor = mMotionProperties != nullptr? mMotionProperties->GetGravityFactor() : 1.0f;
	result.mNumVelocityStepsOverride = mMotionProperties != nullptr? mMotionProperties->GetNumVelocityStepsOverride() : 0;
	result.mNumPositionStepsOverride = mMotionProperties != nullptr? mMotionProperties->GetNumPositionStepsOverride() : 0;
	result.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;

	// Only export a state that survives the actual float mass encode/decode.
	// Oblique null modes cannot in general be encoded by a float Mat44
	// without changing rank. Reject, rather than epsilon-clamping a moment.
	// See Docs/MassRoundTrip.md for the supported domain and error budget.
	// Static bodies without motion properties do not retain mass/inertia;
	// fully locked sectors discard that information too. Use finite canonical
	// values for those unobservable sectors and retain the allowed DOFs.
	result.mMassPropertiesOverride.mMass = 1.0f;
	result.mMassPropertiesOverride.mInertia = Mat44::sZero();
	if (mMotionProperties != nullptr)
	{
		float inv_mass = mMotionProperties->GetInverseMassUnchecked();
		if ((result.mAllowedDOFs & EAllowedDOFs::AllTranslation) != EAllowedDOFs::None)
		{
			if (!(inv_mass > 0.0f) || !std::isfinite(inv_mass) || 1.0 / inv_mass > FLT_MAX)
				return false;
			result.mMassPropertiesOverride.mMass = float(1.0 / inv_mass);
		}

		if ((result.mAllowedDOFs & EAllowedDOFs::AllRotation) != EAllowedDOFs::None)
		{
			// Body is a friend of MotionProperties: read the stored principal
			// diagonal directly, also for static/kinematic bodies that retain it
			// for later motion-type changes. The public getter is dynamic-only.
			Bivec inv_diag = mMotionProperties->mInvInertiaDiagonal;
			double moment[6];
			double trace = 0.0;
			int zero_moments = 0;
			for (int plane = 0; plane < 6; ++plane)
			{
				if (!(inv_diag[plane] >= 0.0f) || !std::isfinite(inv_diag[plane]))
					return false;
				// A zero planar moment has zero inverse by the particle convention,
				// not an infinite physical second moment.
				moment[plane] = inv_diag[plane] > 0.0f? 1.0 / inv_diag[plane] : 0.0;
				zero_moments += inv_diag[plane] == 0.0f;
				trace += moment[plane];
			}
			// I_ij = d_i + d_j. Each d_i occurs in three of the six
			// moments, so trace(M) = sum(I_ij)/3 and
			// d_i = (sum of the three incident moments - trace(M))/2.
			trace /= 3.0;
			double d[4] = {
				0.5 * (moment[0] + moment[1] + moment[2] - trace),
				0.5 * (moment[0] + moment[3] + moment[4] - trace),
				0.5 * (moment[1] + moment[3] + moment[5] - trace),
				0.5 * (moment[2] + moment[4] + moment[5] - trace) };
			if (zero_moments != 0)
			{
				// PSD plus I_ab = d_a + d_b = 0 implies d_a = d_b = 0
				// EXACTLY. Solve from an incident moment instead of subtracting
				// nearly equal reciprocals. This is not an epsilon rank decision.
				constexpr int a[] = { 0, 0, 0, 1, 1, 2 }, b[] = { 1, 2, 3, 2, 3, 3 };
				for (int plane = 0; plane < 6; ++plane)
					if (inv_diag[plane] == 0.0f)
					{
						int axis = a[plane];
						d[axis] = d[b[plane]] = 0.0;
						for (int incident = 0; incident < 6; ++incident)
						{
							if (a[incident] == axis)
								d[b[incident]] = moment[incident];
							else if (b[incident] == axis)
								d[a[incident]] = moment[incident];
						}
						break;
					}
			}
			for (double value : d)
				if (!(value >= 0.0) || value > FLT_MAX)
					return false;
			// This is body-local inertia; body rotation is exported separately.
			Mat44 rotation = Mat44::sRotation(mMotionProperties->GetInertiaRotation());
			Mat44 &encoded = result.mMassPropertiesOverride.mInertia;
			for (int i = 0; i < 4; ++i)
				for (int j = i; j < 4; ++j)
				{
					double value = 0.0;
					for (int k = 0; k < 4; ++k)
						value += double(rotation(i, k)) * d[k] * rotation(j, k);
					if (!std::isfinite(value) || std::abs(value) > FLT_MAX)
						return false;
					encoded(i, j) = encoded(j, i) = float(value);
					// Exact, deliberately conservative supported domain for null
					// modes. An oblique float tensor cannot certify their rank.
					if (zero_moments != 0 && i != j && encoded(i, j) != 0.0f)
						return false;
				}

			// Preflight the eigensolver and reconstruction diagnostics before
			// entering the real decoder. Checked eigenpairs use the SAME solver
			// and tolerances, but return false in all builds rather than asserting.
			// The diagonal path in SetMassProperties does not use that solver.
			bool is_diagonal = true;
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					is_diagonal = is_diagonal && (i == j || encoded(i, j) == 0.0f);
			if (!is_diagonal)
			{
				Matrix<4, 4> inertia;
				inertia.CopyPart(encoded, 0, 0, 4, 4, 0, 0);
				Matrix<4, 4> eigen_vec = Matrix<4, 4>::sIdentity();
				Vector<4> eigen_val;
				if (!EigenValueSymmetric(inertia, eigen_vec, eigen_val, true))
					return false;
				int indices[] = { 0, 1, 2, 3 };
				InsertionSort(indices, indices + 4, [&eigen_val](int a, int b) { return eigen_val[a] > eigen_val[b]; });
				Mat44 principal_rotation, principal_diagonal = Mat44::sZero();
				for (int i = 0; i < 4; ++i)
				{
					if (!(eigen_val[i] >= 0.0f) || !std::isfinite(eigen_val[i]))
						return false;
					const auto &col = eigen_vec.GetColumn(indices[i]);
					principal_rotation.SetColumn4(i, Lane4(col[0], col[1], col[2], col[3]));
					principal_diagonal(i, i) = eigen_val[indices[i]];
				}
				if (principal_rotation.GetDeterminant() < 0.0f)
					principal_rotation.SetColumn4(3, -principal_rotation.GetColumn4(3));
				Mat44 reconstructed = principal_rotation * principal_diagonal * principal_rotation.Transposed();
				for (int i = 0; i < 4; ++i)
				{
					Vec4 original_column(encoded.GetColumn4(i));
					if (!Vec4(reconstructed.GetColumn4(i)).IsClose(original_column, max(original_column.LengthSq(), 1.0f) * 1.0e-6f))
						return false;
				}
			}

			// Dry-run the real decoder only after checked preflight succeeds.
			// No second independently editable dynamical state is retained.
			// Accept only finite inertia, unchanged nullity and response budget.
			MotionProperties decoded;
			decoded.SetMassProperties(EAllowedDOFs::All, result.mMassPropertiesOverride);
			int decoded_zeros = 0;
			for (int k = 0; k < 6; ++k)
			{
				float value = decoded.mInvInertiaDiagonal[k];
				if (!(value >= 0.0f) || !std::isfinite(value))
					return false;
				decoded_zeros += value == 0.0f;
			}
			if (decoded_zeros != zero_moments)
				return false;
			double error_sq = 0.0, reference_sq = 0.0;
			for (int input = 0; input < 6; ++input)
			{
				float v[6] = {};
				v[input] = 1.0f;
				Bivec basis(v[0], v[1], v[2], v[3], v[4], v[5]);
				Rotor original_rotation = mMotionProperties->GetInertiaRotation();
				Bivec expected = original_rotation.RotateBivec(inv_diag * original_rotation.InverseRotateBivec(basis));
				Bivec actual = decoded.mInertiaRotation.RotateBivec(decoded.mInvInertiaDiagonal * decoded.mInertiaRotation.InverseRotateBivec(basis));
				for (int output = 0; output < 6; ++output)
				{
					if (!std::isfinite(actual[output]) || !std::isfinite(expected[output]))
						return false;
					error_sq += Square(double(actual[output]) - expected[output]);
					reference_sq += Square(double(expected[output]));
				}
			}
			if (!(error_sq <= Square(64.0 * FLT_EPSILON) * reference_sq))
				return false;
		}
	}

	result.SetShape(GetShape());

	outSettings = std::move(result);
	return true;
}

#ifdef JPH_ENABLE_SOFT_BODY
SoftBodyCreationSettings Body::GetSoftBodyCreationSettings() const
{
	JPH_ASSERT(IsSoftBody());

	SoftBodyCreationSettings result;

	result.mPosition = GetPosition();
	result.mRotation = GetRotation();
	result.mUserData = mUserData;
	result.mObjectLayer = GetObjectLayer();
	result.mCollisionGroup = GetCollisionGroup();
	result.mFriction = GetFriction();
	result.mRestitution = GetRestitution();
	const SoftBodyMotionProperties *mp = static_cast<const SoftBodyMotionProperties *>(mMotionProperties);
	result.mNumIterations = mp->GetNumIterations();
	result.mLinearDamping = mp->GetLinearDamping();
	result.mMaxLinearVelocity = mp->GetMaxLinearVelocity();
	result.mGravityFactor = mp->GetGravityFactor();
	result.mPressure = mp->GetPressure();
	result.mUpdatePosition = mp->GetUpdatePosition();
	result.mVertexRadius = mp->GetVertexRadius();
	result.mAllowSleeping = mp->GetAllowSleeping();
	result.mFacesDoubleSided = mp->GetFacesDoubleSided();
	result.mSettings = mp->GetSettings();

	return result;
}
#endif // JPH_ENABLE_SOFT_BODY

JPH_NAMESPACE_END
