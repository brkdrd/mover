// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

RMat44 Body::GetWorldTransform() const
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::Read));

	return RMat44::sRotationTranslation(mRotation, mPosition).PreTranslated(-mShape->GetCenterOfMass());
}

RMat44 Body::GetCenterOfMassTransform() const
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::Read));

	return RMat44::sRotationTranslation(mRotation, mPosition);
}

RMat44 Body::GetInverseCenterOfMassTransform() const
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::Read));

	return RMat44::sInverseRotationTranslation(mRotation, mPosition);
}

inline bool Body::sFindCollidingPairsCanCollide(const Body &inBody1, const Body &inBody2)
{
	// First body should never be a soft body
	JPH_ASSERT(!inBody1.IsSoftBody());

	// One of these conditions must be true
	// - We always allow detecting collisions between kinematic and non-dynamic bodies
	// - One of the bodies must be dynamic to collide
	// - A kinematic object can collide with a sensor
	if (!inBody1.GetCollideKinematicVsNonDynamic()
		&& !inBody2.GetCollideKinematicVsNonDynamic()
		&& (!inBody1.IsDynamic() && !inBody2.IsDynamic())
		&& !(inBody1.IsKinematic() && inBody2.IsSensor())
		&& !(inBody2.IsKinematic() && inBody1.IsSensor()))
		return false;

	// Check that body 1 is active
	uint32 body1_index_in_active_bodies = inBody1.GetIndexInActiveBodiesInternal();
	JPH_ASSERT(!inBody1.IsStatic() && body1_index_in_active_bodies != Body::cInactiveIndex, "This function assumes that Body 1 is active");

	// See original comment for ordering explanation
	static_assert(Body::cInactiveIndex == 0xffffffff, "The algorithm below uses this value");
	if (!inBody2.IsSoftBody() && body1_index_in_active_bodies >= inBody2.GetIndexInActiveBodiesInternal())
		return false;
	JPH_ASSERT(inBody1.GetID() != inBody2.GetID(), "Read the comment above, A and B are the same body which should not be possible!");

	// Check collision group filter
	if (!inBody1.GetCollisionGroup().CanCollide(inBody2.GetCollisionGroup()))
		return false;

	return true;
}

void Body::AddRotationStep(BivecArg inAngularVelocityTimesDeltaTime)
{
	JPH_ASSERT(IsRigidBody());
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::ReadWrite));

	// In 4D, angular velocity is a bivector. Integrate rotation as
	//   R(t + dt) = exp(Omega * dt / 2) * R(t)
	// using the general rotor exponential, which is correct for double / isoclinic rotations (the
	// old inline "cos|B| + sin|B|/|B| * B" formula dropped the e1234 term and only handled simple
	// rotations). inAngularVelocityTimesDeltaTime = Omega * dt.
	Rotor delta = Rotor::sExp(0.5f * inAngularVelocityTimesDeltaTime);
	mRotation = (delta * mRotation).Normalized();
	JPH_ASSERT(!mRotation.IsNaN());
}

void Body::SubRotationStep(BivecArg inAngularVelocityTimesDeltaTime)
{
	JPH_ASSERT(IsRigidBody());
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::ReadWrite));

	// Same as AddRotationStep but with negated bivector
	Rotor delta = Rotor::sExp(-0.5f * inAngularVelocityTimesDeltaTime);
	mRotation = (delta * mRotation).Normalized();
	JPH_ASSERT(!mRotation.IsNaN());
}

Vec4 Body::GetWorldSpaceSurfaceNormal(const SubShapeID &inSubShapeID, RVec4Arg inPosition) const
{
	RMat44 inv_com = GetInverseCenterOfMassTransform();
	return inv_com.Multiply3x3Transposed(mShape->GetSurfaceNormal(inSubShapeID, Vec4(inv_com * inPosition))).Normalized();
}

void Body::AddForce(Vec4Arg inForce, RVec4Arg inPosition)
{
	AddForce(inForce);
	// Torque = F ^ r, the bivector conjugate to the angular velocity. NOT r ^ F: in 4D the point
	// velocity operator Bivec::Contract is the *negative* adjoint of the wedge, so the 4D analog of
	// the 3D torque r x F has its operands swapped. See UnitTests/Physics/DynamicsConventionsTest.
	AddTorque(Bivec::sWedge(inForce, Vec4(inPosition - mPosition)));
}

void Body::AddImpulse(Vec4Arg inImpulse)
{
	JPH_ASSERT(IsDynamic());

	SetLinearVelocityClamped(mMotionProperties->GetLinearVelocity() + inImpulse * mMotionProperties->GetInverseMass());
}

void Body::AddImpulse(Vec4Arg inImpulse, RVec4Arg inPosition)
{
	JPH_ASSERT(IsDynamic());

	SetLinearVelocityClamped(mMotionProperties->GetLinearVelocity() + inImpulse * mMotionProperties->GetInverseMass());

	// Angular impulse = impulse ^ r (see the operand order note in AddForce above)
	Bivec angular_impulse = Bivec::sWedge(inImpulse, Vec4(inPosition - mPosition));
	SetAngularVelocityClamped(mMotionProperties->GetAngularVelocity() + mMotionProperties->MultiplyWorldSpaceInverseInertiaByVector(mRotation, angular_impulse));
}

void Body::AddAngularImpulse(BivecArg inAngularImpulse)
{
	JPH_ASSERT(IsDynamic());

	SetAngularVelocityClamped(mMotionProperties->GetAngularVelocity() + mMotionProperties->MultiplyWorldSpaceInverseInertiaByVector(mRotation, inAngularImpulse));
}

void Body::GetSleepTestPoints(RVec4 *outPoints) const
{
	JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sPositionAccess(), BodyAccess::EAccess::Read));

	// Center of mass is the first position
	outPoints[0] = mPosition;

	// The second and third position are on the largest axis of the bounding box
	Vec4 extent = mShape->GetLocalBounds().GetExtent();
	int lowest_component = extent.GetLowestComponentIndex();
	Mat44 rotation = Mat44::sRotation(mRotation);

	// In 4D we have 4 axes; pick the 2 largest (skip the smallest and second-smallest)
	// For backward compatibility we still use 3 test points
	// Find the second smallest component index
	int second_lowest = -1;
	float second_lowest_val = FLT_MAX;
	for (int i = 0; i < 4; ++i)
	{
		if (i != lowest_component && extent[i] < second_lowest_val)
		{
			second_lowest_val = extent[i];
			second_lowest = i;
		}
	}

	// Place points along the two largest axes
	int point_idx = 1;
	for (int i = 0; i < 4 && point_idx < 3; ++i)
	{
		if (i != lowest_component && i != second_lowest)
			outPoints[point_idx++] = mPosition + RVec4(extent[i] * Vec4(rotation.GetColumn4(i)));
	}
}

void Body::ResetSleepTimer()
{
	RVec4 points[3];
	GetSleepTestPoints(points);
	mMotionProperties->ResetSleepTestSpheres(points);
}

JPH_NAMESPACE_END
