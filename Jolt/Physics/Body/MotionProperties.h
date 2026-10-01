// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Bivec.h>
#include <Jolt/Geometry/Sphere.h>
#include <Jolt/Physics/Body/AllowedDOFs.h>
#include <Jolt/Physics/Body/MotionQuality.h>
#include <Jolt/Physics/Body/BodyAccess.h>
#include <Jolt/Physics/Body/MotionType.h>
#include <Jolt/Physics/Body/BodyType.h>
#include <Jolt/Physics/Body/MassProperties.h>
#include <Jolt/Physics/DeterminismLog.h>

JPH_NAMESPACE_BEGIN

class StateRecorder;

/// TODO(4D): the world-space inverse inertia, as the 6x6 symmetric operator on bivectors that it
/// actually is. Upstream caches the world-space inverse inertia as a Mat44 once per constraint and
/// reuses it for every row and every solver iteration; the 4D port replaced that with a rotor-based
/// recomputation inside MultiplyWorldSpaceInverseInertiaByVector, which costs four full Cl(4,0)
/// sandwich products plus a rotor multiply - a few hundred scalar flops - *per constraint row per
/// body*. A contact has 24 rows (6 manifold points x normal + 3 friction), so that was ~50x more
/// work than building the operator once. This is the 4D analog of upstream's cached Mat44.
///
/// Stored as the 6 images of the basis bivectors, i.e. the columns of Ad(R) D^-1 Ad(R)^T. The
/// operator is symmetric positive semi-definite; the columns are kept rather than a packed upper
/// triangle so that applying it is 6 Vec8 multiply-adds.
class InverseInertiaOperator
{
public:
							InverseInertiaOperator() = default;

	/// All-zero operator, for static / kinematic bodies
	static JPH_INLINE InverseInertiaOperator sZero()
	{
		InverseInertiaOperator op;
		for (int i = 0; i < 6; ++i)
			op.mColumn[i] = Bivec::sZero();
		return op;
	}

	/// Apply the operator to a bivector
	JPH_INLINE Bivec		Multiply(BivecArg inV) const
	{
		Vec8 result = mColumn[0].GetVec8() * Vec8::sReplicate(inV[0]);
		for (int i = 1; i < 6; ++i)
			result = result + mColumn[i].GetVec8() * Vec8::sReplicate(inV[i]);
		return Bivec(result);
	}

	Bivec					mColumn[6];
};

/// Enum that determines if an object can go to sleep
enum class ECanSleep
{
	CannotSleep = 0,																		///< Object cannot go to sleep
	CanSleep = 1,																			///< Object can go to sleep
};

/// The Body class only keeps track of state for static bodies, the MotionProperties class keeps the additional state needed for a moving Body. It has a 1-on-1 relationship with the body.
class JPH_EXPORT MotionProperties
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Motion quality, or how well it detects collisions when it has a high velocity
	EMotionQuality			GetMotionQuality() const										{ return mMotionQuality; }

	/// Get the allowed degrees of freedom that this body has (this can be changed by calling SetMassProperties)
	inline EAllowedDOFs		GetAllowedDOFs() const											{ return mAllowedDOFs; }

	/// If this body can go to sleep.
	inline bool				GetAllowSleeping() const										{ return mAllowSleeping; }

	/// Get world space linear velocity of the center of mass (4D vector)
	inline Vec4				GetLinearVelocity() const										{ JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::Read)); return mLinearVelocity; }

	/// Set world space linear velocity of the center of mass
	void					SetLinearVelocity(Vec4Arg inLinearVelocity)						{ JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite)); JPH_ASSERT(inLinearVelocity.Length() <= mMaxLinearVelocity); mLinearVelocity = LockTranslation(inLinearVelocity); }

	/// Set world space linear velocity of the center of mass, will make sure the value is clamped against the maximum linear velocity
	void					SetLinearVelocityClamped(Vec4Arg inLinearVelocity)				{ mLinearVelocity = LockTranslation(inLinearVelocity); ClampLinearVelocity(); }

	/// Get world space angular velocity (bivector with 6 components for the 6 rotation planes)
	inline Bivec			GetAngularVelocity() const										{ JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::Read)); return mAngularVelocity; }

	/// Set world space angular velocity (bivector)
	void					SetAngularVelocity(BivecArg inAngularVelocity)					{ JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite)); JPH_ASSERT(inAngularVelocity.Length() <= mMaxAngularVelocity); mAngularVelocity = LockAngular(inAngularVelocity); }

	/// Set world space angular velocity, will make sure the value is clamped against the maximum angular velocity
	void					SetAngularVelocityClamped(BivecArg inAngularVelocity)			{ mAngularVelocity = LockAngular(inAngularVelocity); ClampAngularVelocity(); }

	/// Set velocity of body such that it will be rotate/translate by inDeltaPosition/Rotation in inDeltaTime seconds.
	inline void				MoveKinematic(Vec4Arg inDeltaPosition, RotorArg inDeltaRotation, float inDeltaTime);

	///@name Velocity limits
	///@{

	/// Maximum linear velocity that a body can achieve. Used to prevent the system from exploding.
	inline float			GetMaxLinearVelocity() const									{ return mMaxLinearVelocity; }
	inline void				SetMaxLinearVelocity(float inLinearVelocity)					{ JPH_ASSERT(inLinearVelocity >= 0.0f); mMaxLinearVelocity = inLinearVelocity; }

	/// Maximum angular velocity that a body can achieve. Used to prevent the system from exploding.
	inline float			GetMaxAngularVelocity() const									{ return mMaxAngularVelocity; }
	inline void				SetMaxAngularVelocity(float inAngularVelocity)					{ JPH_ASSERT(inAngularVelocity >= 0.0f); mMaxAngularVelocity = inAngularVelocity; }
	///@}

	/// Clamp velocity according to limit
	inline void				ClampLinearVelocity();
	inline void				ClampAngularVelocity();

	/// Get linear damping: dv/dt = -c * v. c. Value should be zero or positive and is usually close to 0.
	inline float			GetLinearDamping() const										{ return mLinearDamping; }
	void					SetLinearDamping(float inLinearDamping)							{ JPH_ASSERT(inLinearDamping >= 0.0f); mLinearDamping = inLinearDamping; }

	/// Get angular damping: dw/dt = -c * w. c. Value should be zero or positive and is usually close to 0.
	inline float			GetAngularDamping() const										{ return mAngularDamping; }
	void					SetAngularDamping(float inAngularDamping)						{ JPH_ASSERT(inAngularDamping >= 0.0f); mAngularDamping = inAngularDamping; }

	/// Get gravity factor (1 = normal gravity, 0 = no gravity)
	inline float			GetGravityFactor() const										{ return mGravityFactor; }
	void					SetGravityFactor(float inGravityFactor)							{ mGravityFactor = inGravityFactor; }

	/// Set the mass and inertia tensor
	void					SetMassProperties(EAllowedDOFs inAllowedDOFs, const MassProperties &inMassProperties);

	/// Get inverse mass (1 / mass). Should only be called on a dynamic object (static or kinematic bodies have infinite mass so should be treated as 1 / mass = 0)
	inline float			GetInverseMass() const											{ JPH_ASSERT(mCachedMotionType == EMotionType::Dynamic); return mInvMass; }
	inline float			GetInverseMassUnchecked() const									{ return mInvMass; }

	/// Set the inverse mass (1 / mass).
	void					SetInverseMass(float inInverseMass)								{ mInvMass = inInverseMass; }

	/// Diagonal of inverse inertia (6 values for the 6 rotation planes). Should only be called on a dynamic object.
	inline Bivec			GetInverseInertiaDiagonal() const								{ JPH_ASSERT(mCachedMotionType == EMotionType::Dynamic); return mInvInertiaDiagonal; }

	/// Rotation (R) that takes inverse inertia diagonal to local space: \f$I_{body}^{-1} = Ad(R) \: D \: Ad(R)^{-1}\f$
	inline Rotor			GetInertiaRotation() const										{ return mInertiaRotation; }

	/// Set the inverse inertia tensor in local space by setting the diagonal and the rotation.
	void					SetInverseInertia(BivecArg inDiagonal, RotorArg inRot)			{ mInvInertiaDiagonal = inDiagonal; mInertiaRotation = inRot; }

	/// Sets the mass to inMass and scale the inertia tensor based on the ratio between the old and new mass.
	void					ScaleToMass(float inMass);

	/// Multiply a bivector with the inverse world space inertia tensor (\f$I_{world}^{-1}\f$). Zero if object is static or kinematic.
	JPH_INLINE Bivec		MultiplyWorldSpaceInverseInertiaByVector(RotorArg inBodyRotation, BivecArg inV) const;

	/// Build the world-space inverse inertia operator for a given body rotation, optionally scaled.
	/// Build this once per body per constraint and apply it with Multiply(), rather than calling
	/// MultiplyWorldSpaceInverseInertiaByVector per row - see InverseInertiaOperator.
	JPH_INLINE InverseInertiaOperator GetInverseInertiaForRotation(RotorArg inBodyRotation, float inInvInertiaScale = 1.0f) const;

	/// Velocity of point inPoint (in center of mass space) of the body (unit: m/s)
	/// In 4D: v = v_linear + Omega._|_r (contraction of angular velocity bivector with position)
	JPH_INLINE Vec4			GetPointVelocityCOM(Vec4Arg inPointRelativeToCOM) const			{ return mLinearVelocity + mAngularVelocity.Contract(inPointRelativeToCOM); }

	// Get the total amount of force applied to the center of mass this time step. Note that it will reset to zero after PhysicsSystem::Update.
	JPH_INLINE Vec4			GetAccumulatedForce() const										{ return Vec4(mForce); }

	// Get the total amount of torque applied to the center of mass this time step (bivector). Note that it will reset to zero after PhysicsSystem::Update.
	JPH_INLINE Bivec		GetAccumulatedTorque() const									{ return Bivec(mTorque); }

	// Reset the total accumulated force, note that this will be done automatically after every time step.
	JPH_INLINE void			ResetForce()													{ mForce = Float4(0, 0, 0, 0); }

	// Reset the total accumulated torque, note that this will be done automatically after every time step.
	JPH_INLINE void			ResetTorque()													{ mTorque = Float8(0, 0, 0, 0, 0, 0, 0, 0); }

	// Reset the current velocity and accumulated force and torque.
	JPH_INLINE void			ResetMotion()
	{
		JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite));
		mLinearVelocity = Vec4::sZero();
		mAngularVelocity = Bivec::sZero();
		mForce = Float4(0, 0, 0, 0);
		mTorque = Float8(0, 0, 0, 0, 0, 0, 0, 0);
	}

	/// Returns a UVec4 mask where the linear components allowed by mAllowedDOFs are 0xffffffff and the rest are 0
	JPH_INLINE UVec4		GetLinearDOFsMask() const
	{
		UVec4 mask(uint32(EAllowedDOFs::TranslationX), uint32(EAllowedDOFs::TranslationY), uint32(EAllowedDOFs::TranslationZ), uint32(EAllowedDOFs::TranslationW));
		return UVec4::sEquals(UVec4::sAnd(UVec4::sReplicate(uint32(mAllowedDOFs)), mask), mask);
	}

	/// Takes a translation vector inV and returns a vector where the components that are not allowed by mAllowedDOFs are set to 0
	JPH_INLINE Vec4			LockTranslation(Vec4Arg inV) const
	{
		return Vec4::sAnd(inV, Vec4(GetLinearDOFsMask().ReinterpretAsFloat()));
	}

	/// Returns a Bivec mask where the angular components allowed by mAllowedDOFs are all-bits-set and the rest are 0.
	/// Component order: [e12(XY), e13(XZ), e14(XW), e23(YZ), e24(YW), e34(ZW)]
	JPH_INLINE Bivec		GetAngularDOFsBivecMask() const
	{
		// TODO(4D): built the same way as GetLinearDOFsMask, in two UVec4 halves recombined into a
		// Vec8. The scalar version this replaces round-tripped eight floats through memory on every
		// call, and this is called twice per MultiplyWorldSpaceInverseInertiaByVector - i.e. twice
		// per constraint row per body per solver iteration.
		UVec4 dofs = UVec4::sReplicate(uint32(mAllowedDOFs));
		UVec4 mask_lo(uint32(EAllowedDOFs::RotationXY), uint32(EAllowedDOFs::RotationXZ), uint32(EAllowedDOFs::RotationXW), uint32(EAllowedDOFs::RotationYZ));
		UVec4 mask_hi(uint32(EAllowedDOFs::RotationYW), uint32(EAllowedDOFs::RotationZW), 0, 0);
		Lane4 lo = UVec4::sEquals(UVec4::sAnd(dofs, mask_lo), mask_lo).ReinterpretAsFloat();

		// The last two lanes are unused by Bivec and must stay zero; sEquals would report 0 == 0 as
		// a match, so mask them off explicitly
		Lane4 hi = UVec4::sAnd(UVec4::sEquals(UVec4::sAnd(dofs, mask_hi), mask_hi), UVec4(0xffffffff, 0xffffffff, 0, 0)).ReinterpretAsFloat();
		return Bivec(Vec8(lo, hi));
	}

	/// Takes an angular velocity / torque bivector inV and returns one where the components that are not allowed by mAllowedDOFs are set to 0
	JPH_INLINE Bivec		LockAngular(BivecArg inV) const
	{
		return Bivec::sAnd(inV, GetAngularDOFsBivecMask());
	}

	/// Used only when this body is dynamic and colliding. Override for the number of solver velocity iterations to run, 0 means use the default in PhysicsSettings::mNumVelocitySteps. The number of iterations to use is the max of all contacts and constraints in the island.
	void					SetNumVelocityStepsOverride(uint inN)							{ JPH_ASSERT(inN < 256); mNumVelocityStepsOverride = uint8(inN); }
	uint					GetNumVelocityStepsOverride() const								{ return mNumVelocityStepsOverride; }

	/// Used only when this body is dynamic and colliding. Override for the number of solver position iterations to run, 0 means use the default in PhysicsSettings::mNumPositionSteps. The number of iterations to use is the max of all contacts and constraints in the island.
	void					SetNumPositionStepsOverride(uint inN)							{ JPH_ASSERT(inN < 256); mNumPositionStepsOverride = uint8(inN); }
	uint					GetNumPositionStepsOverride() const								{ return mNumPositionStepsOverride; }

#ifdef JPH_TRACK_SIMULATION_STATS
	/// Stats for this body. These are average for the simulation island the body was part of.
	struct SimulationStats
	{
		void				Reset()															{ mBroadPhaseTicks = 0; mNarrowPhaseTicks.store(0, memory_order_relaxed); mVelocityConstraintTicks = 0; mPositionConstraintTicks = 0; mUpdateBoundsTicks = 0; mCCDTicks.store(0, memory_order_relaxed); mNumContactConstraints.store(0, memory_order_relaxed); mNumVelocitySteps = 0; mNumPositionSteps = 0; mIsLargeIsland = false; }

		uint64				mBroadPhaseTicks = 0;											///< Number of processor ticks spent doing broad phase collision detection
		atomic<uint64>		mNarrowPhaseTicks = 0;											///< Number of ticks spent doing narrow phase collision detection
		uint64				mVelocityConstraintTicks = 0;									///< Number of ticks spent solving velocity constraints
		uint64				mPositionConstraintTicks = 0;									///< Number of ticks spent solving position constraints
		uint64				mUpdateBoundsTicks = 0;											///< Number of ticks spent updating the broadphase and checking if the body should go to sleep
		atomic<uint64>		mCCDTicks = 0;													///< Number of ticks spent doing CCD
		atomic<uint32>		mNumContactConstraints = 0;										///< Number of contact constraints created for this body
		uint8				mNumVelocitySteps = 0;											///< Number of velocity iterations performed
		uint8				mNumPositionSteps = 0;											///< Number of position iterations performed
		bool				mIsLargeIsland = false;											///< If this body was part of a large island
	};

	const SimulationStats &	GetSimulationStats() const										{ return mSimulationStats; }
	SimulationStats &		GetSimulationStats()											{ return mSimulationStats; }
#endif // JPH_TRACK_SIMULATION_STATS

	////////////////////////////////////////////////////////////
	// FUNCTIONS BELOW THIS LINE ARE FOR INTERNAL USE ONLY
	////////////////////////////////////////////////////////////

	///@name Update linear and angular velocity (used during constraint solving)
	///@{
	inline void				AddLinearVelocityStep(Vec4Arg inLinearVelocityChange)			{ JPH_DET_LOG("AddLinearVelocityStep: " << inLinearVelocityChange); JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite)); mLinearVelocity = LockTranslation(mLinearVelocity + inLinearVelocityChange); JPH_ASSERT(!mLinearVelocity.IsNaN()); }
	inline void				SubLinearVelocityStep(Vec4Arg inLinearVelocityChange)			{ JPH_DET_LOG("SubLinearVelocityStep: " << inLinearVelocityChange); JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite)); mLinearVelocity = LockTranslation(mLinearVelocity - inLinearVelocityChange); JPH_ASSERT(!mLinearVelocity.IsNaN()); }
	inline void				AddAngularVelocityStep(BivecArg inAngularVelocityChange)		{ JPH_DET_LOG("AddAngularVelocityStep: " << inAngularVelocityChange); JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite)); mAngularVelocity += inAngularVelocityChange; JPH_ASSERT(!mAngularVelocity.IsNaN()); }
	inline void				SubAngularVelocityStep(BivecArg inAngularVelocityChange)		{ JPH_DET_LOG("SubAngularVelocityStep: " << inAngularVelocityChange); JPH_ASSERT(BodyAccess::sCheckRights(BodyAccess::sVelocityAccess(), BodyAccess::EAccess::ReadWrite)); mAngularVelocity -= inAngularVelocityChange; JPH_ASSERT(!mAngularVelocity.IsNaN()); }
	///@}

	/// Apply the gyroscopic force (aka Dzhanibekov effect, see https://en.wikipedia.org/wiki/Tennis_racket_theorem)
	inline void				ApplyGyroscopicForceInternal(RotorArg inBodyRotation, float inDeltaTime);

	/// Apply all accumulated forces, torques and drag (should only be called by the PhysicsSystem)
	inline void				ApplyForceTorqueAndDragInternal(RotorArg inBodyRotation, Vec4Arg inGravity, float inDeltaTime);

	/// Access to the island index
	uint32					GetIslandIndexInternal() const									{ return mIslandIndex; }
	void					SetIslandIndexInternal(uint32 inIndex)							{ mIslandIndex = inIndex; }

	/// Access to the index in the active bodies array
	uint32					GetIndexInActiveBodiesInternal() const							{ return mIndexInActiveBodies; }

#ifdef JPH_DOUBLE_PRECISION
	inline DVec4			GetSleepTestOffset() const										{ return DVec4::sLoadDouble3Unsafe(mSleepTestOffset); }
#endif // JPH_DOUBLE_PRECISION

	/// Reset spheres to center around inPoints with radius 0
	inline void				ResetSleepTestSpheres(const RVec4 *inPoints);

	/// Reset the sleep test timer without resetting the sleep test spheres
	inline void				ResetSleepTestTimer()											{ mSleepTestTimer = 0.0f; }

	/// Accumulate sleep time and return if a body can go to sleep
	inline ECanSleep		AccumulateSleepTime(float inDeltaTime, float inTimeBeforeSleep);

	/// Saving state for replay
	void					SaveState(StateRecorder &inStream) const;

	/// Restoring state for replay
	void					RestoreState(StateRecorder &inStream);

	static constexpr uint32	cInactiveIndex = uint32(-1);									///< Constant indicating that body is not active

private:
	friend class BodyManager;
	friend class Body;

	// 16 byte aligned
	Vec4					mLinearVelocity { Vec4::sZero() };								///< World space linear velocity of the center of mass (m/s)
	// 32 byte aligned
	Bivec					mAngularVelocity { Bivec::sZero() };							///< World space angular velocity (bivector, rad/s in each rotation plane)
	Bivec					mInvInertiaDiagonal { Bivec::sZero() };							///< Diagonal of inverse inertia: 6 values for the 6 rotation planes
	Rotor					mInertiaRotation { Rotor::sIdentity() };						///< Rotation (R) that takes inverse inertia diagonal to local space via adjoint action

	// 4 byte aligned
	Float4					mForce { 0, 0, 0, 0 };											///< Accumulated world space force (N)
	Float8					mTorque { 0, 0, 0, 0, 0, 0, 0, 0 };							///< Accumulated world space torque bivector (N m)
	float					mInvMass;														///< Inverse mass of the object (1/kg)
	float					mLinearDamping;													///< Linear damping: dv/dt = -c * v. Value should be zero or positive and is usually close to 0.
	float					mAngularDamping;												///< Angular damping: dw/dt = -c * w. Value should be zero or positive and is usually close to 0.
	float					mMaxLinearVelocity;												///< Maximum linear velocity that this body can reach (m/s)
	float					mMaxAngularVelocity;											///< Maximum angular velocity that this body can reach (rad/s)
	float					mGravityFactor;													///< Factor to multiply gravity with
	uint32					mIndexInActiveBodies = cInactiveIndex;							///< If the body is active, this is the index in the active body list or cInactiveIndex if it is not active
	uint32					mIslandIndex = cInactiveIndex;									///< Index of the island that this body is part of

	// 1 byte aligned
	EMotionQuality			mMotionQuality;													///< Motion quality, or how well it detects collisions when it has a high velocity
	bool					mAllowSleeping;													///< If this body can go to sleep
	EAllowedDOFs			mAllowedDOFs = EAllowedDOFs::All;								///< Allowed degrees of freedom for this body (uint16)
	uint8					mNumVelocityStepsOverride = 0;									///< Override for the number of solver velocity iterations
	uint8					mNumPositionStepsOverride = 0;									///< Override for the number of solver position iterations

	// Sleep test data
#ifdef JPH_DOUBLE_PRECISION
	Double3					mSleepTestOffset;												///< mSleepTestSpheres are relative to this offset to prevent floating point inaccuracies
#endif // JPH_DOUBLE_PRECISION
	Sphere					mSleepTestSpheres[3];											///< Measure motion for 3 points on the body to see if it is resting
	float					mSleepTestTimer;												///< How long this body has been within the movement tolerance

#ifdef JPH_ENABLE_ASSERTS
	EBodyType				mCachedBodyType;												///< Copied from Body::mBodyType and cached for asserting purposes
	EMotionType				mCachedMotionType;												///< Copied from Body::mMotionType and cached for asserting purposes
#endif

#ifdef JPH_TRACK_SIMULATION_STATS
	SimulationStats			mSimulationStats;
#endif // JPH_TRACK_SIMULATION_STATS
};

JPH_NAMESPACE_END

#include "MotionProperties.inl"
