// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Constraints/ContactConstraintManager.h>
#include <Jolt/Physics/Constraints/CalculateSolverSteps.h>
#include <Jolt/Physics/Body/Body.h>
#include <Jolt/Physics/PhysicsUpdateContext.h>
#include <Jolt/Physics/PhysicsSettings.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/Physics/IslandBuilder.h>
#include <Jolt/Physics/DeterminismLog.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Core/QuickSort.h>
#include <Jolt/Core/Prefetch.h>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_DEBUG_RENDERER

JPH_NAMESPACE_BEGIN
#ifdef JPH_WORK_COUNTERS
#include <atomic>
#include <cstdio>
struct WC { std::atomic<unsigned long long> mConstraints{0}, mPoints{0};
  ~WC(){ if(mConstraints.load()) fprintf(stderr,"\nWORK: contact constraints=%llu manifold points=%llu (%.2f pts/constraint)\n",
    (unsigned long long)mConstraints.load(),(unsigned long long)mPoints.load(),double(mPoints.load())/mConstraints.load()); } };
inline WC gWC;
#define JPH_WC(n) { gWC.mConstraints++; gWC.mPoints += (n); }
#else
#define JPH_WC(n)
#endif


using namespace literals;

#ifdef JPH_DEBUG_RENDERER
bool ContactConstraintManager::sDrawContactPoint = false;
bool ContactConstraintManager::sDrawSupportingFaces = false;
bool ContactConstraintManager::sDrawContactPointReduction = false;
bool ContactConstraintManager::sDrawContactManifolds = false;
#endif // JPH_DEBUG_RENDERER

//#define JPH_MANIFOLD_CACHE_DEBUG

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::WorldContactPoint
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::WorldContactPoint::CalculateNonPenetrationConstraintProperties(const Body &inBody1, float inInvMass1, float inInvInertiaScale1, const Body &inBody2, float inInvMass2, float inInvInertiaScale2, RVec4Arg inWorldSpacePosition1, RVec4Arg inWorldSpacePosition2, Vec4Arg inWorldSpaceNormal)
{
	// Calculate collision points relative to body
	RVec4 p = 0.5_r * (inWorldSpacePosition1 + inWorldSpacePosition2);
	Vec4 r1 = Vec4(p - inBody1.GetCenterOfMassPosition());
	Vec4 r2 = Vec4(p - inBody2.GetCenterOfMassPosition());

	mNonPenetrationConstraint.CalculateConstraintPropertiesWithMassOverride(inBody1, inInvMass1, inInvInertiaScale1, r1, inBody2, inInvMass2, inInvInertiaScale2, r2, inWorldSpaceNormal);
}

template <EMotionType Type1, EMotionType Type2>
JPH_INLINE void ContactConstraintManager::WorldContactPoint::TemplatedCalculateFrictionAndNonPenetrationConstraintProperties(float inDeltaTime, Vec4Arg inGravity, const Body &inBody1, const Body &inBody2, float inInvM1, float inInvM2, const InverseInertiaOperator &inInvI1, const InverseInertiaOperator &inInvI2, RVec4Arg inWorldSpacePosition1, RVec4Arg inWorldSpacePosition2, Vec4Arg inWorldSpaceNormal, Vec4Arg inWorldSpaceTangent1, Vec4Arg inWorldSpaceTangent2, Vec4Arg inWorldSpaceTangent3, const ContactSettings &inSettings, float inMinVelocityForRestitution)
{
	JPH_DET_LOG("TemplatedCalculateFrictionAndNonPenetrationConstraintProperties: p1: " << inWorldSpacePosition1 << " p2: " << inWorldSpacePosition2
		<< " normal: " << inWorldSpaceNormal << " tangent1: " << inWorldSpaceTangent1 << " tangent2: " << inWorldSpaceTangent2 << " tangent3: " << inWorldSpaceTangent3
		<< " restitution: " << inSettings.mCombinedRestitution << " friction: " << inSettings.mCombinedFriction << " minv: " << inMinVelocityForRestitution
		<< " surface_vel: " << inSettings.mRelativeLinearSurfaceVelocity << " surface_ang: " << inSettings.mRelativeAngularSurfaceVelocity);

	// Calculate collision points relative to body
	RVec4 p = 0.5_r * (inWorldSpacePosition1 + inWorldSpacePosition2);
	Vec4 r1 = Vec4(p - inBody1.GetCenterOfMassPosition());
	Vec4 r2 = Vec4(p - inBody2.GetCenterOfMassPosition());

	const MotionProperties *mp1 = inBody1.GetMotionPropertiesUnchecked();
	const MotionProperties *mp2 = inBody2.GetMotionPropertiesUnchecked();
	r1.StoreFloat4(&mR1);
	r2.StoreFloat4(&mR2);

	// Calculate velocity of collision points
	Vec4 relative_velocity;
	if constexpr (Type1 != EMotionType::Static && Type2 != EMotionType::Static)
		relative_velocity = mp2->GetPointVelocityCOM(r2) - mp1->GetPointVelocityCOM(r1);
	else if constexpr (Type1 != EMotionType::Static)
		relative_velocity = -mp1->GetPointVelocityCOM(r1);
	else if constexpr (Type2 != EMotionType::Static)
		relative_velocity = mp2->GetPointVelocityCOM(r2);
	else
	{
		JPH_ASSERT(false, "Static vs static makes no sense");
		relative_velocity = Vec4::sZero();
	}
	float normal_velocity = relative_velocity.Dot(inWorldSpaceNormal);

	// How much the shapes are penetrating (> 0 if penetrating, < 0 if separated)
	float penetration = Vec4(inWorldSpacePosition1 - inWorldSpacePosition2).Dot(inWorldSpaceNormal);

	// If there is no penetration, this is a speculative contact and we will apply a bias to the contact constraint
	// so that the constraint becomes relative_velocity . contact normal > -penetration / delta_time
	// instead of relative_velocity . contact normal > 0
	// See: GDC 2013: "Physics for Game Programmers; Continuous Collision" - Erin Catto
	float speculative_contact_velocity_bias = max(0.0f, -penetration / inDeltaTime);

	// Determine if the velocity is big enough for restitution
	float normal_velocity_bias;
	if (inSettings.mCombinedRestitution > 0.0f && normal_velocity < -inMinVelocityForRestitution)
	{
		// We have a velocity that is big enough for restitution. This is where speculative contacts don't work
		// great as we have to decide now if we're going to apply the restitution or not. If the relative
		// velocity is big enough for a hit, we apply the restitution (in the end, due to other constraints,
		// the objects may actually not collide and we will have applied restitution incorrectly). Another
		// artifact that occurs because of this approximation is that the object will bounce from its current
		// position rather than from a position where it is touching the other object. This causes the object
		// to appear to move faster for 1 frame (the opposite of time stealing).
		if (normal_velocity < -speculative_contact_velocity_bias)
		{
			// The gravity / constant forces are applied in the beginning of the time step.
			// If we get here, there was a collision at the beginning of the time step, so we've applied too much force.
			// This means that our calculated restitution can be too high resulting in an increase in energy.
			// So, when we apply restitution, we cancel the added velocity due to these forces.
			Vec4 relative_acceleration;

			// Calculate effect of gravity
			if constexpr (Type1 != EMotionType::Static && Type2 != EMotionType::Static)
				relative_acceleration = inGravity * (mp2->GetGravityFactor() - mp1->GetGravityFactor());
			else if constexpr (Type1 != EMotionType::Static)
				relative_acceleration = -inGravity * mp1->GetGravityFactor();
			else if constexpr (Type2 != EMotionType::Static)
				relative_acceleration = inGravity * mp2->GetGravityFactor();
			else
			{
				JPH_ASSERT(false, "Static vs static makes no sense");
				relative_acceleration = Vec4::sZero();
			}

			// Calculate effect of accumulated forces
			if constexpr (Type1 == EMotionType::Dynamic)
				relative_acceleration -= mp1->GetAccumulatedForce() * mp1->GetInverseMass();
			if constexpr (Type2 == EMotionType::Dynamic)
				relative_acceleration += mp2->GetAccumulatedForce() * mp2->GetInverseMass();

			// We only compensate forces towards the contact normal.
			float force_delta_velocity = min(0.0f, relative_acceleration.Dot(inWorldSpaceNormal) * inDeltaTime);

			normal_velocity_bias = inSettings.mCombinedRestitution * (normal_velocity - force_delta_velocity);
		}
		else
		{
			// In this case we have predicted that we don't hit the other object, but if we do (due to other constraints changing velocities)
			// the speculative contact will prevent penetration but will not apply restitution leading to another artifact.
			normal_velocity_bias = speculative_contact_velocity_bias;
		}
	}
	else
	{
		// No restitution. We can safely apply our contact velocity bias.
		normal_velocity_bias = speculative_contact_velocity_bias;
	}

	mNonPenetrationConstraint.TemplatedCalculateConstraintProperties<Type1, Type2>(inInvM1, inInvI1, r1, inInvM2, inInvI2, r2, inWorldSpaceNormal, normal_velocity_bias);

	// Calculate friction part
	if (inSettings.mCombinedFriction > 0.0f)
	{
		// Get surface velocity relative to tangents. The angular surface velocity is a bivector, so the
		// velocity it induces at r1 is a contraction rather than a cross product.
		Vec4 ws_surface_velocity = inSettings.mRelativeLinearSurfaceVelocity + inSettings.mRelativeAngularSurfaceVelocity.Contract(r1);
		float surface_velocity1 = inWorldSpaceTangent1.Dot(ws_surface_velocity);
		float surface_velocity2 = inWorldSpaceTangent2.Dot(ws_surface_velocity);
		float surface_velocity3 = inWorldSpaceTangent3.Dot(ws_surface_velocity);

		// TODO(4D): implement friction as 3 AxisConstraintParts -- the tangent space of a 4D contact is
		// 3 dimensional, so Coulomb friction needs 3 directions instead of 2.
		mFrictionConstraint1.TemplatedCalculateConstraintProperties<Type1, Type2>(inInvM1, inInvI1, r1, inInvM2, inInvI2, r2, inWorldSpaceTangent1, surface_velocity1);
		mFrictionConstraint2.TemplatedCalculateConstraintProperties<Type1, Type2>(inInvM1, inInvI1, r1, inInvM2, inInvI2, r2, inWorldSpaceTangent2, surface_velocity2);
		mFrictionConstraint3.TemplatedCalculateConstraintProperties<Type1, Type2>(inInvM1, inInvI1, r1, inInvM2, inInvI2, r2, inWorldSpaceTangent3, surface_velocity3);
	}
	else
	{
		// Turn off friction constraint
		mFrictionConstraint1.Deactivate();
		mFrictionConstraint2.Deactivate();
		mFrictionConstraint3.Deactivate();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::ContactConstraint
////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef JPH_DEBUG_RENDERER
void ContactConstraintManager::ContactConstraint::Draw(DebugRenderer *inRenderer, ColorArg inManifoldColor) const
{
	if (mContactPoints.empty())
		return;

	// Get body transforms
	RMat44 transform_body1 = mBody1->GetCenterOfMassTransform();
	RMat44 transform_body2 = mBody2->GetCenterOfMassTransform();

	RVec4 prev_point = transform_body1 * Vec4::sLoadFloat4(&mContactPoints.back().mContactPoint->mPosition1);
	for (const WorldContactPoint &wcp : mContactPoints)
	{
		// Test if any lambda from the previous frame was transferred
		float radius = wcp.mNonPenetrationConstraint.GetTotalLambda() == 0.0f
					&& wcp.mFrictionConstraint1.GetTotalLambda() == 0.0f
					&& wcp.mFrictionConstraint2.GetTotalLambda() == 0.0f? 0.1f :  0.2f;

		RVec4 next_point = transform_body1 * Vec4::sLoadFloat4(&wcp.mContactPoint->mPosition1);
		inRenderer->DrawMarker(next_point, Color::sCyan, radius);
		inRenderer->DrawMarker(transform_body2 * Vec4::sLoadFloat4(&wcp.mContactPoint->mPosition2), Color::sPurple, radius);

		// Draw edge
		inRenderer->DrawArrow(prev_point, next_point, inManifoldColor, 0.05f);
		prev_point = next_point;
	}

	// Draw normal
	RVec4 wp = transform_body1 * Vec4::sLoadFloat4(&mContactPoints[0].mContactPoint->mPosition1);
	inRenderer->DrawArrow(wp, wp + GetWorldSpaceNormal(), Color::sRed, 0.05f);

	// Get tangents (3 of them in 4D)
	Vec4 t1, t2, t3;
	GetTangents(t1, t2, t3);

	// Draw tangents
	inRenderer->DrawLine(wp, wp + t1, Color::sGreen);
	inRenderer->DrawLine(wp, wp + t2, Color::sBlue);
	inRenderer->DrawLine(wp, wp + t3, Color::sYellow);
}
#endif // JPH_DEBUG_RENDERER

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::CachedContactPoint
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::CachedContactPoint::SaveState(StateRecorder &inStream) const
{
	inStream.Write(mPosition1);
	inStream.Write(mPosition2);
	inStream.Write(mNonPenetrationLambda);
	inStream.Write(mFrictionLambda);
}

void ContactConstraintManager::CachedContactPoint::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mPosition1);
	inStream.Read(mPosition2);
	inStream.Read(mNonPenetrationLambda);
	inStream.Read(mFrictionLambda);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::CachedManifold
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::CachedManifold::SaveState(StateRecorder &inStream) const
{
	inStream.Write(mContactNormal);
}

void ContactConstraintManager::CachedManifold::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mContactNormal);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::CachedBodyPair
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::CachedBodyPair::SaveState(StateRecorder &inStream) const
{
	inStream.Write(mDeltaPosition);
	inStream.Write(mDeltaRotation);
}

void ContactConstraintManager::CachedBodyPair::RestoreState(StateRecorder &inStream)
{
	inStream.Read(mDeltaPosition);
	inStream.Read(mDeltaRotation);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager::ManifoldCache
////////////////////////////////////////////////////////////////////////////////////////////////////////

void ContactConstraintManager::ManifoldCache::Init(uint inMaxBodyPairs, uint inMaxContactConstraints, uint inCachedManifoldsSize)
{
	JPH_ASSERT(inMaxContactConstraints <= cMaxContactConstraintsLimit); // Should have been enforced by caller

	uint max_body_pairs = min(inMaxBodyPairs, cMaxBodyPairsLimit);
	JPH_ASSERT(max_body_pairs == inMaxBodyPairs, "Cannot support this many body pairs!");
	max_body_pairs = max(max_body_pairs, 4u); // Because our hash map requires at least 4 buckets, we need to have a minimum number of body pairs

	mAllocator.Init(uint(min(uint64(max_body_pairs) * sizeof(BPKeyValue) + inCachedManifoldsSize, uint64(~uint(0)))));

	mCachedManifolds.Init(GetNextPowerOf2(inMaxContactConstraints));
	mCachedBodyPairs.Init(GetNextPowerOf2(max_body_pairs));
}

void ContactConstraintManager::ManifoldCache::Clear()
{
	JPH_PROFILE_FUNCTION();

	mCachedManifolds.Clear();
	mCachedBodyPairs.Clear();
	mAllocator.Clear();

#ifdef JPH_ENABLE_ASSERTS
	// Mark as incomplete
	mIsFinalized = false;
#endif
}

void ContactConstraintManager::ManifoldCache::Prepare(uint inExpectedNumBodyPairs, uint inExpectedNumManifolds)
{
	// Minimum amount of buckets to use in the hash map
	constexpr uint32 cMinBuckets = 1024;

	// Use the next higher power of 2 of amount of objects in the cache from last frame to determine the amount of buckets in this frame
	mCachedManifolds.SetNumBuckets(min(max(cMinBuckets, GetNextPowerOf2(inExpectedNumManifolds)), mCachedManifolds.GetMaxBuckets()));
	mCachedBodyPairs.SetNumBuckets(min(max(cMinBuckets, GetNextPowerOf2(inExpectedNumBodyPairs)), mCachedBodyPairs.GetMaxBuckets()));
}

const ContactConstraintManager::MKeyValue *ContactConstraintManager::ManifoldCache::Find(const SubShapeIDPair &inKey, uint64 inKeyHash) const
{
	JPH_ASSERT(mIsFinalized);
	return mCachedManifolds.Find(inKey, inKeyHash);
}

ContactConstraintManager::MKeyValue *ContactConstraintManager::ManifoldCache::Create(ContactAllocator &ioContactAllocator, const SubShapeIDPair &inKey, uint64 inKeyHash, int inNumContactPoints)
{
	JPH_ASSERT(!mIsFinalized);
	MKeyValue *kv = mCachedManifolds.Create(ioContactAllocator, inKey, inKeyHash, CachedManifold::sGetRequiredExtraSize(inNumContactPoints));
	if (kv == nullptr)
	{
		ioContactAllocator.mErrors |= EPhysicsUpdateError::ManifoldCacheFull;
		return nullptr;
	}
	kv->GetValue().mNumContactPoints = uint16(inNumContactPoints);
	++ioContactAllocator.mNumManifolds;
	return kv;
}

ContactConstraintManager::MKVAndCreated ContactConstraintManager::ManifoldCache::FindOrCreate(ContactAllocator &ioContactAllocator, const SubShapeIDPair &inKey, uint64 inKeyHash, int inNumContactPoints)
{
	MKeyValue *kv = const_cast<MKeyValue *>(mCachedManifolds.Find(inKey, inKeyHash));
	if (kv != nullptr)
		return { kv, false };

	return { Create(ioContactAllocator, inKey, inKeyHash, inNumContactPoints), true };
}

uint32 ContactConstraintManager::ManifoldCache::ToHandle(const MKeyValue *inKeyValue) const
{
	JPH_ASSERT(!mIsFinalized);
	return mCachedManifolds.ToHandle(inKeyValue);
}

const ContactConstraintManager::MKeyValue *ContactConstraintManager::ManifoldCache::FromHandle(uint32 inHandle) const
{
	JPH_ASSERT(mIsFinalized);
	return mCachedManifolds.FromHandle(inHandle);
}

const ContactConstraintManager::BPKeyValue *ContactConstraintManager::ManifoldCache::Find(const BodyPair &inKey, uint64 inKeyHash) const
{
	JPH_ASSERT(mIsFinalized);
	return mCachedBodyPairs.Find(inKey, inKeyHash);
}

ContactConstraintManager::BPKeyValue *ContactConstraintManager::ManifoldCache::Create(ContactAllocator &ioContactAllocator, const BodyPair &inKey, uint64 inKeyHash)
{
	JPH_ASSERT(!mIsFinalized);
	BPKeyValue *kv = mCachedBodyPairs.Create(ioContactAllocator, inKey, inKeyHash, 0);
	if (kv == nullptr)
	{
		ioContactAllocator.mErrors |= EPhysicsUpdateError::BodyPairCacheFull;
		return nullptr;
	}
	++ioContactAllocator.mNumBodyPairs;
	return kv;
}

void ContactConstraintManager::ManifoldCache::GetAllBodyPairsSorted(Array<const BPKeyValue *> &outAll) const
{
	JPH_ASSERT(mIsFinalized);
	mCachedBodyPairs.GetAllKeyValues(outAll);

	// Sort by key
	QuickSort(outAll.begin(), outAll.end(), [](const BPKeyValue *inLHS, const BPKeyValue *inRHS) {
		return inLHS->GetKey() < inRHS->GetKey();
	});
}

void ContactConstraintManager::ManifoldCache::GetAllManifoldsSorted(const CachedBodyPair &inBodyPair, Array<const MKeyValue *> &outAll) const
{
	JPH_ASSERT(mIsFinalized);

	// Iterate through the attached manifolds
	for (uint32 handle = inBodyPair.mFirstCachedManifold; handle != ManifoldMap::cInvalidHandle; handle = FromHandle(handle)->GetValue().mNextWithSameBodyPair)
	{
		const MKeyValue *kv = mCachedManifolds.FromHandle(handle);
		outAll.push_back(kv);
	}

	// Sort by key
	QuickSort(outAll.begin(), outAll.end(), [](const MKeyValue *inLHS, const MKeyValue *inRHS) {
		return inLHS->GetKey() < inRHS->GetKey();
	});
}

void ContactConstraintManager::ManifoldCache::GetAllCCDManifoldsSorted(Array<const MKeyValue *> &outAll) const
{
	mCachedManifolds.GetAllKeyValues(outAll);

	for (int i = (int)outAll.size() - 1; i >= 0; --i)
		if ((outAll[i]->GetValue().mFlags & (uint16)CachedManifold::EFlags::CCDContact) == 0)
		{
			outAll[i] = outAll.back();
			outAll.pop_back();
		}

	// Sort by key
	QuickSort(outAll.begin(), outAll.end(), [](const MKeyValue *inLHS, const MKeyValue *inRHS) {
		return inLHS->GetKey() < inRHS->GetKey();
	});
}

void ContactConstraintManager::ManifoldCache::ContactPointRemovedCallbacks(ContactListener *inListener)
{
	JPH_PROFILE_FUNCTION();

	for (MKeyValue &kv : mCachedManifolds)
		if ((kv.GetValue().mFlags & uint16(CachedManifold::EFlags::ContactPersisted)) == 0)
			inListener->OnContactRemoved(kv.GetKey());
}

#ifdef JPH_ENABLE_ASSERTS

void ContactConstraintManager::ManifoldCache::Finalize()
{
	mIsFinalized = true;

#ifdef JPH_MANIFOLD_CACHE_DEBUG
	Trace("ManifoldMap:");
	mCachedManifolds.TraceStats();
	Trace("BodyPairMap:");
	mCachedBodyPairs.TraceStats();
#endif // JPH_MANIFOLD_CACHE_DEBUG
}

#endif

void ContactConstraintManager::ManifoldCache::SaveState(StateRecorder &inStream, const StateRecorderFilter *inFilter) const
{
	JPH_ASSERT(mIsFinalized);

	// Get contents of cache
	Array<const BPKeyValue *> all_bp;
	GetAllBodyPairsSorted(all_bp);

	// Determine which ones to save
	Array<const BPKeyValue *> selected_bp;
	if (inFilter == nullptr)
		selected_bp = std::move(all_bp);
	else
	{
		selected_bp.reserve(all_bp.size());
		for (const BPKeyValue *bp_kv : all_bp)
			if (inFilter->ShouldSaveContact(bp_kv->GetKey().mBodyA, bp_kv->GetKey().mBodyB))
				selected_bp.push_back(bp_kv);
	}

	// Write body pairs
	uint32 num_body_pairs = uint32(selected_bp.size());
	inStream.Write(num_body_pairs);
	for (const BPKeyValue *bp_kv : selected_bp)
	{
		// Write body pair key
		inStream.Write(bp_kv->GetKey());

		// Write body pair
		const CachedBodyPair &bp = bp_kv->GetValue();
		bp.SaveState(inStream);

		// Get attached manifolds
		Array<const MKeyValue *> all_m;
		GetAllManifoldsSorted(bp, all_m);

		// Write num manifolds
		uint32 num_manifolds = uint32(all_m.size());
		inStream.Write(num_manifolds);

		// Write all manifolds
		for (const MKeyValue *m_kv : all_m)
		{
			// Write key
			inStream.Write(m_kv->GetKey());
			const CachedManifold &cm = m_kv->GetValue();
			JPH_ASSERT((cm.mFlags & (uint16)CachedManifold::EFlags::CCDContact) == 0);

			// Write amount of contacts
			inStream.Write(cm.mNumContactPoints);

			// Write manifold
			cm.SaveState(inStream);

			// Write contact points
			for (uint32 i = 0; i < cm.mNumContactPoints; ++i)
				cm.mContactPoints[i].SaveState(inStream);
		}
	}

	// Get CCD manifolds
	Array<const MKeyValue *> all_m;
	GetAllCCDManifoldsSorted(all_m);

	// Determine which ones to save
	Array<const MKeyValue *> selected_m;
	if (inFilter == nullptr)
		selected_m = std::move(all_m);
	else
	{
		selected_m.reserve(all_m.size());
		for (const MKeyValue *m_kv : all_m)
			if (inFilter->ShouldSaveContact(m_kv->GetKey().GetBody1ID(), m_kv->GetKey().GetBody2ID()))
				selected_m.push_back(m_kv);
	}

	// Write all CCD manifold keys
	uint32 num_manifolds = uint32(selected_m.size());
	inStream.Write(num_manifolds);
	for (const MKeyValue *m_kv : selected_m)
		inStream.Write(m_kv->GetKey());
}

bool ContactConstraintManager::ManifoldCache::RestoreState(const ManifoldCache &inReadCache, StateRecorder &inStream, const StateRecorderFilter *inFilter)
{
	JPH_ASSERT(!mIsFinalized);

	bool success = true;

	// Create a contact allocator for restoring the contact cache
	ContactAllocator contact_allocator(GetContactAllocator());

	// When validating, get all existing body pairs
	Array<const BPKeyValue *> all_bp;
	if (inStream.IsValidating())
		inReadCache.GetAllBodyPairsSorted(all_bp);

	// Read amount of body pairs
	uint32 num_body_pairs;
	if (inStream.IsValidating())
		num_body_pairs = uint32(all_bp.size());
	inStream.Read(num_body_pairs);

	// Read entire cache
	for (uint32 i = 0; i < num_body_pairs; ++i)
	{
		// Read key
		BodyPair body_pair_key;
		if (inStream.IsValidating() && i < all_bp.size())
			body_pair_key = all_bp[i]->GetKey();
		inStream.Read(body_pair_key);

		// Check if we want to restore this contact
		if (inFilter == nullptr || inFilter->ShouldRestoreContact(body_pair_key.mBodyA, body_pair_key.mBodyB))
		{
			// Create new entry for this body pair
			uint64 body_pair_hash = body_pair_key.GetHash();
			BPKeyValue *bp_kv = Create(contact_allocator, body_pair_key, body_pair_hash);
			if (bp_kv == nullptr)
			{
				// Out of cache space
				success = false;
				break;
			}
			CachedBodyPair &bp = bp_kv->GetValue();

			// Read body pair
			if (inStream.IsValidating() && i < all_bp.size())
				memcpy(&bp, &all_bp[i]->GetValue(), sizeof(CachedBodyPair));
			bp.RestoreState(inStream);

			// When validating, get all existing manifolds
			Array<const MKeyValue *> all_m;
			if (inStream.IsValidating())
				inReadCache.GetAllManifoldsSorted(all_bp[i]->GetValue(), all_m);

			// Read amount of manifolds
			uint32 num_manifolds = 0;
			if (inStream.IsValidating())
				num_manifolds = uint32(all_m.size());
			inStream.Read(num_manifolds);

			uint32 handle = ManifoldMap::cInvalidHandle;
			for (uint32 j = 0; j < num_manifolds; ++j)
			{
				// Read key
				SubShapeIDPair sub_shape_key;
				if (inStream.IsValidating() && j < all_m.size())
					sub_shape_key = all_m[j]->GetKey();
				inStream.Read(sub_shape_key);
				uint64 sub_shape_key_hash = sub_shape_key.GetHash();

				// Read amount of contact points
				uint16 num_contact_points = 0;
				if (inStream.IsValidating() && j < all_m.size())
					num_contact_points = all_m[j]->GetValue().mNumContactPoints;
				inStream.Read(num_contact_points);

				// Read manifold
				MKeyValue *m_kv = Create(contact_allocator, sub_shape_key, sub_shape_key_hash, num_contact_points);
				if (m_kv == nullptr)
				{
					// Out of cache space
					success = false;
					break;
				}
				CachedManifold &cm = m_kv->GetValue();
				if (inStream.IsValidating() && j < all_m.size())
				{
					memcpy(&cm, &all_m[j]->GetValue(), CachedManifold::sGetRequiredTotalSize(num_contact_points));
					cm.mNumContactPoints = uint16(num_contact_points); // Restore num contact points
				}
				cm.RestoreState(inStream);
				cm.mNextWithSameBodyPair = handle;
				handle = ToHandle(m_kv);

				// Read contact points
				for (uint32 k = 0; k < num_contact_points; ++k)
					cm.mContactPoints[k].RestoreState(inStream);
			}
			bp.mFirstCachedManifold = handle;
		}
		else
		{
			// Skip the contact
			CachedBodyPair bp;
			bp.RestoreState(inStream);
			uint32 num_manifolds = 0;
			inStream.Read(num_manifolds);
			for (uint32 j = 0; j < num_manifolds; ++j)
			{
				SubShapeIDPair sub_shape_key;
				inStream.Read(sub_shape_key);
				uint16 num_contact_points;
				inStream.Read(num_contact_points);
				CachedManifold cm;
				cm.RestoreState(inStream);
				for (uint32 k = 0; k < num_contact_points; ++k)
					cm.mContactPoints[0].RestoreState(inStream);
			}
		}
	}

	// When validating, get all existing CCD manifolds
	Array<const MKeyValue *> all_m;
	if (inStream.IsValidating())
		inReadCache.GetAllCCDManifoldsSorted(all_m);

	// Read amount of CCD manifolds
	uint32 num_manifolds;
	if (inStream.IsValidating())
		num_manifolds = uint32(all_m.size());
	inStream.Read(num_manifolds);

	for (uint32 j = 0; j < num_manifolds; ++j)
	{
		// Read key
		SubShapeIDPair sub_shape_key;
		if (inStream.IsValidating() && j < all_m.size())
			sub_shape_key = all_m[j]->GetKey();
		inStream.Read(sub_shape_key);

		// Check if we want to restore this contact
		if (inFilter == nullptr || inFilter->ShouldRestoreContact(sub_shape_key.GetBody1ID(), sub_shape_key.GetBody2ID()))
		{
			// Create CCD manifold
			uint64 sub_shape_key_hash = sub_shape_key.GetHash();
			MKeyValue *m_kv = Create(contact_allocator, sub_shape_key, sub_shape_key_hash, 0);
			if (m_kv == nullptr)
			{
				// Out of cache space
				success = false;
				break;
			}
			CachedManifold &cm = m_kv->GetValue();
			cm.mFlags |= (uint16)CachedManifold::EFlags::CCDContact;
		}
	}

#ifdef JPH_ENABLE_ASSERTS
	// We don't finalize until the last part is restored
	if (inStream.IsLastPart())
		mIsFinalized = true;
#endif

	return success;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContactConstraintManager
////////////////////////////////////////////////////////////////////////////////////////////////////////

ContactConstraintManager::ContactConstraintManager(const PhysicsSettings &inPhysicsSettings) :
	mPhysicsSettings(inPhysicsSettings)
{
#ifdef JPH_ENABLE_ASSERTS
	// For the first frame mark this empty buffer as finalized
	mCache[mCacheWriteIdx ^ 1].Finalize();
#endif
}

ContactConstraintManager::~ContactConstraintManager()
{
	JPH_ASSERT(mConstraints == nullptr);
}

void ContactConstraintManager::Init(uint inMaxBodyPairs, uint inMaxContactConstraints)
{
	// Limit the number of constraints so that the allocation size fits in an unsigned integer
	mMaxConstraints = min(inMaxContactConstraints, cMaxContactConstraintsLimit);
	JPH_ASSERT(mMaxConstraints == inMaxContactConstraints, "Cannot support this many contact constraints!");
	mMaxConstraints = max(mMaxConstraints, 4u); // Because our hash map requires at least 4 buckets, we need to have a minimum number of constraints

	// Calculate worst case cache usage
	constexpr uint cMaxManifoldSizePerConstraint = sizeof(MKeyValue) + CachedManifold::sGetRequiredExtraSize(MaxContactPoints);
	static_assert(cMaxManifoldSizePerConstraint < sizeof(ContactConstraint)); // If not true, then the next line can overflow
	uint cached_manifolds_size = mMaxConstraints * cMaxManifoldSizePerConstraint;

	// Init the caches
	mCache[0].Init(inMaxBodyPairs, mMaxConstraints, cached_manifolds_size);
	mCache[1].Init(inMaxBodyPairs, mMaxConstraints, cached_manifolds_size);
}

void ContactConstraintManager::PrepareConstraintBuffer(PhysicsUpdateContext *inContext)
{
	// Store context
	mUpdateContext = inContext;

	// Allocate temporary constraint buffer
	JPH_ASSERT(mConstraints == nullptr);
	mConstraints = (ContactConstraint *)inContext->mTempAllocator->Allocate(mMaxConstraints * sizeof(ContactConstraint));
}

template <EMotionType Type1, EMotionType Type2>
JPH_INLINE void ContactConstraintManager::TemplatedCalculateFrictionAndNonPenetrationConstraintProperties(ContactConstraint &ioConstraint, const ContactSettings &inSettings, float inDeltaTime, Vec4Arg inGravity, RMat44Arg inTransformBody1, RMat44Arg inTransformBody2, const Body &inBody1, const Body &inBody2)
{
	// Calculate scaled mass and the world space inverse inertia. TODO(4D): the inverse inertia is a
	// 6x6 operator on bivectors rather than upstream's Mat44, but like upstream's it is built once
	// here and reused by all 24 rows of the contact - see InverseInertiaOperator.
	InverseInertiaOperator inv_i1 = Type1 == EMotionType::Dynamic? inBody1.GetMotionPropertiesUnchecked()->GetInverseInertiaForRotation(inBody1.GetRotation(), inSettings.mInvInertiaScale1) : InverseInertiaOperator::sZero();
	if constexpr (Type1 == EMotionType::Dynamic)
		ioConstraint.mInvMass1 = inSettings.mInvMassScale1 * inBody1.GetMotionPropertiesUnchecked()->GetInverseMass();
	else
		ioConstraint.mInvMass1 = 0.0f;

	InverseInertiaOperator inv_i2 = Type2 == EMotionType::Dynamic? inBody2.GetMotionPropertiesUnchecked()->GetInverseInertiaForRotation(inBody2.GetRotation(), inSettings.mInvInertiaScale2) : InverseInertiaOperator::sZero();
	if constexpr (Type2 == EMotionType::Dynamic)
		ioConstraint.mInvMass2 = inSettings.mInvMassScale2 * inBody2.GetMotionPropertiesUnchecked()->GetInverseMass();
	else
		ioConstraint.mInvMass2 = 0.0f;

	// Calculate tangents (3 of them in 4D)
	Vec4 t1, t2, t3;
	ioConstraint.GetTangents(t1, t2, t3);

	Vec4 ws_normal = ioConstraint.GetWorldSpaceNormal();

	// Setup velocity constraint properties
	float min_velocity_for_restitution = mPhysicsSettings.mMinVelocityForRestitution;
	for (WorldContactPoint &wcp : ioConstraint.mContactPoints)
	{
		RVec4 p1 = inTransformBody1 * Vec4::sLoadFloat4(&wcp.mContactPoint->mPosition1);
		RVec4 p2 = inTransformBody2 * Vec4::sLoadFloat4(&wcp.mContactPoint->mPosition2);
		wcp.TemplatedCalculateFrictionAndNonPenetrationConstraintProperties<Type1, Type2>(inDeltaTime, inGravity, inBody1, inBody2, ioConstraint.mInvMass1, ioConstraint.mInvMass2, inv_i1, inv_i2, p1, p2, ws_normal, t1, t2, t3, inSettings, min_velocity_for_restitution);
	}
}

inline void ContactConstraintManager::CalculateFrictionAndNonPenetrationConstraintProperties(ContactConstraint &ioConstraint, const ContactSettings &inSettings, float inDeltaTime, Vec4Arg inGravity, RMat44Arg inTransformBody1, RMat44Arg inTransformBody2, const Body &inBody1, const Body &inBody2)
{
	// Dispatch to the correct templated form
	switch (inBody1.GetMotionType())
	{
	case EMotionType::Dynamic:
		switch (inBody2.GetMotionType())
		{
		case EMotionType::Dynamic:
			TemplatedCalculateFrictionAndNonPenetrationConstraintProperties<EMotionType::Dynamic, EMotionType::Dynamic>(ioConstraint, inSettings, inDeltaTime, inGravity, inTransformBody1, inTransformBody2, inBody1, inBody2);
			break;

		case EMotionType::Kinematic:
			TemplatedCalculateFrictionAndNonPenetrationConstraintProperties<EMotionType::Dynamic, EMotionType::Kinematic>(ioConstraint, inSettings, inDeltaTime, inGravity, inTransformBody1, inTransformBody2, inBody1, inBody2);
			break;

		case EMotionType::Static:
			TemplatedCalculateFrictionAndNonPenetrationConstraintProperties<EMotionType::Dynamic, EMotionType::Static>(ioConstraint, inSettings, inDeltaTime, inGravity, inTransformBody1, inTransformBody2, inBody1, inBody2);
			break;

		default:
			JPH_ASSERT(false);
			break;
		}
		break;

	case EMotionType::Kinematic:
		JPH_ASSERT(inBody2.IsDynamic());
		TemplatedCalculateFrictionAndNonPenetrationConstraintProperties<EMotionType::Kinematic, EMotionType::Dynamic>(ioConstraint, inSettings, inDeltaTime, inGravity, inTransformBody1, inTransformBody2, inBody1, inBody2);
		break;

	case EMotionType::Static:
		JPH_ASSERT(inBody2.IsDynamic());
		TemplatedCalculateFrictionAndNonPenetrationConstraintProperties<EMotionType::Static, EMotionType::Dynamic>(ioConstraint, inSettings, inDeltaTime, inGravity, inTransformBody1, inTransformBody2, inBody1, inBody2);
		break;

	default:
		JPH_ASSERT(false);
		break;
	}
}

void ContactConstraintManager::GetContactsFromCache(ContactAllocator &ioContactAllocator, Body &inBody1, Body &inBody2, bool &outPairHandled, bool &outConstraintCreated)
{
	// Start with nothing found and not handled
	outConstraintCreated = false;
	outPairHandled = false;

	// Swap bodies so that body 1 id < body 2 id
	Body *body1, *body2;
	if (inBody1.GetID() < inBody2.GetID())
	{
		body1 = &inBody1;
		body2 = &inBody2;
	}
	else
	{
		body1 = &inBody2;
		body2 = &inBody1;
	}

	// Find the cached body pair
	BodyPair body_pair_key(body1->GetID(), body2->GetID());
	uint64 body_pair_hash = body_pair_key.GetHash();
	const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
	const BPKeyValue *kv = read_cache.Find(body_pair_key, body_pair_hash);
	if (kv == nullptr)
		return;
	const CachedBodyPair &input_cbp = kv->GetValue();

	// Get relative translation
	Rotor inv_r1 = body1->GetRotation().Reversed();
	Vec4 delta_position = inv_r1 * Vec4(body2->GetCenterOfMassPosition() - body1->GetCenterOfMassPosition());

	// Get old position delta
	Vec4 old_delta_position = Vec4::sLoadFloat4(&input_cbp.mDeltaPosition);

	// Check if bodies are still roughly in the same relative position
	if ((delta_position - old_delta_position).LengthSq() > mPhysicsSettings.mBodyPairCacheMaxDeltaPositionSq)
		return;

	// Determine relative orientation
	Rotor delta_rotation = inv_r1 * body2->GetRotation();

	// Load the old rotor delta (stored uncompressed, see CachedBodyPair::mDeltaRotation)
	Rotor old_delta_rotation(input_cbp.mDeltaRotation);

	// Check if bodies are still roughly in the same relative orientation.
	// TODO(4D): the 3D test is |p . q| < cos(angle / 2), using that the scalar part of the relative
	// quaternion is cos(angle / 2). The 8 component rotor dot product is likewise the scalar part of
	// p ~q, which under the SU(2)xSU(2) split equals (cos(a/2) + cos(b/2)) / 2 for a rotation with
	// angles a and b. It is still 1 at zero rotation and decreases as either angle grows, so the same
	// comparison remains a valid (and for double rotations, conservative) small-rotation test.
	if (abs(delta_rotation.Dot(old_delta_rotation)) < mPhysicsSettings.mBodyPairCacheCosMaxDeltaRotationDiv2)
		return;

	// The cache is valid, return that we've handled this body pair
	outPairHandled = true;

	// Copy the cached body pair to this frame
	ManifoldCache &write_cache = mCache[mCacheWriteIdx];
	BPKeyValue *output_bp_kv = write_cache.Create(ioContactAllocator, body_pair_key, body_pair_hash);
	if (output_bp_kv == nullptr)
		return; // Out of cache space
	CachedBodyPair *output_cbp = &output_bp_kv->GetValue();
	memcpy(output_cbp, &input_cbp, sizeof(CachedBodyPair));

	// If there were no contacts, we have handled the contact
	if (input_cbp.mFirstCachedManifold == ManifoldMap::cInvalidHandle)
		return;

	// Get body transforms
	RMat44 transform_body1 = body1->GetCenterOfMassTransform();
	RMat44 transform_body2 = body2->GetCenterOfMassTransform();

	// Get time step and gravity
	float delta_time = mUpdateContext->mStepDeltaTime;
	Vec4 gravity = mUpdateContext->mPhysicsSystem->GetGravity();

	// Copy manifolds
	uint32 output_handle = ManifoldMap::cInvalidHandle;
	uint32 input_handle = input_cbp.mFirstCachedManifold;
	do
	{
		JPH_PROFILE("Add Constraint From Cached Manifold");

		// Find the existing manifold
		const MKeyValue *input_kv = read_cache.FromHandle(input_handle);
		const SubShapeIDPair &input_key = input_kv->GetKey();
		const CachedManifold &input_cm = input_kv->GetValue();
		JPH_ASSERT(input_cm.mNumContactPoints > 0); // There should be contact points in this manifold!

		// Create room for manifold in write buffer and copy data
		uint64 input_hash = input_key.GetHash();
		MKeyValue *output_kv = write_cache.Create(ioContactAllocator, input_key, input_hash, input_cm.mNumContactPoints);
		if (output_kv == nullptr)
			break; // Out of cache space
		CachedManifold *output_cm = &output_kv->GetValue();
		memcpy(output_cm, &input_cm, CachedManifold::sGetRequiredTotalSize(input_cm.mNumContactPoints));

		// Link the object under the body pairs
		output_cm->mNextWithSameBodyPair = output_handle;
		output_handle = write_cache.ToHandle(output_kv);

		// Calculate default contact settings
		ContactSettings settings;
		settings.mCombinedFriction = mCombineFriction(*body1, input_key.GetSubShapeID1(), *body2, input_key.GetSubShapeID2());
		settings.mCombinedRestitution = mCombineRestitution(*body1, input_key.GetSubShapeID1(), *body2, input_key.GetSubShapeID2());
		settings.mIsSensor = body1->IsSensor() || body2->IsSensor();

		// Calculate world space contact normal
		Vec4 world_space_normal = transform_body2.Multiply3x3(Vec4::sLoadFloat4(&output_cm->mContactNormal)).Normalized();

		// Call contact listener to update settings
		if (mContactListener != nullptr)
		{
			// Convert constraint to manifold structure for callback
			ContactManifold manifold;
			manifold.mWorldSpaceNormal = world_space_normal;
			manifold.mSubShapeID1 = input_key.GetSubShapeID1();
			manifold.mSubShapeID2 = input_key.GetSubShapeID2();
			manifold.mBaseOffset = transform_body1.GetTranslation();
			manifold.mRelativeContactPointsOn1.resize(output_cm->mNumContactPoints);
			manifold.mRelativeContactPointsOn2.resize(output_cm->mNumContactPoints);
			Mat44 local_transform_body2 = transform_body2.PostTranslated(-manifold.mBaseOffset).ToMat44();
			float penetration_depth = -FLT_MAX;
			for (uint32 i = 0; i < output_cm->mNumContactPoints; ++i)
			{
				const CachedContactPoint &ccp = output_cm->mContactPoints[i];
				manifold.mRelativeContactPointsOn1[i] = transform_body1.Multiply3x3(Vec4::sLoadFloat4(&ccp.mPosition1));
				manifold.mRelativeContactPointsOn2[i] = local_transform_body2 * Vec4::sLoadFloat4(&ccp.mPosition2);
				penetration_depth = max(penetration_depth, (manifold.mRelativeContactPointsOn1[i] - manifold.mRelativeContactPointsOn2[i]).Dot(world_space_normal));
			}
			manifold.mPenetrationDepth = penetration_depth; // We don't have the penetration depth anymore, estimate it

			// Notify callback
			mContactListener->OnContactPersisted(*body1, *body2, manifold, settings);
		}

		JPH_ASSERT(settings.mIsSensor || !(body1->IsSensor() || body2->IsSensor()), "Sensors cannot be converted into regular bodies by a contact callback!");
		if (!settings.mIsSensor // If one of the bodies is a sensor, don't actually create the constraint
			&& ((body1->IsDynamic() && settings.mInvMassScale1 != 0.0f) // One of the bodies must have mass to be able to create a contact constraint
				|| (body2->IsDynamic() && settings.mInvMassScale2 != 0.0f)))
		{
			// Add contact constraint in world space for the solver
			uint32 constraint_idx = mNumConstraints++;
			if (constraint_idx >= mMaxConstraints)
			{
				ioContactAllocator.mErrors |= EPhysicsUpdateError::ContactConstraintsFull;
				break;
			}

			// A constraint will be created
			outConstraintCreated = true;

			ContactConstraint &constraint = mConstraints[constraint_idx];
			new (&constraint) ContactConstraint();
			constraint.mBody1 = body1;
			constraint.mBody2 = body2;
			constraint.mSortKey = input_hash;
			constraint.SetWorldSpaceNormal(world_space_normal);
			constraint.mCombinedFriction = settings.mCombinedFriction;
			constraint.mInvInertiaScale1 = settings.mInvInertiaScale1;
			constraint.mInvInertiaScale2 = settings.mInvInertiaScale2;
			constraint.mContactPoints.resize(output_cm->mNumContactPoints);
			JPH_WC(output_cm->mNumContactPoints)
			for (uint32 i = 0; i < output_cm->mNumContactPoints; ++i)
			{
				CachedContactPoint &ccp = output_cm->mContactPoints[i];
				WorldContactPoint &wcp = constraint.mContactPoints[i];
				wcp.mNonPenetrationConstraint.SetTotalLambda(ccp.mNonPenetrationLambda);
				wcp.mFrictionConstraint1.SetTotalLambda(ccp.mFrictionLambda[0]);
				wcp.mFrictionConstraint2.SetTotalLambda(ccp.mFrictionLambda[1]);
				wcp.mFrictionConstraint3.SetTotalLambda(ccp.mFrictionLambda[2]);
				wcp.mContactPoint = &ccp;
			}

			JPH_DET_LOG("GetContactsFromCache: id1: " << constraint.mBody1->GetID() << " id2: " << constraint.mBody2->GetID() << " key: " << constraint.mSortKey);

			// Calculate friction and non-penetration constraint properties for all contact points
			CalculateFrictionAndNonPenetrationConstraintProperties(constraint, settings, delta_time, gravity, transform_body1, transform_body2, *body1, *body2);

			// Notify island builder
			mUpdateContext->mIslandBuilder->LinkContact(constraint_idx, body1->GetIndexInActiveBodiesInternal(), body2->GetIndexInActiveBodiesInternal());

		#ifdef JPH_DEBUG_RENDERER
			// Draw the manifold
			if (sDrawContactManifolds)
				constraint.Draw(DebugRenderer::sInstance, Color::sYellow);
		#endif // JPH_DEBUG_RENDERER

		#ifdef JPH_TRACK_SIMULATION_STATS
			// Track new contact constraints
			if (!body1->IsStatic())
				body1->GetMotionPropertiesUnchecked()->GetSimulationStats().mNumContactConstraints.fetch_add(1, memory_order_relaxed);
			if (!body2->IsStatic())
				body2->GetMotionPropertiesUnchecked()->GetSimulationStats().mNumContactConstraints.fetch_add(1, memory_order_relaxed);
		#endif
		}

		// Mark contact as persisted so that we won't fire OnContactRemoved callbacks
		input_cm.mFlags |= (uint16)CachedManifold::EFlags::ContactPersisted;

		// Fetch the next manifold
		input_handle = input_cm.mNextWithSameBodyPair;
	}
	while (input_handle != ManifoldMap::cInvalidHandle);
	output_cbp->mFirstCachedManifold = output_handle;
}

ContactConstraintManager::BodyPairHandle ContactConstraintManager::AddBodyPair(ContactAllocator &ioContactAllocator, const Body &inBody1, const Body &inBody2)
{
	// Swap bodies so that body 1 id < body 2 id
	const Body *body1, *body2;
	if (inBody1.GetID() < inBody2.GetID())
	{
		body1 = &inBody1;
		body2 = &inBody2;
	}
	else
	{
		body1 = &inBody2;
		body2 = &inBody1;
	}

	// Add an entry
	BodyPair body_pair_key(body1->GetID(), body2->GetID());
	uint64 body_pair_hash = body_pair_key.GetHash();
	BPKeyValue *body_pair_kv = mCache[mCacheWriteIdx].Create(ioContactAllocator, body_pair_key, body_pair_hash);
	if (body_pair_kv == nullptr)
		return nullptr; // Out of cache space
	CachedBodyPair *cbp = &body_pair_kv->GetValue();
	cbp->mFirstCachedManifold = ManifoldMap::cInvalidHandle;

	// Get relative translation
	Rotor inv_r1 = body1->GetRotation().Reversed();
	Vec4 delta_position = inv_r1 * Vec4(body2->GetCenterOfMassPosition() - body1->GetCenterOfMassPosition());

	// Store it
	delta_position.StoreFloat4(&cbp->mDeltaPosition);

	// Determine relative orientation
	Rotor delta_rotation = inv_r1 * body2->GetRotation();

	// Store it
	delta_rotation.StoreFloat8(&cbp->mDeltaRotation);

	return cbp;
}

template <EMotionType Type1, EMotionType Type2>
bool ContactConstraintManager::TemplatedAddContactConstraint(ContactAllocator &ioContactAllocator, BodyPairHandle inBodyPairHandle, Body &inBody1, Body &inBody2, const ContactManifold &inManifold)
{
	// Calculate hash
	SubShapeIDPair key { inBody1.GetID(), inManifold.mSubShapeID1, inBody2.GetID(), inManifold.mSubShapeID2 };
	uint64 key_hash = key.GetHash();

	// Determine number of contact points
	int num_contact_points = (int)inManifold.mRelativeContactPointsOn1.size();
	JPH_ASSERT(num_contact_points <= MaxContactPoints);
	JPH_ASSERT(num_contact_points == (int)inManifold.mRelativeContactPointsOn2.size());

	// Reserve space for new contact cache entry
	// Note that for dynamic vs dynamic we always require the first body to have a lower body id to get a consistent key
	// under which to look up the contact
	ManifoldCache &write_cache = mCache[mCacheWriteIdx];
	MKeyValue *new_manifold_kv = write_cache.Create(ioContactAllocator, key, key_hash, num_contact_points);
	if (new_manifold_kv == nullptr)
		return false; // Out of cache space
	CachedManifold *new_manifold = &new_manifold_kv->GetValue();

	// Transform the world space normal to the space of body 2 (this is usually the static body)
	RMat44 inverse_transform_body2 = inBody2.GetInverseCenterOfMassTransform();
	inverse_transform_body2.Multiply3x3(inManifold.mWorldSpaceNormal).Normalized().StoreFloat4(&new_manifold->mContactNormal);

	// Settings object that gets passed to the callback
	ContactSettings settings;
	settings.mCombinedFriction = mCombineFriction(inBody1, inManifold.mSubShapeID1, inBody2, inManifold.mSubShapeID2);
	settings.mCombinedRestitution = mCombineRestitution(inBody1, inManifold.mSubShapeID1, inBody2, inManifold.mSubShapeID2);
	settings.mIsSensor = inBody1.IsSensor() || inBody2.IsSensor();

	// Get the contact points for the old cache entry
	const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
	const MKeyValue *old_manifold_kv = read_cache.Find(key, key_hash);
	const CachedContactPoint *ccp_start;
	const CachedContactPoint *ccp_end;
	if (old_manifold_kv != nullptr)
	{
		// Call point persisted listener
		if (mContactListener != nullptr)
			mContactListener->OnContactPersisted(inBody1, inBody2, inManifold, settings);

		// Fetch the contact points from the old manifold
		const CachedManifold *old_manifold = &old_manifold_kv->GetValue();
		ccp_start = old_manifold->mContactPoints;
		ccp_end = ccp_start + old_manifold->mNumContactPoints;

		// Mark contact as persisted so that we won't fire OnContactRemoved callbacks
		old_manifold->mFlags |= (uint16)CachedManifold::EFlags::ContactPersisted;
	}
	else
	{
		// Call point added listener
		if (mContactListener != nullptr)
			mContactListener->OnContactAdded(inBody1, inBody2, inManifold, settings);

		// No contact points available from old manifold
		ccp_start = nullptr;
		ccp_end = nullptr;
	}

	// Get inverse transform for body 1
	RMat44 inverse_transform_body1 = inBody1.GetInverseCenterOfMassTransform();

	bool contact_constraint_created = false;

	// If one of the bodies is a sensor, don't actually create the constraint
	JPH_ASSERT(settings.mIsSensor || !(inBody1.IsSensor() || inBody2.IsSensor()), "Sensors cannot be converted into regular bodies by a contact callback!");
	if (!settings.mIsSensor
		&& ((inBody1.IsDynamic() && settings.mInvMassScale1 != 0.0f) // One of the bodies must have mass to be able to create a contact constraint
			|| (inBody2.IsDynamic() && settings.mInvMassScale2 != 0.0f)))
	{
		// Add contact constraint
		uint32 constraint_idx = mNumConstraints++;
		if (constraint_idx >= mMaxConstraints)
		{
			ioContactAllocator.mErrors |= EPhysicsUpdateError::ContactConstraintsFull;

			// Manifold has been created already, we're not filling it in, so we need to reset the contact number of points.
			// Note that we don't hook it up to the body pair cache so that it won't be used as a cache during the next simulation.
			new_manifold->mNumContactPoints = 0;
			return false;
		}

		// We will create a contact constraint
		contact_constraint_created = true;

		ContactConstraint &constraint = mConstraints[constraint_idx];
		new (&constraint) ContactConstraint();
		constraint.mBody1 = &inBody1;
		constraint.mBody2 = &inBody2;
		constraint.mSortKey = key_hash;
		constraint.SetWorldSpaceNormal(inManifold.mWorldSpaceNormal);
		constraint.mCombinedFriction = settings.mCombinedFriction;
		constraint.mInvInertiaScale1 = settings.mInvInertiaScale1;
		constraint.mInvInertiaScale2 = settings.mInvInertiaScale2;

		JPH_DET_LOG("TemplatedAddContactConstraint: id1: " << constraint.mBody1->GetID() << " id2: " << constraint.mBody2->GetID() << " key: " << constraint.mSortKey);

		// Notify island builder
		mUpdateContext->mIslandBuilder->LinkContact(constraint_idx, inBody1.GetIndexInActiveBodiesInternal(), inBody2.GetIndexInActiveBodiesInternal());

		// Get time step and gravity
		float delta_time = mUpdateContext->mStepDeltaTime;
		Vec4 gravity = mUpdateContext->mPhysicsSystem->GetGravity();

		// Calculate scaled mass and the world space inverse inertia; see the note in
		// TemplatedCalculateFrictionAndNonPenetrationConstraintProperties
		InverseInertiaOperator inv_i1 = Type1 == EMotionType::Dynamic? inBody1.GetMotionPropertiesUnchecked()->GetInverseInertiaForRotation(inBody1.GetRotation(), settings.mInvInertiaScale1) : InverseInertiaOperator::sZero();
		if constexpr (Type1 == EMotionType::Dynamic)
			constraint.mInvMass1 = settings.mInvMassScale1 * inBody1.GetMotionPropertiesUnchecked()->GetInverseMass();
		else
			constraint.mInvMass1 = 0.0f;

		InverseInertiaOperator inv_i2 = Type2 == EMotionType::Dynamic? inBody2.GetMotionPropertiesUnchecked()->GetInverseInertiaForRotation(inBody2.GetRotation(), settings.mInvInertiaScale2) : InverseInertiaOperator::sZero();
		if constexpr (Type2 == EMotionType::Dynamic)
			constraint.mInvMass2 = settings.mInvMassScale2 * inBody2.GetMotionPropertiesUnchecked()->GetInverseMass();
		else
			constraint.mInvMass2 = 0.0f;

		// Calculate tangents (3 of them in 4D)
		Vec4 t1, t2, t3;
		constraint.GetTangents(t1, t2, t3);

		constraint.mContactPoints.resize(num_contact_points);
		JPH_WC(num_contact_points)
		for (int i = 0; i < num_contact_points; ++i)
		{
			// Convert to world space and set positions
			WorldContactPoint &wcp = constraint.mContactPoints[i];
			RVec4 p1_ws = inManifold.mBaseOffset + inManifold.mRelativeContactPointsOn1[i];
			RVec4 p2_ws = inManifold.mBaseOffset + inManifold.mRelativeContactPointsOn2[i];

			// Convert to local space to the body
			Vec4 p1_ls = Vec4(inverse_transform_body1 * p1_ws);
			Vec4 p2_ls = Vec4(inverse_transform_body2 * p2_ws);

			// Check if we have a close contact point from last update
			bool lambda_set = false;
			for (const CachedContactPoint *ccp = ccp_start; ccp < ccp_end; ccp++)
				if (Vec4::sLoadFloat4(&ccp->mPosition1).IsClose(p1_ls, mPhysicsSettings.mContactPointPreserveLambdaMaxDistSq)
					&& Vec4::sLoadFloat4(&ccp->mPosition2).IsClose(p2_ls, mPhysicsSettings.mContactPointPreserveLambdaMaxDistSq))
				{
					// Get lambdas from previous frame
					wcp.mNonPenetrationConstraint.SetTotalLambda(ccp->mNonPenetrationLambda);
					wcp.mFrictionConstraint1.SetTotalLambda(ccp->mFrictionLambda[0]);
					wcp.mFrictionConstraint2.SetTotalLambda(ccp->mFrictionLambda[1]);
					wcp.mFrictionConstraint3.SetTotalLambda(ccp->mFrictionLambda[2]);
					lambda_set = true;
					break;
				}
			if (!lambda_set)
			{
				wcp.mNonPenetrationConstraint.SetTotalLambda(0.0f);
				wcp.mFrictionConstraint1.SetTotalLambda(0.0f);
				wcp.mFrictionConstraint2.SetTotalLambda(0.0f);
				wcp.mFrictionConstraint3.SetTotalLambda(0.0f);
			}

			// Create new contact point
			CachedContactPoint &cp = new_manifold->mContactPoints[i];
			p1_ls.StoreFloat4(&cp.mPosition1);
			p2_ls.StoreFloat4(&cp.mPosition2);
			wcp.mContactPoint = &cp;

			// Setup velocity constraint
			wcp.TemplatedCalculateFrictionAndNonPenetrationConstraintProperties<Type1, Type2>(delta_time, gravity, inBody1, inBody2, constraint.mInvMass1, constraint.mInvMass2, inv_i1, inv_i2, p1_ws, p2_ws, inManifold.mWorldSpaceNormal, t1, t2, t3, settings, mPhysicsSettings.mMinVelocityForRestitution);
		}

	#ifdef JPH_DEBUG_RENDERER
		// Draw the manifold
		if (sDrawContactManifolds)
			constraint.Draw(DebugRenderer::sInstance, Color::sOrange);
	#endif // JPH_DEBUG_RENDERER

	#ifdef JPH_TRACK_SIMULATION_STATS
		// Track new contact constraints
		if constexpr (Type1 != EMotionType::Static)
			inBody1.GetMotionPropertiesUnchecked()->GetSimulationStats().mNumContactConstraints.fetch_add(1, memory_order_relaxed);
		if constexpr (Type2 != EMotionType::Static)
			inBody2.GetMotionPropertiesUnchecked()->GetSimulationStats().mNumContactConstraints.fetch_add(1, memory_order_relaxed);
	#endif
	}
	else
	{
		// Store the contact manifold in the cache
		for (int i = 0; i < num_contact_points; ++i)
		{
			// Convert to local space to the body
			Vec4 p1 = Vec4(inverse_transform_body1 * (inManifold.mBaseOffset + inManifold.mRelativeContactPointsOn1[i]));
			Vec4 p2 = Vec4(inverse_transform_body2 * (inManifold.mBaseOffset + inManifold.mRelativeContactPointsOn2[i]));

			// Create new contact point
			CachedContactPoint &cp = new_manifold->mContactPoints[i];
			p1.StoreFloat4(&cp.mPosition1);
			p2.StoreFloat4(&cp.mPosition2);

			// Reset contact impulses, we haven't applied any
			cp.mNonPenetrationLambda = 0.0f;
			cp.mFrictionLambda[0] = 0.0f;
			cp.mFrictionLambda[1] = 0.0f;
			cp.mFrictionLambda[2] = 0.0f;
		}
	}

	// Store cached contact point in body pair cache
	CachedBodyPair *cbp = reinterpret_cast<CachedBodyPair *>(inBodyPairHandle);
	new_manifold->mNextWithSameBodyPair = cbp->mFirstCachedManifold;
	cbp->mFirstCachedManifold = write_cache.ToHandle(new_manifold_kv);

	// A contact constraint was added
	return contact_constraint_created;
}

bool ContactConstraintManager::AddContactConstraint(ContactAllocator &ioContactAllocator, BodyPairHandle inBodyPairHandle, Body &inBody1, Body &inBody2, const ContactManifold &inManifold)
{
	JPH_PROFILE_FUNCTION();

	JPH_DET_LOG("AddContactConstraint: id1: " << inBody1.GetID() << " id2: " << inBody2.GetID()
		<< " subshape1: " << inManifold.mSubShapeID1 << " subshape2: " << inManifold.mSubShapeID2
		<< " normal: " << inManifold.mWorldSpaceNormal << " pendepth: " << inManifold.mPenetrationDepth);

	JPH_ASSERT(inManifold.mWorldSpaceNormal.IsNormalized());

	// Swap bodies so that body 1 id < body 2 id
	const ContactManifold *manifold;
	Body *body1, *body2;
	ContactManifold temp;
	if (inBody2.GetID() < inBody1.GetID())
	{
		body1 = &inBody2;
		body2 = &inBody1;
		temp = inManifold.SwapShapes();
		manifold = &temp;
	}
	else
	{
		body1 = &inBody1;
		body2 = &inBody2;
		manifold = &inManifold;
	}

	// Dispatch to the correct templated form
	// Note: Non-dynamic vs non-dynamic can happen in this case due to one body being a sensor, so we need to have an extended switch case here
	switch (body1->GetMotionType())
	{
	case EMotionType::Dynamic:
		{
			switch (body2->GetMotionType())
			{
			case EMotionType::Dynamic:
				return TemplatedAddContactConstraint<EMotionType::Dynamic, EMotionType::Dynamic>(ioContactAllocator, inBodyPairHandle, *body1, *body2, *manifold);

			case EMotionType::Kinematic:
				return TemplatedAddContactConstraint<EMotionType::Dynamic, EMotionType::Kinematic>(ioContactAllocator, inBodyPairHandle, *body1, *body2, *manifold);

			case EMotionType::Static:
				return TemplatedAddContactConstraint<EMotionType::Dynamic, EMotionType::Static>(ioContactAllocator, inBodyPairHandle, *body1, *body2, *manifold);

			default:
				JPH_ASSERT(false);
				break;
			}
			break;
		}

	case EMotionType::Kinematic:
		switch (body2->GetMotionType())
		{
		case EMotionType::Dynamic:
			return TemplatedAddContactConstraint<EMotionType::Kinematic, EMotionType::Dynamic>(ioContactAllocator, inBodyPairHandle, *body1, *body2, *manifold);

		case EMotionType::Kinematic:
			return TemplatedAddContactConstraint<EMotionType::Kinematic, EMotionType::Kinematic>(ioContactAllocator, inBodyPairHandle, *body1, *body2, *manifold);

		case EMotionType::Static:
			return TemplatedAddContactConstraint<EMotionType::Kinematic, EMotionType::Static>(ioContactAllocator, inBodyPairHandle, *body1, *body2, *manifold);

		default:
			JPH_ASSERT(false);
			break;
		}
		break;

	case EMotionType::Static:
		switch (body2->GetMotionType())
		{
		case EMotionType::Dynamic:
			return TemplatedAddContactConstraint<EMotionType::Static, EMotionType::Dynamic>(ioContactAllocator, inBodyPairHandle, *body1, *body2, *manifold);

		case EMotionType::Kinematic:
			return TemplatedAddContactConstraint<EMotionType::Static, EMotionType::Kinematic>(ioContactAllocator, inBodyPairHandle, *body1, *body2, *manifold);

		case EMotionType::Static: // Static vs static not possible
		default:
			JPH_ASSERT(false);
			break;
		}
		break;

	default:
		JPH_ASSERT(false);
		break;
	}

	return false;
}

void ContactConstraintManager::OnCCDContactAdded(ContactAllocator &ioContactAllocator, const Body &inBody1, const Body &inBody2, const ContactManifold &inManifold, ContactSettings &outSettings)
{
	JPH_ASSERT(inManifold.mWorldSpaceNormal.IsNormalized());

	// Calculate contact settings
	outSettings.mCombinedFriction = mCombineFriction(inBody1, inManifold.mSubShapeID1, inBody2, inManifold.mSubShapeID2);
	outSettings.mCombinedRestitution = mCombineRestitution(inBody1, inManifold.mSubShapeID1, inBody2, inManifold.mSubShapeID2);
	outSettings.mIsSensor = false; // For now, no sensors are supported during CCD

	// The remainder of this function only deals with calling contact callbacks, if there's no contact callback we also don't need to do this work
	if (mContactListener != nullptr)
	{
		// Swap bodies so that body 1 id < body 2 id
		const ContactManifold *manifold;
		const Body *body1, *body2;
		ContactManifold temp;
		if (inBody2.GetID() < inBody1.GetID())
		{
			body1 = &inBody2;
			body2 = &inBody1;
			temp = inManifold.SwapShapes();
			manifold = &temp;
		}
		else
		{
			body1 = &inBody1;
			body2 = &inBody2;
			manifold = &inManifold;
		}

		// Calculate hash
		SubShapeIDPair key { body1->GetID(), manifold->mSubShapeID1, body2->GetID(), manifold->mSubShapeID2 };
		uint64 key_hash = key.GetHash();

		// Check if we already created this contact this physics update
		ManifoldCache &write_cache = mCache[mCacheWriteIdx];
		MKVAndCreated new_manifold_kv = write_cache.FindOrCreate(ioContactAllocator, key, key_hash, 0);
		if (new_manifold_kv.second)
		{
			// This contact is new for this physics update, check if previous update we already had this contact.
			const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
			const MKeyValue *old_manifold_kv = read_cache.Find(key, key_hash);
			if (old_manifold_kv == nullptr)
			{
				// New contact
				mContactListener->OnContactAdded(*body1, *body2, *manifold, outSettings);
			}
			else
			{
				// Existing contact
				mContactListener->OnContactPersisted(*body1, *body2, *manifold, outSettings);

				// Mark contact as persisted so that we won't fire OnContactRemoved callbacks
				old_manifold_kv->GetValue().mFlags |= (uint16)CachedManifold::EFlags::ContactPersisted;
			}

			// Check if the cache is full
			if (new_manifold_kv.first != nullptr)
			{
				// We don't store any contact points in this manifold as it is not for caching impulses, we only need to know that the contact was created
				CachedManifold &new_manifold = new_manifold_kv.first->GetValue();
				new_manifold.mContactNormal = { 0, 0, 0, 0 };
				new_manifold.mFlags |= (uint16)CachedManifold::EFlags::CCDContact;
			}
		}
		else
		{
			// Already found this contact this physics update.
			// Note that we can trigger OnContactPersisted multiple times per physics update, but otherwise we have no way of obtaining the settings
			mContactListener->OnContactPersisted(*body1, *body2, *manifold, outSettings);
		}

		// If we swapped body1 and body2 we need to swap the mass scales back
		if (manifold == &temp)
		{
			std::swap(outSettings.mInvMassScale1, outSettings.mInvMassScale2);
			std::swap(outSettings.mInvInertiaScale1, outSettings.mInvInertiaScale2);
			// Note we do not need to negate the relative surface velocity as it is not applied by the CCD collision constraint
		}
	}

	JPH_ASSERT(outSettings.mIsSensor || !(inBody1.IsSensor() || inBody2.IsSensor()), "Sensors cannot be converted into regular bodies by a contact callback!");
}

void ContactConstraintManager::SortContacts(uint32 *inConstraintIdxBegin, uint32 *inConstraintIdxEnd) const
{
	JPH_PROFILE_FUNCTION();

	QuickSort(inConstraintIdxBegin, inConstraintIdxEnd, [this](uint32 inLHS, uint32 inRHS) {
		const ContactConstraint &lhs = mConstraints[inLHS];
		const ContactConstraint &rhs = mConstraints[inRHS];

		// Most of the time the sort key will be different so we sort on that
		if (lhs.mSortKey != rhs.mSortKey)
			return lhs.mSortKey < rhs.mSortKey;

		// If they're equal we use the IDs of body 1 to order
		if (lhs.mBody1 != rhs.mBody1)
			return lhs.mBody1->GetID() < rhs.mBody1->GetID();

		// If they're still equal we use the IDs of body 2 to order
		if (lhs.mBody2 != rhs.mBody2)
			return lhs.mBody2->GetID() < rhs.mBody2->GetID();

		JPH_ASSERT(inLHS == inRHS, "Hash collision, ordering will be inconsistent");
		return false;
	});
}

void ContactConstraintManager::FinalizeContactCacheAndCallContactPointRemovedCallbacks(uint inExpectedNumBodyPairs, uint inExpectedNumManifolds)
{
	JPH_PROFILE_FUNCTION();

#ifdef JPH_ENABLE_ASSERTS
	// Mark cache as finalized
	ManifoldCache &old_write_cache = mCache[mCacheWriteIdx];
	old_write_cache.Finalize();

	// Check that the count of body pairs and manifolds that we tracked outside of the cache (to avoid contention on an atomic) is correct
	JPH_ASSERT(old_write_cache.GetNumBodyPairs() == inExpectedNumBodyPairs);
	JPH_ASSERT(old_write_cache.GetNumManifolds() == inExpectedNumManifolds);
#endif

	// Buffers are now complete, make write buffer the read buffer
	mCacheWriteIdx ^= 1;

	// Get the old read cache / new write cache
	ManifoldCache &old_read_cache = mCache[mCacheWriteIdx];

	// Call the contact point removal callbacks
	if (mContactListener != nullptr)
		old_read_cache.ContactPointRemovedCallbacks(mContactListener);

	// We're done with the old read cache now
	old_read_cache.Clear();

	// Use the amount of contacts from the last iteration to determine the amount of buckets to use in the hash map for the next iteration
	old_read_cache.Prepare(inExpectedNumBodyPairs, inExpectedNumManifolds);
}

bool ContactConstraintManager::WereBodiesInContact(const BodyID &inBody1ID, const BodyID &inBody2ID) const
{
	// The body pair needs to be in the cache and it needs to have a manifold (otherwise it's just a record indicating that there are no collisions)
	const ManifoldCache &read_cache = mCache[mCacheWriteIdx ^ 1];
	BodyPair key;
	if (inBody1ID < inBody2ID)
		key = BodyPair(inBody1ID, inBody2ID);
	else
		key = BodyPair(inBody2ID, inBody1ID);
	uint64 key_hash = key.GetHash();
	const BPKeyValue *kv = read_cache.Find(key, key_hash);
	return kv != nullptr && kv->GetValue().mFirstCachedManifold != ManifoldMap::cInvalidHandle;
}

template <EMotionType Type1, EMotionType Type2>
JPH_INLINE void ContactConstraintManager::sWarmStartConstraint(ContactConstraint &ioConstraint, MotionProperties *ioMotionProperties1, MotionProperties *ioMotionProperties2, float inWarmStartImpulseRatio)
{
	// Calculate tangents (3 of them in 4D)
	Vec4 t1, t2, t3;
	ioConstraint.GetTangents(t1, t2, t3);

	Vec4 ws_normal = ioConstraint.GetWorldSpaceNormal();

	for (WorldContactPoint &wcp : ioConstraint.mContactPoints)
	{
		// Warm starting: Apply impulse from last frame
		if (wcp.mFrictionConstraint1.IsActive() || wcp.mFrictionConstraint2.IsActive() || wcp.mFrictionConstraint3.IsActive())
		{
			wcp.mFrictionConstraint1.TemplatedWarmStart<Type1, Type2>(ioMotionProperties1, ioConstraint.mInvMass1, ioMotionProperties2, ioConstraint.mInvMass2, t1, inWarmStartImpulseRatio);
			wcp.mFrictionConstraint2.TemplatedWarmStart<Type1, Type2>(ioMotionProperties1, ioConstraint.mInvMass1, ioMotionProperties2, ioConstraint.mInvMass2, t2, inWarmStartImpulseRatio);
			wcp.mFrictionConstraint3.TemplatedWarmStart<Type1, Type2>(ioMotionProperties1, ioConstraint.mInvMass1, ioMotionProperties2, ioConstraint.mInvMass2, t3, inWarmStartImpulseRatio);
		}
		wcp.mNonPenetrationConstraint.TemplatedWarmStart<Type1, Type2>(ioMotionProperties1, ioConstraint.mInvMass1, ioMotionProperties2, ioConstraint.mInvMass2, ws_normal, inWarmStartImpulseRatio);
	}
}

template <class MotionPropertiesCallback>
void ContactConstraintManager::WarmStartVelocityConstraints(const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inWarmStartImpulseRatio, MotionPropertiesCallback &ioCallback)
{
	JPH_PROFILE_FUNCTION();

	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		const uint32 *next_constraint = constraint_idx + 1;
		if (next_constraint < inConstraintIdxEnd)
			PrefetchL1(&mConstraints[*next_constraint]);

		ContactConstraint &constraint = mConstraints[*constraint_idx];

		// Fetch bodies
		Body &body1 = *constraint.mBody1;
		EMotionType motion_type1 = body1.GetMotionType();
		MotionProperties *motion_properties1 = body1.GetMotionPropertiesUnchecked();

		Body &body2 = *constraint.mBody2;
		EMotionType motion_type2 = body2.GetMotionType();
		MotionProperties *motion_properties2 = body2.GetMotionPropertiesUnchecked();

		// Dispatch to the correct templated form
		// Note: Warm starting doesn't differentiate between kinematic/static bodies so we handle both as static bodies
		if (motion_type1 == EMotionType::Dynamic)
		{
			if (motion_type2 == EMotionType::Dynamic)
			{
				sWarmStartConstraint<EMotionType::Dynamic, EMotionType::Dynamic>(constraint, motion_properties1, motion_properties2, inWarmStartImpulseRatio);

				ioCallback(motion_properties2);
			}
			else
				sWarmStartConstraint<EMotionType::Dynamic, EMotionType::Static>(constraint, motion_properties1, motion_properties2, inWarmStartImpulseRatio);

			ioCallback(motion_properties1);
		}
		else
		{
			JPH_ASSERT(motion_type2 == EMotionType::Dynamic);

			sWarmStartConstraint<EMotionType::Static, EMotionType::Dynamic>(constraint, motion_properties1, motion_properties2, inWarmStartImpulseRatio);

			ioCallback(motion_properties2);
		}
	}
}

// Specialize for the two body callback types
template void ContactConstraintManager::WarmStartVelocityConstraints<CalculateSolverSteps>(const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inWarmStartImpulseRatio, CalculateSolverSteps &ioCallback);
template void ContactConstraintManager::WarmStartVelocityConstraints<DummyCalculateSolverSteps>(const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd, float inWarmStartImpulseRatio, DummyCalculateSolverSteps &ioCallback);

template <EMotionType Type1, EMotionType Type2>
JPH_INLINE bool ContactConstraintManager::sSolveVelocityConstraint(ContactConstraint &ioConstraint, MotionProperties *ioMotionProperties1, MotionProperties *ioMotionProperties2)
{
	// A manifold owns both bodies while it is solved. Keep their velocities local
	// across its rows, avoiding repeated aliased loads/stores and DOF-mask creation.
	struct LocalVelocity
	{
		explicit LocalVelocity(const MotionProperties *inProperties)
		{
			if (inProperties != nullptr)
			{
				mLinear = inProperties->GetLinearVelocity();
				mAngular = inProperties->GetAngularVelocity();
				mMask = inProperties->GetLinearDOFsMask().ReinterpretAsFloat();
			}
		}
		Vec4 GetPointVelocityCOM(Vec4Arg inR) const { return mLinear + mAngular.Contract(inR); }
		void AddLinearVelocityStep(Vec4Arg inDelta) { mLinear = Vec4::sAnd(mLinear + inDelta, mMask); }
		void SubLinearVelocityStep(Vec4Arg inDelta) { mLinear = Vec4::sAnd(mLinear - inDelta, mMask); }
		void AddAngularVelocityStep(BivecArg inDelta) { mAngular += inDelta; }
		void SubAngularVelocityStep(BivecArg inDelta) { mAngular -= inDelta; }
		Vec4 mLinear = Vec4::sZero();
		Bivec mAngular = Bivec::sZero();
		Vec4 mMask = Vec4::sZero();
	};
	LocalVelocity state1(ioMotionProperties1), state2(ioMotionProperties2);
	bool any_impulse_applied = false;
	auto relative_point_velocity = [&](const WorldContactPoint &point) {
		Vec4 velocity = Vec4::sZero();
		if constexpr (Type1 != EMotionType::Static)
			velocity = state1.GetPointVelocityCOM(Vec4::sLoadFloat4(&point.mR1));
		if constexpr (Type2 != EMotionType::Static)
			velocity -= state2.GetPointVelocityCOM(Vec4::sLoadFloat4(&point.mR2));
		return velocity;
	};

	// Calculate tangents (3 of them in 4D)
	Vec4 t1, t2, t3;
	ioConstraint.GetTangents(t1, t2, t3);

	// First apply all friction constraints (non-penetration is more important than friction)
	// Evaluate the three tangent projections together. Transpose once per manifold,
	// then retain Dot's pairwise sum order for every point.
	Mat44 tangent_projection = Mat44(t1, t2, t3, Vec4::sZero()).Transposed();
	for (WorldContactPoint &wcp : ioConstraint.mContactPoints)
	{
		// Check if friction is enabled
		if (wcp.mFrictionConstraint1.IsActive() || wcp.mFrictionConstraint2.IsActive() || wcp.mFrictionConstraint3.IsActive())
		{
			// With no normal impulse the friction cone has zero radius. If there is
			// no accumulated friction to undo, this point cannot change either body.
			float max_lambda_f = ioConstraint.mCombinedFriction * wcp.mNonPenetrationConstraint.GetTotalLambda();
			if (max_lambda_f == 0.0f
				&& wcp.mFrictionConstraint1.GetTotalLambda() == 0.0f
				&& wcp.mFrictionConstraint2.GetTotalLambda() == 0.0f
				&& wcp.mFrictionConstraint3.GetTotalLambda() == 0.0f)
				continue;

			// Calculate impulse to stop motion in tangential direction
			Vec4 velocity = relative_point_velocity(wcp);
			Lane4 projected = (tangent_projection.GetColumn4(0) * velocity.SplatX() + tangent_projection.GetColumn4(1) * velocity.SplatY())
				+ (tangent_projection.GetColumn4(2) * velocity.SplatZ() + tangent_projection.GetColumn4(3) * velocity.SplatW());
			float lambda1 = wcp.mFrictionConstraint1.GetTotalLambdaFromVelocity(projected.GetX());
			float lambda2 = wcp.mFrictionConstraint2.GetTotalLambdaFromVelocity(projected.GetY());
			float lambda3 = wcp.mFrictionConstraint3.GetTotalLambdaFromVelocity(projected.GetZ());

			// TODO(4D): the Coulomb friction cone constrains the *magnitude* of the tangential impulse,
			// which now has 3 components instead of 2. The cone is still |lambda_t| <= mu * lambda_n; only
			// the dimension of the disc it is clamped to changes.
			float total_lambda_sq = Square(lambda1) + Square(lambda2) + Square(lambda3);

			// Calculate max impulse that can be applied. Note that we're using the non-penetration impulse from the previous iteration here.
			// We do this because non-penetration is more important so is solved last (the last things that are solved in an iterative solver
			// contribute the most).
			// If the total lambda that we will apply is too large, scale it back
			if (total_lambda_sq > Square(max_lambda_f))
			{
				float scale = max_lambda_f / sqrt(total_lambda_sq);
				lambda1 *= scale;
				lambda2 *= scale;
				lambda3 *= scale;
			}

			// Apply the friction impulse
			if (wcp.mFrictionConstraint1.TemplatedSolveVelocityConstraintApplyLambda<Type1, Type2>(&state1, ioConstraint.mInvMass1, &state2, ioConstraint.mInvMass2, t1, lambda1))
				any_impulse_applied = true;
			if (wcp.mFrictionConstraint2.TemplatedSolveVelocityConstraintApplyLambda<Type1, Type2>(&state1, ioConstraint.mInvMass1, &state2, ioConstraint.mInvMass2, t2, lambda2))
				any_impulse_applied = true;
			if (wcp.mFrictionConstraint3.TemplatedSolveVelocityConstraintApplyLambda<Type1, Type2>(&state1, ioConstraint.mInvMass1, &state2, ioConstraint.mInvMass2, t3, lambda3))
				any_impulse_applied = true;
		}
	}

	Vec4 ws_normal = ioConstraint.GetWorldSpaceNormal();

	// Then apply all non-penetration constraints
	for (WorldContactPoint &wcp : ioConstraint.mContactPoints)
	{
		// Solve non penetration velocities
		float lambda = Clamp(wcp.mNonPenetrationConstraint.GetTotalLambdaFromVelocity(ws_normal.Dot(relative_point_velocity(wcp))), 0.0f, FLT_MAX);
		if (wcp.mNonPenetrationConstraint.TemplatedSolveVelocityConstraintApplyLambda<Type1, Type2>(&state1, ioConstraint.mInvMass1, &state2, ioConstraint.mInvMass2, ws_normal, lambda))
			any_impulse_applied = true;
	}

	if (any_impulse_applied)
	{
		if constexpr (Type1 == EMotionType::Dynamic)
			ioMotionProperties1->SetVelocityStep(state1.mLinear, state1.mAngular);
		if constexpr (Type2 == EMotionType::Dynamic)
			ioMotionProperties2->SetVelocityStep(state2.mLinear, state2.mAngular);
	}
	return any_impulse_applied;
}

bool ContactConstraintManager::SolveVelocityConstraints(const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd)
{
	JPH_PROFILE_FUNCTION();

	bool any_impulse_applied = false;

	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		const uint32 *next_constraint = constraint_idx + 1;
		if (next_constraint < inConstraintIdxEnd)
			PrefetchL1(&mConstraints[*next_constraint]);

		ContactConstraint &constraint = mConstraints[*constraint_idx];

		// Fetch bodies
		Body &body1 = *constraint.mBody1;
		EMotionType motion_type1 = body1.GetMotionType();
		MotionProperties *motion_properties1 = body1.GetMotionPropertiesUnchecked();

		Body &body2 = *constraint.mBody2;
		EMotionType motion_type2 = body2.GetMotionType();
		MotionProperties *motion_properties2 = body2.GetMotionPropertiesUnchecked();

		// Dispatch to the correct templated form
		switch (motion_type1)
		{
		case EMotionType::Dynamic:
			switch (motion_type2)
			{
			case EMotionType::Dynamic:
				any_impulse_applied |= sSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Dynamic>(constraint, motion_properties1, motion_properties2);
				break;

			case EMotionType::Kinematic:
				any_impulse_applied |= sSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Kinematic>(constraint, motion_properties1, motion_properties2);
				break;

			case EMotionType::Static:
				any_impulse_applied |= sSolveVelocityConstraint<EMotionType::Dynamic, EMotionType::Static>(constraint, motion_properties1, motion_properties2);
				break;

			default:
				JPH_ASSERT(false);
				break;
			}
			break;

		case EMotionType::Kinematic:
			JPH_ASSERT(motion_type2 == EMotionType::Dynamic);
			any_impulse_applied |= sSolveVelocityConstraint<EMotionType::Kinematic, EMotionType::Dynamic>(constraint, motion_properties1, motion_properties2);
			break;

		case EMotionType::Static:
			JPH_ASSERT(motion_type2 == EMotionType::Dynamic);
			any_impulse_applied |= sSolveVelocityConstraint<EMotionType::Static, EMotionType::Dynamic>(constraint, motion_properties1, motion_properties2);
			break;

		default:
			JPH_ASSERT(false);
			break;
		}
	}

	return any_impulse_applied;
}

void ContactConstraintManager::StoreAppliedImpulses(const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd) const
{
	// Copy back total applied impulse to cache for the next frame
	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		const uint32 *next_constraint = constraint_idx + 1;
		if (next_constraint < inConstraintIdxEnd)
			PrefetchL1(&mConstraints[*next_constraint]);

		const ContactConstraint &constraint = mConstraints[*constraint_idx];

		for (const WorldContactPoint &wcp : constraint.mContactPoints)
		{
			wcp.mContactPoint->mNonPenetrationLambda = wcp.mNonPenetrationConstraint.GetTotalLambda();
			wcp.mContactPoint->mFrictionLambda[0] = wcp.mFrictionConstraint1.GetTotalLambda();
			wcp.mContactPoint->mFrictionLambda[1] = wcp.mFrictionConstraint2.GetTotalLambda();
			wcp.mContactPoint->mFrictionLambda[2] = wcp.mFrictionConstraint3.GetTotalLambda();
		}
	}
}

bool ContactConstraintManager::SolvePositionConstraints(const uint32 *inConstraintIdxBegin, const uint32 *inConstraintIdxEnd)
{
	JPH_PROFILE_FUNCTION();

	bool any_impulse_applied = false;

	for (const uint32 *constraint_idx = inConstraintIdxBegin; constraint_idx < inConstraintIdxEnd; ++constraint_idx)
	{
		const uint32 *next_constraint = constraint_idx + 1;
		if (next_constraint < inConstraintIdxEnd)
			PrefetchL1(&mConstraints[*next_constraint]);

		ContactConstraint &constraint = mConstraints[*constraint_idx];

		// Fetch bodies
		Body &body1 = *constraint.mBody1;
		Body &body2 = *constraint.mBody2;

		// Get transforms
		RMat44 transform1 = body1.GetCenterOfMassTransform();
		RMat44 transform2 = body2.GetCenterOfMassTransform();

		Vec4 ws_normal = constraint.GetWorldSpaceNormal();

		for (WorldContactPoint &wcp : constraint.mContactPoints)
		{
			// Calculate new contact point positions in world space (the bodies may have moved)
			RVec4 p1 = transform1 * Vec4::sLoadFloat4(&wcp.mContactPoint->mPosition1);
			RVec4 p2 = transform2 * Vec4::sLoadFloat4(&wcp.mContactPoint->mPosition2);

			// Calculate separation along the normal (negative if interpenetrating)
			// Allow a little penetration by default (PhysicsSettings::mPenetrationSlop) to avoid jittering between contact/no-contact which wipes out the contact cache and warm start impulses
			// Clamp penetration to a max PhysicsSettings::mMaxPenetrationDistance so that we don't apply a huge impulse if we're penetrating a lot
			float separation = max(Vec4(p2 - p1).Dot(ws_normal) + mPhysicsSettings.mPenetrationSlop, -mPhysicsSettings.mMaxPenetrationDistance);

			// Only enforce constraint when separation < 0 (otherwise we're apart)
			if (separation < 0.0f)
			{
				// Update constraint properties (bodies may have moved)
				wcp.CalculateNonPenetrationConstraintProperties(body1, constraint.mInvMass1, constraint.mInvInertiaScale1, body2, constraint.mInvMass2, constraint.mInvInertiaScale2, p1, p2, ws_normal);

				// Solve position errors
				if (wcp.mNonPenetrationConstraint.SolvePositionConstraintWithMassOverride(body1, constraint.mInvMass1, body2, constraint.mInvMass2, ws_normal, separation, mPhysicsSettings.mBaumgarte))
					any_impulse_applied = true;
			}
		}
	}

	return any_impulse_applied;
}

void ContactConstraintManager::RecycleConstraintBuffer()
{
	// Reset constraint array
	mNumConstraints = 0;
}

void ContactConstraintManager::FinishConstraintBuffer()
{
	// Free constraints buffer
	mUpdateContext->mTempAllocator->Free(mConstraints, mMaxConstraints * sizeof(ContactConstraint));
	mConstraints = nullptr;
	mNumConstraints = 0;

	// Reset update context
	mUpdateContext = nullptr;
}

void ContactConstraintManager::SaveState(StateRecorder &inStream, const StateRecorderFilter *inFilter) const
{
	mCache[mCacheWriteIdx ^ 1].SaveState(inStream, inFilter);
}

bool ContactConstraintManager::RestoreState(StateRecorder &inStream, const StateRecorderFilter *inFilter)
{
	bool success = mCache[mCacheWriteIdx].RestoreState(mCache[mCacheWriteIdx ^ 1], inStream, inFilter);

	// If this is the last part, the cache is finalized
	if (inStream.IsLastPart())
	{
		mCacheWriteIdx ^= 1;
		mCache[mCacheWriteIdx].Clear();
	}

	return success;
}

JPH_NAMESPACE_END
