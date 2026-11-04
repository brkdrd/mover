// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT
#pragma once
#include <Jolt/Jolt.h>
#include <Jolt/Math/Bivec.h>
#include <Jolt/Math/Rotor.h>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>

JPH_NAMESPACE_BEGIN
namespace SpaceSim {

/// Combat query acceptance envelope (metres, seconds; caller supplies coherent units):
/// - Authoritative union of sharp oriented 4-boxes; no convex rounding/render meshes.
/// - Constant world acceleration and angular velocity; duration [0, 1] seconds.
/// - Position, speed, acceleration norms <= 1e6; angular Frobenius speed <= 10000/s
///   and angular travel <= 16 radians. These bounds are admission limits, not a
///   promise that every admitted input converges within the configured work budget.
/// - Half extents >= .005, component COM radius |localCenter|+|halfExtent| <= 256;
///   sphere radius [1e-4, 256]. Physical feature width is at least .01.
/// - Distance tolerance [1e-6, .05]; numerical allowance (0, .05], default .001 each.
///   A Hit may be early by this spatial allowance; no uniform time error exists at
///   tangency. Narrow phase uses double arithmetic around float Rotor evaluation;
///   64 float eps * (1 + relative distance + component radius) guards its error.
///   Float world witness quantization is included in the reported allowance. Rebase
///   if PrecisionExceeded: raising tolerance is an explicit accuracy-policy change.
/// - <= 1024 bodies, 256 boxes/body, 4096 total boxes, 1e6 distance evaluations.
///   Default evaluation budget 4096 is GLOBAL, including simultaneous components.
/// - Exact finite epoch/frame agreement; unique body IDs and unique component IDs
///   within a body (zero is permitted). Preserve IDs across proxy rebuild/reordering.
///
/// This is a standalone, read-only projectile/proxy query, NOT general rotating
/// convex/compound-pair CCD, a PhysicsSystem integrator, a spatial index, or a
/// production workload-capacity guarantee. Core Jolt still owns finite-mass ramming
/// and coupled contact response. Geometry/state arrays must remain alive and immutable
/// throughout the query; stale frame/topology snapshots must be invalidated by owner.

/// These are authoritative sharp 4-boxes, NOT visibility or end-pose proxies.
/// Accepted trajectories: x(t)=x0+v0*t+0.5*a*t^2, R(t)=exp(Omega*t/2)*R0.
/// a and WORLD Omega are constant over the entire query. Gyroscopic/steered motion
/// must be resampled by the physics owner into coherent intervals; it is not inferred.
struct CombatMotion
{
	uint64 mFrameID = 0;
	double mEpoch = 0;
	Vec4 mPosition = Vec4::sZero();
	Vec4 mVelocity = Vec4::sZero();
	Vec4 mAcceleration = Vec4::sZero();
	Rotor mRotation = Rotor::sIdentity();
	Bivec mAngularVelocity = Bivec::sZero();
};
struct CombatBox
{
	uint64 mComponentID = 0;
	Vec4 mLocalCenter = Vec4::sZero();
	Rotor mLocalRotation = Rotor::sIdentity();
	Vec4 mHalfExtent = Vec4::sReplicate(0.5f);
};
struct CombatCompound
{
	uint64 mBodyID = 0;
	CombatMotion mMotion;
	const CombatBox *mBoxes = nullptr;
	uint32 mBoxCount = 0;
};
struct CombatProjectile
{
	CombatMotion mMotion;
	double mRadius = 0.1;
};
enum class CombatQueryStatus { NoHit, Hit, Indeterminate, InvalidInput };
enum class CombatQueryFailure { None, InvalidData, FrameEpochMismatch, OutsideEnvelope, BudgetExceeded, NumericalStall, PrecisionExceeded };
struct CombatQuery
{
	uint64 mFrameID = 0;
	double mEpoch = 0;
	double mDuration = 1.0 / 120.0;
	double mDistanceTolerance = 0.001;
	double mMaxNumericalAllowance = 0.001;
	uint32 mMaxDistanceEvaluations = 4096;
};
struct CombatQueryTelemetry
{
	uint32 mDistanceEvaluations = 0;
	uint32 mBroadPhaseRejected = 0;
	uint32 mCandidateBoxes = 0;
	uint32 mBudgetFailures = 0;
};
struct CombatHit
{
	CombatQueryStatus mStatus = CombatQueryStatus::Indeterminate;
	CombatQueryFailure mFailure = CombatQueryFailure::None;
	uint64 mBodyID = 0;
	uint64 mComponentID = 0;
	double mTime = 0; ///< Hit: seconds after epoch. Failure: diagnostic stop time ONLY, not a globally cleared interval.
	double mSeparation = 0; ///< Point-to-box set distance minus radius; not signed depth for an embedded center.
	double mNumericalAllowance = 0; ///< Distance evaluation + world-space float witness allowance.
	Vec4 mPosition = Vec4::sZero(); ///< Point on target surface.
	Vec4 mNormal = Vec4::sAxisX(); ///< Outward, target toward projectile.
	CombatQueryTelemetry mTelemetry;
};

namespace CombatDetail {
using D4 = std::array<double, 4>;
inline D4 D(Vec4Arg v) { return { v[0], v[1], v[2], v[3] }; }
inline D4 Add(const D4 &a, const D4 &b) { return { a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3] }; }
inline D4 Sub(const D4 &a, const D4 &b) { return { a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3] }; }
inline D4 Mul(const D4 &a, double s) { return { a[0]*s, a[1]*s, a[2]*s, a[3]*s }; }
inline double Dot(const D4 &a, const D4 &b) { double s = 0; for (uint i=0; i<4; ++i) s += a[i]*b[i]; return s; }
inline double Length(const D4 &a) { return std::sqrt(Dot(a,a)); }
inline Vec4 V(const D4 &a) { return Vec4(float(a[0]),float(a[1]),float(a[2]),float(a[3])); }
inline D4 Position(const CombatMotion &m, double t) { return Add(Add(D(m.mPosition), Mul(D(m.mVelocity),t)), Mul(D(m.mAcceleration),0.5*t*t)); }
inline Rotor Orientation(const CombatMotion &m, double t) { return Rotor::sExp(m.mAngularVelocity * float(0.5*t)) * m.mRotation; }
inline double AngularSpeed(const CombatMotion &m)
{
	double sum = 0;
	for (uint i=0; i<6; ++i) sum += double(m.mAngularVelocity[i])*m.mAngularVelocity[i];
	return std::sqrt(2*sum); // Frobenius norm of the skew matrix bounds its operator norm.
}
inline D4 Rotate(RotorArg r, const D4 &v)
{
	D4 out = {0,0,0,0};
	for (uint i=0; i<4; ++i) { Vec4 axis = Vec4::sZero(); axis[i] = 1; out = Add(out, Mul(D(r * axis), v[i])); }
	return out;
}
struct Distance
{
	double mGap;
	double mRoundoff;
	D4 mPoint;
	D4 mNormal;
};
// Exact Euclidean point-to-sharp-box distance: clamp each of FOUR local coordinates.
// No GJK, no inflated-box corner approximation.
inline Distance PointBox(const CombatProjectile &p, const CombatCompound &b, const CombatBox &box, double t)
{
	Rotor body = Orientation(b.mMotion,t).Normalized();
	Rotor rotation = (body * box.mLocalRotation).Normalized();
	D4 center = Add(Position(b.mMotion,t), Rotate(body,D(box.mLocalCenter)));
	D4 delta = Sub(Position(p.mMotion,t),center), local, closest;
	for (uint i=0; i<4; ++i)
	{
		Vec4 axis = Vec4::sZero(); axis[i] = 1;
		local[i] = Dot(delta,D(rotation * axis));
		closest[i] = std::max(-double(box.mHalfExtent[i]), std::min(double(box.mHalfExtent[i]),local[i]));
	}
	D4 diff = Sub(local,closest);
	double distance = Length(diff);
	D4 normal;
	if (distance > 0) normal = Mul(diff,1.0/distance);
	else
	{
		// Initially embedded: deterministic nearest face, including exact center.
		uint face = 0;
		for (uint i=1; i<4; ++i)
			if (box.mHalfExtent[i]-std::abs(local[i]) < box.mHalfExtent[face]-std::abs(local[face])) face = i;
		normal = {0,0,0,0}; normal[face] = local[face] < 0? -1.0 : 1.0;
		closest[face] = normal[face]*box.mHalfExtent[face];
	}
	return { distance-p.mRadius, 64.0*std::numeric_limits<float>::epsilon()*(1+Length(delta)+Length(D(box.mHalfExtent))+Length(D(box.mLocalCenter))),
		Add(center,Rotate(rotation,closest)), Rotate(rotation,normal) };
}
struct Bounds { D4 mMin, mMax; };
inline Bounds SweptBounds(const CombatMotion &m, double radius, double duration)
{
	D4 start = Position(m,0), end = Position(m,duration);
	Bounds bounds;
	for (uint i=0; i<4; ++i)
	{
		// A parabola differs from its endpoint chord by at most |a_i| h^2/8.
		double pad = radius + std::abs(double(m.mAcceleration[i]))*duration*duration/8;
		bounds.mMin[i] = std::min(start[i],end[i])-pad;
		bounds.mMax[i] = std::max(start[i],end[i])+pad;
	}
	return bounds;
}
inline bool Candidate(const CombatProjectile &p, const CombatCompound &target, const CombatBox &box, const CombatQuery &q)
{
	// A COM-centred sphere bounds EVERY orientation of this authoritative component,
	// including orbital motion of localCenter. No endpoint OBB shortcut is used.
	double radius = Length(D(box.mLocalCenter))+Length(D(box.mHalfExtent));
	double separation_bound = Length(Sub(D(p.mMotion.mPosition),D(target.mMotion.mPosition)))
		+ q.mDuration*Length(Sub(D(p.mMotion.mVelocity),D(target.mMotion.mVelocity)))
		+ 0.5*q.mDuration*q.mDuration*Length(Sub(D(p.mMotion.mAcceleration),D(target.mMotion.mAcceleration)));
	double roundoff = 64.0*std::numeric_limits<float>::epsilon()*(1+separation_bound+2*radius);
	Bounds a = SweptBounds(p.mMotion,p.mRadius+q.mDistanceTolerance+roundoff,q.mDuration);
	Bounds b = SweptBounds(target.mMotion,radius,q.mDuration);
	for (uint i=0; i<4; ++i) if (a.mMax[i] < b.mMin[i] || b.mMax[i] < a.mMin[i]) return false;
	return true;
}
inline bool Finite(Vec4Arg v)
{
	for (uint i=0; i<4; ++i) if (!std::isfinite(v[i])) return false;
	return true;
}
inline bool ValidRotor(RotorArg r)
{
	for (uint i=0; i<8; ++i) if (!std::isfinite(r.mValue[i])) return false;
	return r.IsNormalized(1.0e-5f);
}
inline CombatQueryFailure ValidateMotion(const CombatMotion &m, const CombatQuery &q)
{
	if (!std::isfinite(m.mEpoch) || !Finite(m.mPosition) || !Finite(m.mVelocity) || !Finite(m.mAcceleration) || !ValidRotor(m.mRotation))
		return CombatQueryFailure::InvalidData;
	for (uint i=0; i<6; ++i) if (!std::isfinite(m.mAngularVelocity[i])) return CombatQueryFailure::InvalidData;
	if (m.mAngularVelocity.GetVec8()[6] != 0 || m.mAngularVelocity.GetVec8()[7] != 0) return CombatQueryFailure::InvalidData;
	if (m.mFrameID != q.mFrameID || m.mEpoch != q.mEpoch) return CombatQueryFailure::FrameEpochMismatch;
	if (Length(D(m.mPosition)) > 1.0e6 || Length(D(m.mVelocity)) > 1.0e6 || Length(D(m.mAcceleration)) > 1.0e6
		|| AngularSpeed(m) > 10000 || AngularSpeed(m)*q.mDuration > 16)
		return CombatQueryFailure::OutsideEnvelope;
	return CombatQueryFailure::None;
}
inline CombatQueryFailure Validate(const CombatProjectile &p, const CombatCompound *targets, uint32 count, const CombatQuery &q)
{
	if (!std::isfinite(q.mEpoch) || !std::isfinite(q.mDuration) || q.mDuration < 0 || !std::isfinite(q.mDistanceTolerance) || q.mDistanceTolerance <= 0
		|| !std::isfinite(q.mMaxNumericalAllowance) || q.mMaxNumericalAllowance <= 0
		|| !std::isfinite(p.mRadius) || p.mRadius <= 0 || (count != 0 && targets == nullptr)) return CombatQueryFailure::InvalidData;
	if (q.mDuration > 1 || q.mDistanceTolerance < 1.0e-6 || q.mDistanceTolerance > 0.05 || p.mRadius < 1.0e-4 || p.mRadius > 256
		|| q.mMaxNumericalAllowance > 0.05 || count > 1024 || q.mMaxDistanceEvaluations > 1000000) return CombatQueryFailure::OutsideEnvelope;
	CombatQueryFailure failure = ValidateMotion(p.mMotion,q);
	if (failure != CombatQueryFailure::None) return failure;
	uint32 total = 0;
	for (uint32 b=0; b<count; ++b)
	{
		const CombatCompound &target = targets[b];
		if (target.mBoxCount != 0 && target.mBoxes == nullptr) return CombatQueryFailure::InvalidData;
		if (target.mBoxCount > 256 || target.mBoxCount > 4096-total) return CombatQueryFailure::OutsideEnvelope;
		total += target.mBoxCount;
		failure = ValidateMotion(target.mMotion,q);
		if (failure != CombatQueryFailure::None) return failure;
		for (uint32 prev=0; prev<b; ++prev) if (targets[prev].mBodyID == target.mBodyID) return CombatQueryFailure::InvalidData;
		for (uint32 c=0; c<target.mBoxCount; ++c)
		{
			const CombatBox &box = target.mBoxes[c];
			if (!Finite(box.mHalfExtent) || !Finite(box.mLocalCenter) || !ValidRotor(box.mLocalRotation)) return CombatQueryFailure::InvalidData;
			for (uint i=0; i<4; ++i) if (box.mHalfExtent[i] <= 0) return CombatQueryFailure::InvalidData;
			for (uint i=0; i<4; ++i) if (box.mHalfExtent[i] < 0.005f) return CombatQueryFailure::OutsideEnvelope;
			if (Length(D(box.mHalfExtent))+Length(D(box.mLocalCenter)) > 256) return CombatQueryFailure::OutsideEnvelope;
			for (uint32 prev=0; prev<c; ++prev) if (target.mBoxes[prev].mComponentID == box.mComponentID) return CombatQueryFailure::InvalidData;
		}
	}
	return CombatQueryFailure::None;
}
inline CombatHit Failure(CombatHit result, CombatQueryFailure why, double time = 0)
{
	result.mStatus = why == CombatQueryFailure::InvalidData || why == CombatQueryFailure::FrameEpochMismatch? CombatQueryStatus::InvalidInput : CombatQueryStatus::Indeterminate;
	result.mFailure = why;
	result.mTime = time;
	result.mBodyID = 0; result.mComponentID = 0;
	result.mPosition = Vec4::sZero(); result.mNormal = Vec4::sZero(); result.mSeparation = 0; result.mNumericalAllowance = 0;
	if (why == CombatQueryFailure::BudgetExceeded) ++result.mTelemetry.mBudgetFailures;
	return result;
}
} // namespace CombatDetail

/// Bounded, allocation-free query. Only Hit is consumable as an impact. A Hit is a
/// conservative spatial-tolerance event (gap <= tolerance + numerical allowance),
/// NOT an exact analytic TOI. Grazing can have large time uncertainty. Indeterminate
/// requires subdivision/resampling or authoritative fail-closed policy, never NoHit.
/// Exact time ties use (body ID, component ID), independent of array order.
inline CombatHit SweepProjectile(const CombatProjectile &p, const CombatCompound *targets, uint32 count, const CombatQuery &q)
{
	using namespace CombatDetail;
	CombatHit result;
	CombatQueryFailure validation = Validate(p,targets,count,q);
	if (validation != CombatQueryFailure::None) return Failure(result,validation);
	result.mStatus = CombatQueryStatus::NoHit;
	result.mTime = q.mDuration;
	for (uint32 b=0; b<count; ++b)
		for (uint32 c=0; c<targets[b].mBoxCount; ++c)
		{
			const CombatCompound &target = targets[b];
			const CombatBox &box = target.mBoxes[c];
			if (!Candidate(p,target,box,q)) { ++result.mTelemetry.mBroadPhaseRejected; continue; }
			++result.mTelemetry.mCandidateBoxes;
			D4 relative_v = Sub(D(p.mMotion.mVelocity),D(target.mMotion.mVelocity));
			D4 relative_a = Sub(D(p.mMotion.mAcceleration),D(target.mMotion.mAcceleration));
			// Norm of affine relative velocity is convex: its maximum is at an endpoint.
			// Hausdorff speed of the box adds ||Omega_skew|| * max radius about COM.
			double speed = std::max(Length(relative_v),Length(Add(relative_v,Mul(relative_a,q.mDuration))))
				+ AngularSpeed(target.mMotion)*(Length(D(box.mLocalCenter))+Length(D(box.mHalfExtent)));
			double time = 0;
			for (;;)
			{
				if (result.mTelemetry.mDistanceEvaluations >= q.mMaxDistanceEvaluations)
					return Failure(result,CombatQueryFailure::BudgetExceeded,time);
				++result.mTelemetry.mDistanceEvaluations;
				Distance d = PointBox(p,target,box,time);
				if (d.mGap <= q.mDistanceTolerance+d.mRoundoff)
				{
					double allowance = d.mRoundoff+2*std::numeric_limits<float>::epsilon()*(1+Length(d.mPoint));
					if (allowance > q.mMaxNumericalAllowance) return Failure(result,CombatQueryFailure::PrecisionExceeded,time);
					if (result.mStatus != CombatQueryStatus::Hit || time < result.mTime
						|| (time == result.mTime && (target.mBodyID < result.mBodyID || (target.mBodyID == result.mBodyID && box.mComponentID < result.mComponentID))))
					{
						result.mStatus = CombatQueryStatus::Hit;
						result.mTime = time;
						result.mBodyID = target.mBodyID;
						result.mComponentID = box.mComponentID;
						result.mNormal = V(d.mNormal);
						result.mPosition = V(d.mPoint);
						result.mSeparation = d.mGap;
						result.mNumericalAllowance = allowance;
					}
					break;
				}
				if (speed == 0) break;
				double step = 0.9*(d.mGap-d.mRoundoff)/speed;
				// Do not truncate at the current winner: subtracting/re-adding that floating
				// time can discard an identical simultaneous component by one ULP.
				if (step > q.mDuration-time) break;
				if (!(time+step > time)) return Failure(result,CombatQueryFailure::NumericalStall,time);
				time += step;
			}
		}
	return result;
}

enum class CombatResponseOwner { PhysicsSolver, QueryEvent };
struct CombatLaunch
{
	bool mValid = false;
	Vec4 mPosition = Vec4::sZero();
	Vec4 mProjectileVelocity = Vec4::sZero();
	Vec4 mRecoilImpulse = Vec4::sZero();
	Bivec mRecoilAngularImpulse = Bivec::sZero();
};
struct CombatImpactOutcome
{
	bool mValid = false;
	bool mApplyImpulse = false;
	Vec4 mTargetImpulse = Vec4::sZero();
	Bivec mTargetAngularImpulse = Bivec::sZero();
	double mProjectileEnergyLoss = 0;
};
/// Inputs are world-frame COM offsets and an ejection DELTA velocity, not a final
/// muzzle-relative speed after recoil. Assumes carried projectile mass has already
/// been detached with inherited point momentum by the mass/topology owner. It then
/// accelerates that mass by deltaV and supplies equal/opposite impulse J and J wedge r.
/// Does not remove mass, mutate a Body, or claim a variable-mass rocket solve.
inline CombatLaunch PrepareCombatLaunch(const CombatMotion &ship, Vec4Arg offset, Vec4Arg delta_velocity, double mass)
{
	using namespace CombatDetail;
	CombatLaunch launch;
	CombatQuery q; q.mFrameID = ship.mFrameID; q.mEpoch = ship.mEpoch; q.mDuration = 0;
	if (ValidateMotion(ship,q) != CombatQueryFailure::None || !Finite(offset) || !Finite(delta_velocity) || !std::isfinite(mass) || mass <= 0 || mass > 1.0e12
		|| Length(D(offset)) > 256 || Length(D(delta_velocity)) > 1.0e6) return launch;
	launch.mPosition = ship.mPosition+offset;
	launch.mProjectileVelocity = ship.mVelocity+ship.mAngularVelocity.Contract(offset)+delta_velocity;
	launch.mRecoilImpulse = V(Mul(D(delta_velocity),-mass));
	launch.mRecoilAngularImpulse = Bivec::sWedge(launch.mRecoilImpulse,offset);
	launch.mValid = Finite(launch.mProjectileVelocity) && Finite(launch.mRecoilImpulse) && !launch.mRecoilAngularImpulse.IsNaN();
	return launch;
}

/// Pure response plan for already co-framed event-time velocities and COM offset.
/// Solver-owned contacts expose energy accounting but ZERO extra impulse, so damage
/// cannot accidentally reapply the solver response. Query events transfer the modeled
/// incoming-minus-outgoing projectile momentum. EnergyLoss is frame-dependent and is
/// NOT damage energy: target kinetic work/material rules still belong to the owner.
/// The coordinator must deduplicate event IDs and invalidate trajectories after applying
/// this plan once. This helper is not an event ledger, ramming solver or ricochet model.
inline CombatImpactOutcome PrepareCombatImpact(CombatResponseOwner owner, double mass, Vec4Arg incoming, Vec4Arg outgoing, Vec4Arg offset)
{
	using namespace CombatDetail;
	CombatImpactOutcome impact;
	if ((owner != CombatResponseOwner::PhysicsSolver && owner != CombatResponseOwner::QueryEvent) || !std::isfinite(mass) || mass <= 0 || mass > 1.0e12
		|| !Finite(incoming) || !Finite(outgoing) || !Finite(offset) || Length(D(incoming)) > 1.0e6 || Length(D(outgoing)) > 1.0e6 || Length(D(offset)) > 256) return impact;
	impact.mValid = true;
	impact.mProjectileEnergyLoss = 0.5*mass*(Dot(D(incoming),D(incoming))-Dot(D(outgoing),D(outgoing)));
	if (owner == CombatResponseOwner::QueryEvent)
	{
		impact.mApplyImpulse = true;
		impact.mTargetImpulse = V(Mul(Sub(D(incoming),D(outgoing)),mass));
		impact.mTargetAngularImpulse = Bivec::sWedge(impact.mTargetImpulse,offset);
	}
	return impact;
}

} // namespace SpaceSim
JPH_NAMESPACE_END
