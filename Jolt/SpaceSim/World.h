// SPDX-License-Identifier: MIT
#pragma once

#include <Jolt/Jolt.h>
#include <Jolt/Math/Vec4.h>
#include <Jolt/Math/Rotor.h>
#include <Jolt/Math/Bivec.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace JPH::SpaceSim
{
using ObjectID = uint64_t;
using FrameID = uint64_t;
using Tick = uint64_t;
using Vector = std::array<double, 4>;

// This is a residency/trajectory foundation, not a second rigid-body integrator.
// All interfaces are boundary operations; callers must quiesce the physics owner.
// No global coordinate is flattened to a floating-point scalar.
// Contract/limits:
// - IEEE binary64, round-to-nearest, no fast-math/reassociation or denormal flush.
// - Positions are sector integers plus bounded offsets; local snapshots and
//   frame-origin differences are limited to +/-65536 on EACH of four axes.
// - A translation is analytic force-free drift, not a gravity/orbit model.
//   Attitude/angular state remains stamped with GetSourceEpoch(): reconstruct
//   nontrivial attitude in the physics owner before detailed solving. Promotion
//   alone does not claim to integrate an asymmetric body's rotational dynamics.
// - Reject overflow, nonfinite values, intervals >2^53 ticks, and unrepresentable
//   solve frames. No clamping, freezing, dropped pairs or duplicate body owners.
// - Representation revisions invalidate caller-owned collision/history caches.
//   This class neither owns Jolt bodies/caches nor promises allocation recovery
//   or a production object-capacity/performance envelope.
static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53, "SpaceSim requires IEEE binary64");
inline constexpr double cMaxLocalRange = 65536.0;

inline bool IsFinite(const Vector &v)
{
	for (double x : v)
		if (!std::isfinite(x))
			return false;
	return true;
}

struct GlobalPosition
{
	static constexpr double cSectorSize = 1048576.0; // Power of two, canonical offset [0, size).
	std::array<int64_t, 4> mSector {};
	Vector mOffset {};

	bool IsValid() const
	{
		for (double x : mOffset)
			if (!std::isfinite(x) || x < 0 || x >= cSectorSize)
				return false;
		return true;
	}

	static bool AddSector(int64_t a, int64_t b, int64_t &out)
	{
		if ((b > 0 && a > INT64_MAX - b) || (b < 0 && a < INT64_MIN - b))
			return false;
		out = a + b;
		return true;
	}

	// Checked and transactional, including aliasing out == *this.
	bool Translate(const Vector &delta, GlobalPosition &out) const
	{
		if (!IsValid() || !IsFinite(delta))
			return false;
		GlobalPosition result = *this;
		for (size_t i = 0; i < 4; ++i)
		{
			double whole = std::trunc(delta[i] / cSectorSize);
			if (whole < -0x1p63 || whole >= 0x1p63)
				return false;
			int64_t sector;
			if (!AddSector(mSector[i], static_cast<int64_t>(whole), sector))
				return false;
			double offset = mOffset[i] + std::fmod(delta[i], cSectorSize);
			if (offset < 0)
			{
				if (!AddSector(sector, -1, sector)) return false;
				offset += cSectorSize;
			}
			if (offset >= cSectorSize)
			{
				if (!AddSector(sector, 1, sector)) return false;
				offset -= cSectorSize;
			}
			result.mSector[i] = sector;
			result.mOffset[i] = offset;
		}
		out = result;
		return true;
	}

	// Only bounded differences enter local doubles. The unsigned subtraction
	// handles INT64_MIN vs INT64_MAX without signed overflow or lossy conversion.
	bool Difference(const GlobalPosition &reference, double range, Vector &out) const
	{
		if (!IsValid() || !reference.IsValid() || !std::isfinite(range) || range < 0 || range > cMaxLocalRange)
			return false;
		Vector result;
		for (size_t i = 0; i < 4; ++i)
		{
			bool positive = mSector[i] >= reference.mSector[i];
			uint64_t magnitude = positive? uint64_t(mSector[i]) - uint64_t(reference.mSector[i])
				: uint64_t(reference.mSector[i]) - uint64_t(mSector[i]);
			if (magnitude > 1) return false;
			double d = magnitude == 0? mOffset[i] - reference.mOffset[i]
				: positive? (cSectorSize - reference.mOffset[i]) + mOffset[i]
				: -((cSectorSize - mOffset[i]) + reference.mOffset[i]);
			if (std::abs(d) > range) return false;
			result[i] = d;
		}
		out = result;
		return true;
	}
};

// Individual double operations are bounded, not arbitrary precision. Transfers
// reject an addition that discards more than this absolute length/speed budget.
inline constexpr double cTransferError = 1.0e-9;
inline bool SplitSum(double a, double b, double &head, double &tail)
{
	head = a + b;
	if (!std::isfinite(head)) return false;
	double z = head - a;
	tail = (a - (head - z)) + (b - z);
	return std::isfinite(tail);
}
inline bool CheckedSum(double a, double b, double &out)
{
	double head, tail;
	if (!SplitSum(a, b, head, tail) || std::abs(tail) > cTransferError) return false;
	out = head;
	return true;
}

// Tick differences, not absolute ticks, are converted. Very old/new snapshots
// are rejected if their interval cannot be represented exactly as an integer.
inline bool Elapsed(Tick from, Tick to, double seconds_per_tick, double &out)
{
	uint64_t count = to >= from? to - from : from - to;
	if (!std::isfinite(seconds_per_tick) || seconds_per_tick <= 0 || count > (UINT64_C(1) << 53)) return false;
	out = double(count) * seconds_per_tick;
	if (to < from) out = -out;
	return std::isfinite(out);
}

inline double ElapsedRemainder(Tick from, Tick to, double step, double elapsed)
{
	uint64_t count = to >= from? to - from : from - to;
	double tail = std::fma(double(count), step, -std::abs(elapsed));
	return to < from? -tail : tail;
}

// Retain the product residual in the bounded offset instead of discarding it
// in a huge double displacement (e.g. (1e18 + 128) * 1.5 loses 64 otherwise).
inline bool TranslateProduct(const GlobalPosition &position, const Vector &velocity, double time, GlobalPosition &out)
{
	Vector head, tail;
	for (size_t i = 0; i < 4; ++i)
	{
		head[i] = velocity[i] * time;
		if (!std::isfinite(head[i])) return false;
		tail[i] = std::fma(velocity[i], time, -head[i]);
	}
	GlobalPosition next;
	if (!position.Translate(head, next) || !next.Translate(tail, next)) return false;
	out = next;
	return true;
}

struct Frame
{
	FrameID mID = 0;
	GlobalPosition mOrigin;
	Vector mVelocity {};
	Tick mEpoch = 0;

	bool IsValid() const { return mID != 0 && mOrigin.IsValid() && IsFinite(mVelocity); }

	bool OriginAt(Tick tick, double step, GlobalPosition &out) const
	{
		double dt;
		if (!IsValid() || !Elapsed(mEpoch, tick, step, dt)) return false;
		GlobalPosition next;
		if (!TranslateProduct(mOrigin, mVelocity, dt, next)
			|| !TranslateProduct(next, mVelocity, ElapsedRemainder(mEpoch, tick, step, dt), next)) return false;
		out = next;
		return true;
	}

	// A global-sector position in a co-moving index whose reference origin stays
	// at its own epoch. Subtract reference VELOCITY before multiplying by time.
	// Equal enormous velocities therefore cancel without swallowing local drift.
	bool IndexOrigin(const Frame &reference, Tick tick, double step, GlobalPosition &out) const
	{
		if (!reference.IsValid()) return false;
		GlobalPosition anchor;
		double dt;
		if (!OriginAt(reference.mEpoch, step, anchor) || !Elapsed(reference.mEpoch, tick, step, dt)) return false;
		Vector relative, relative_tail;
		for (size_t i = 0; i < 4; ++i)
			if (!SplitSum(mVelocity[i], -reference.mVelocity[i], relative[i], relative_tail[i])) return false;
		double time_tail = ElapsedRemainder(reference.mEpoch, tick, step, dt);
		GlobalPosition next;
		if (!TranslateProduct(anchor, relative, dt, next)
			|| !TranslateProduct(next, relative_tail, dt, next)
			|| !TranslateProduct(next, relative, time_tail, next)
			|| !TranslateProduct(next, relative_tail, time_tail, next)) return false;
		out = next;
		return true;
	}
};

enum class Tier : uint8_t { Detailed, Trajectory, Distant };
enum DemotionBlocker : uint32_t
{
	ActiveContact = 1,
	PendingImpulse = 2,
	WeaponInteraction = 4,
	TopologyTransaction = 8
};

struct ObjectState
{
	ObjectID mID = 0;
	FrameID mFrame = 0;
	Tick mEpoch = 0; // Translation epoch.
	// Lazy source stamp survives analytic resampling/frame transfers. Attitude
	// and the supplied acceleration envelope still refer to this source epoch.
	Tick mSourceEpoch = 0;
	bool mHasSourceEpoch = false;
	Tick GetSourceEpoch() const { return mHasSourceEpoch? mSourceEpoch : mEpoch; }
	Vector mPosition {};
	Vector mVelocity {};
	Rotor mOrientation = Rotor::sIdentity();
	Bivec mAngularVelocity = Bivec::sZero();
	double mMass = 1;
	double mHealth = 1;
	uint64_t mTopologyVersion = 0;
	uint64_t mHistoryRevision = 0;
	uint64_t mPopulationSeed = 0;
	Tier mTier = Tier::Trajectory;
	Tick mLastPromotion = 0;
	uint32_t mDemotionBlockers = 0;
	double mRotationRadius = 0; // Full hull radius, covers every attitude during interval.
	double mAccelerationBound = 0; // Supplied bound, not a gravitational/force law.
	double mWeaponReach = 0;
	Tick mSchedulingInterval = 0;

	bool IsValid() const
	{
		if (mID == 0 || mFrame == 0 || !IsFinite(mPosition) || !IsFinite(mVelocity)
			|| !std::isfinite(mMass) || mMass <= 0 || !std::isfinite(mHealth)
			|| mHealth < 0 || !mOrientation.IsNormalized()
			|| !std::isfinite(mRotationRadius) || mRotationRadius < 0
			|| !std::isfinite(mAccelerationBound) || mAccelerationBound < 0
			|| !std::isfinite(mWeaponReach) || mWeaponReach < 0
			|| (mTier != Tier::Detailed && mTier != Tier::Trajectory && mTier != Tier::Distant)) return false;
		for (size_t i = 0; i < 4; ++i)
			if (std::abs(mPosition[i]) > cMaxLocalRange) return false;
		for (uint i = 0; i < 6; ++i)
			if (!std::isfinite(mAngularVelocity[i])) return false;
		return true;
	}

	// Analytic force-free TRANSLATION, one correctly rounded fma from the stored
	// epoch. Attitude/angular state is retained as an epoch snapshot, NOT falsely
	// integrated as constant omega for an asymmetric torque-free rigid body.
	bool Evaluate(Tick tick, double step, ObjectState &out) const
	{
		double dt;
		if (!IsValid() || !Elapsed(mEpoch, tick, step, dt)) return false;
		ObjectState result = *this;
		result.mSourceEpoch = GetSourceEpoch();
		result.mHasSourceEpoch = true;
		double time_tail = ElapsedRemainder(mEpoch, tick, step, dt);
		for (size_t i = 0; i < 4; ++i)
		{
			result.mPosition[i] = std::fma(mVelocity[i], dt, mPosition[i]);
			if (!CheckedSum(result.mPosition[i], mVelocity[i] * time_tail, result.mPosition[i])) return false;
			if (!std::isfinite(result.mPosition[i]) || std::abs(result.mPosition[i]) > cMaxLocalRange) return false;
		}
		result.mEpoch = tick;
		out = result;
		return true;
	}
};

// Jolt stores local velocities and (in the default build) positions as floats.
// Conversion is explicit and budgeted; DOUBLE_PRECISION positions do not make
// velocities double. No global or huge common velocity may enter this boundary.
struct PhysicsLimits
{
	double mMaxPosition = 8192;
	double mMaxVelocity = 8192;
	double mPositionError = 1.0e-4;
	double mVelocityError = 1.0e-6;
};
struct PhysicsState
{
	ObjectID mID = 0;
	FrameID mFrame = 0;
	Tick mTick = 0;
	Vec4 mPosition = Vec4::sZero();
	Vec4 mVelocity = Vec4::sZero();
};
inline bool TryPhysicsState(const ObjectState &object, const PhysicsLimits &limits, PhysicsState &out)
{
	if (!object.IsValid() || !std::isfinite(limits.mMaxPosition) || limits.mMaxPosition <= 0 || limits.mMaxPosition > cMaxLocalRange
		|| !std::isfinite(limits.mMaxVelocity) || limits.mMaxVelocity <= 0 || limits.mMaxVelocity > std::numeric_limits<float>::max()
		|| !std::isfinite(limits.mPositionError) || limits.mPositionError < 0
		|| !std::isfinite(limits.mVelocityError) || limits.mVelocityError < 0) return false;
	std::array<float, 4> position, velocity;
	for (size_t i = 0; i < 4; ++i)
	{
		if (std::abs(object.mPosition[i]) > limits.mMaxPosition || std::abs(object.mVelocity[i]) > limits.mMaxVelocity) return false;
		position[i] = static_cast<float>(object.mPosition[i]);
		velocity[i] = static_cast<float>(object.mVelocity[i]);
		if (std::abs(double(position[i]) - object.mPosition[i]) > limits.mPositionError
			|| std::abs(double(velocity[i]) - object.mVelocity[i]) > limits.mVelocityError) return false;
	}
	out = {object.mID, object.mFrame, object.mEpoch,
		Vec4(position[0], position[1], position[2], position[3]), Vec4(velocity[0], velocity[1], velocity[2], velocity[3])};
	return true;
}

// Bounds are in one co-moving index coordinate system, never absolute doubles.
struct GlobalBounds
{
	GlobalPosition mMin, mMax;
	static bool Less(const GlobalPosition &a, const GlobalPosition &b, size_t axis)
	{
		return a.mSector[axis] != b.mSector[axis]? a.mSector[axis] < b.mSector[axis] : a.mOffset[axis] < b.mOffset[axis];
	}
	bool Overlaps(const GlobalBounds &other) const
	{
		for (size_t i = 0; i < 4; ++i)
			if (Less(mMax, other.mMin, i) || Less(other.mMax, mMin, i)) return false;
		return true;
	}
	void Encapsulate(const GlobalBounds &other)
	{
		for (size_t i = 0; i < 4; ++i)
		{
			if (Less(other.mMin, mMin, i)) { mMin.mSector[i] = other.mMin.mSector[i]; mMin.mOffset[i] = other.mMin.mOffset[i]; }
			if (Less(mMax, other.mMax, i)) { mMax.mSector[i] = other.mMax.mSector[i]; mMax.mOffset[i] = other.mMax.mOffset[i]; }
		}
	}
	bool Contains(const GlobalPosition &point) const { return Overlaps({point, point}); }
};

namespace Detail
{
struct Interval { double mMin, mMax; };
inline double Down(double x) { return std::nextafter(x, -std::numeric_limits<double>::infinity()); }
inline double Up(double x) { return std::nextafter(x, std::numeric_limits<double>::infinity()); }
inline Interval Product(Interval a, Interval b)
{
	if ((a.mMin == 0 && a.mMax == 0) || (b.mMin == 0 && b.mMax == 0)) return {0, 0};
	std::array<double, 4> p {a.mMin * b.mMin, a.mMin * b.mMax, a.mMax * b.mMin, a.mMax * b.mMax};
	return {Down(*std::min_element(p.begin(), p.end())), Up(*std::max_element(p.begin(), p.end()))};
}
inline Interval Difference(double a, double b)
{
	double s = a - b, z = s - a;
	double error = (a - (s - z)) + (-b - z);
	return {error < 0? Down(s) : s, error > 0? Up(s) : s};
}
inline bool Time(Tick from, Tick to, double step, Interval &out)
{
	double dt;
	if (!Elapsed(from, to, step, dt)) return false;
	out = from == to? Interval{0, 0} : Interval{Down(dt), Up(dt)};
	return true;
}
inline bool Shift(GlobalBounds &box, const std::array<Interval, 4> &delta)
{
	Vector lo, hi;
	for (size_t i = 0; i < 4; ++i) { lo[i] = delta[i].mMin; hi[i] = delta[i].mMax; }
	// Canonical-offset addition has at most a few ulps at sector scale. Pad
	// BEFORE the next operation, including sector carry, not just nextafter(0).
	// 2^-28 > 8 ulps at size 2^20; this avoids asymmetric rounding at boundaries.
	constexpr double pad = 0x1p-28;
	GlobalPosition low, high;
	if (!box.mMin.Translate(lo, low) || !box.mMax.Translate(hi, high)
		|| !low.Translate({-pad, -pad, -pad, -pad}, low)
		|| !high.Translate({pad, pad, pad, pad}, high)) return false;
	box = {low, high};
	return true;
}
inline bool Endpoint(const ObjectState &object, const Frame &frame, const Frame &reference, Tick tick, double step, GlobalBounds &out)
{
	Interval anchor_time, relative_time, local_time;
	if (!Time(frame.mEpoch, reference.mEpoch, step, anchor_time)
		|| !Time(reference.mEpoch, tick, step, relative_time)
		|| !Time(object.mEpoch, tick, step, local_time)) return false;
	GlobalBounds result {frame.mOrigin, frame.mOrigin};
	std::array<Interval, 4> delta;
	for (size_t i = 0; i < 4; ++i) delta[i] = Product({frame.mVelocity[i], frame.mVelocity[i]}, anchor_time);
	if (!Shift(result, delta)) return false;
	for (size_t i = 0; i < 4; ++i) delta[i] = Product(Difference(frame.mVelocity[i], reference.mVelocity[i]), relative_time);
	if (!Shift(result, delta)) return false;
	for (size_t i = 0; i < 4; ++i) delta[i] = {object.mPosition[i], object.mPosition[i]};
	if (!Shift(result, delta)) return false;
	for (size_t i = 0; i < 4; ++i) delta[i] = Product({object.mVelocity[i], object.mVelocity[i]}, local_time);
	if (!Shift(result, delta)) return false;
	out = result;
	return true;
}
}

// Sweeps the nominal force-free path. a is a supplied bound on unmodeled
// acceleration since the STATE epoch: pad 1/2 a T^2 over the whole age/horizon.
// Hull rotation is covered by its full radius, not start/end orientation alone.
// This is a conservative candidate bound, NOT a time-of-impact or hit test.
inline bool BuildSweptBounds(const ObjectState &object, const Frame &frame, const Frame &reference,
	Tick tick, Tick interval, double step, GlobalBounds &out)
{
	Tick horizon = std::max(interval, object.mSchedulingInterval);
	if (!object.IsValid() || !frame.IsValid() || !reference.IsValid() || object.mFrame != frame.mID || tick > UINT64_MAX - horizon) return false;
	Tick end = tick + horizon;
	Detail::Interval start_age, end_age;
	if (!Detail::Time(object.GetSourceEpoch(), tick, step, start_age) || !Detail::Time(object.GetSourceEpoch(), end, step, end_age)) return false;
	double age = std::max({std::abs(start_age.mMin), std::abs(start_age.mMax), std::abs(end_age.mMin), std::abs(end_age.mMax)});
	double acceleration = 0;
	if (object.mAccelerationBound != 0)
	{
		// Multiply by time before halving: 0.5 * a can underflow even when
		// a * T * T is enormous. Reject overflow rather than masking it with 0.
		auto velocity_bound = Detail::Product({object.mAccelerationBound, object.mAccelerationBound}, {age, age});
		auto displacement = Detail::Product(velocity_bound, {age, age});
		acceleration = Detail::Product(displacement, {0.5, 0.5}).mMax;
	}
	double radius = Detail::Up(Detail::Up(object.mRotationRadius + object.mWeaponReach) + acceleration);
	if (!std::isfinite(radius)) return false;
	std::array<Detail::Interval, 4> padding;
	padding.fill({-radius, radius});
	GlobalBounds first, last;
	if (!Detail::Endpoint(object, frame, reference, tick, step, first)
		|| !Detail::Endpoint(object, frame, reference, end, step, last)) return false;
	first.Encapsulate(last);
	if (!Detail::Shift(first, padding)) return false;
	// Validate the physical global range too: a co-moving index must not hide
	// overflow of an actual global position near the int64 sector boundary.
	Frame stationary;
	stationary.mID = 1;
	stationary.mEpoch = reference.mEpoch;
	GlobalBounds absolute_first, absolute_last;
	if (!Detail::Endpoint(object, frame, stationary, tick, step, absolute_first)
		|| !Detail::Endpoint(object, frame, stationary, end, step, absolute_last)) return false;
	absolute_first.Encapsulate(absolute_last);
	if (!Detail::Shift(absolute_first, padding)) return false;
	out = first;
	return true;
}

using InteractionPair = std::pair<ObjectID, ObjectID>;
struct SolveGroup
{
	FrameID mFrame = 0;
	std::vector<ObjectID> mObjects; // Sorted; every connected component has one solve frame.
};
struct ActivationResult
{
	Tick mTick = 0;
	Tick mInterval = 0;
	uint64_t mResidentCount = 0;
	uint64_t mNodeVisits = 0;
	std::vector<InteractionPair> mPairs;
	std::vector<ObjectID> mPromoted;
	std::vector<SolveGroup> mGroups;
};

// Balanced 4D AABB hierarchy, rebuilt per query: O(N log^2 N) construction,
// spatially pruned traversal, O(N) nodes. Dense overlap can necessarily emit
// O(N^2) pairs; there is no unconditional all-pairs loop. Split the widest
// sector/offset axis so a population separated only in W is still pruned.
namespace Detail
{
class SpatialIndex
{
public:
	struct Item { ObjectID mID; GlobalBounds mBounds; };
	explicit SpatialIndex(std::vector<Item> items) : mItems(std::move(items))
	{
		if (!mItems.empty()) Build(0, mItems.size());
	}
	void Find(ActivationResult &out) const
	{
		if (mNodes.empty()) return;
		for (const Item &item : mItems) Query(item, 0, out);
		std::sort(out.mPairs.begin(), out.mPairs.end());
	}
private:
	struct Node
	{
		GlobalBounds mBounds;
		ObjectID mMaxID = 0;
		size_t mLeft = SIZE_MAX, mRight = SIZE_MAX, mItem = SIZE_MAX;
	};
	size_t Build(size_t begin, size_t end)
	{
		Node node;
		node.mBounds = mItems[begin].mBounds;
		for (size_t i = begin; i < end; ++i)
		{
			node.mBounds.Encapsulate(mItems[i].mBounds);
			node.mMaxID = std::max(node.mMaxID, mItems[i].mID);
		}
		size_t index = mNodes.size();
		mNodes.push_back(node);
		if (end - begin == 1) { mNodes[index].mItem = begin; return index; }
		size_t axis = 0;
		uint64_t best_sectors = 0;
		double best_offset = -1;
		for (size_t i = 0; i < 4; ++i)
		{
			uint64_t sectors = uint64_t(node.mBounds.mMax.mSector[i]) - uint64_t(node.mBounds.mMin.mSector[i]);
			double offset = node.mBounds.mMax.mOffset[i] - node.mBounds.mMin.mOffset[i];
			if (sectors > best_sectors || (sectors == best_sectors && offset > best_offset))
			{ axis = i; best_sectors = sectors; best_offset = offset; }
		}
		std::sort(mItems.begin() + begin, mItems.begin() + end, [axis](const Item &a, const Item &b)
		{
			if (GlobalBounds::Less(a.mBounds.mMin, b.mBounds.mMin, axis)) return true;
			if (GlobalBounds::Less(b.mBounds.mMin, a.mBounds.mMin, axis)) return false;
			return a.mID < b.mID;
		});
		size_t middle = begin + (end - begin) / 2;
		size_t left = Build(begin, middle), right = Build(middle, end);
		mNodes[index].mLeft = left;
		mNodes[index].mRight = right;
		return index;
	}
	void Query(const Item &item, size_t index, ActivationResult &out) const
	{
		++out.mNodeVisits;
		const Node &node = mNodes[index];
		if (node.mMaxID <= item.mID || !node.mBounds.Overlaps(item.mBounds)) return;
		if (node.mItem != SIZE_MAX) { out.mPairs.emplace_back(item.mID, mItems[node.mItem].mID); return; }
		Query(item, node.mLeft, out);
		Query(item, node.mRight, out);
	}
	std::vector<Item> mItems;
	std::vector<Node> mNodes;
};
}

class World
{
public:
	explicit World(double seconds_per_tick = 1.0 / 120.0) : mSecondsPerTick(seconds_per_tick) { }
	bool IsValid() const { return std::isfinite(mSecondsPerTick) && mSecondsPerTick > 0; }
	bool AddFrame(const Frame &frame)
	{
		return IsValid() && frame.IsValid() && mFrames.emplace(frame.mID, frame).second;
	}
	bool AddObject(const ObjectState &object)
	{
		if (!IsValid() || !object.IsValid()) return false;
		const Frame *frame = FindFrame(object.mFrame);
		GlobalPosition physical;
		return frame != nullptr && frame->OriginAt(object.mEpoch, mSecondsPerTick, physical)
			&& physical.Translate(object.mPosition, physical) && mObjects.emplace(object.mID, object).second;
	}
	// Only the detailed physics owner publishes an integrated pose. A committed
	// snapshot establishes a NEW attitude/envelope epoch; an analytic resample
	// does not. Tier or frame changes must use their atomic boundary operations.
	bool CommitDetailedState(const ObjectState &state)
	{
		const ObjectState *old = FindObject(state.mID);
		if (!IsValid() || !state.IsValid() || old == nullptr || old->mTier != Tier::Detailed || state.mTier != Tier::Detailed
			|| state.mFrame != old->mFrame || state.mEpoch < old->mEpoch || mRepresentationVersion == UINT64_MAX) return false;
		GlobalPosition physical;
		if (!FindFrame(state.mFrame)->OriginAt(state.mEpoch, mSecondsPerTick, physical) || !physical.Translate(state.mPosition, physical)) return false;
		ObjectState next = state;
		next.mLastPromotion = old->mLastPromotion;
		next.mSourceEpoch = state.mEpoch;
		next.mHasSourceEpoch = false;
		mObjects.find(state.mID)->second = next;
		++mRepresentationVersion;
		return true;
	}

	const ObjectState *FindObject(ObjectID id) const
	{
		auto it = mObjects.find(id);
		return it == mObjects.end()? nullptr : &it->second;
	}
	const Frame *FindFrame(FrameID id) const
	{
		auto it = mFrames.find(id);
		return it == mFrames.end()? nullptr : &it->second;
	}
	uint64_t GetRepresentationVersion() const { return mRepresentationVersion; }

	// Long-range force-free translation is evaluated directly from the stored
	// epoch into sectors. Unlike a local physics snapshot it may cross arbitrarily
	// many sectors (within the checked int64/displacement/time envelope).
	bool EvaluateGlobal(ObjectID id, Tick tick, GlobalPosition &out) const
	{
		const ObjectState *object = FindObject(id);
		if (!IsValid() || object == nullptr) return false;
		const Frame *frame = FindFrame(object->mFrame);
		GlobalPosition position;
		double dt;
		if (!frame->OriginAt(tick, mSecondsPerTick, position)
			|| !position.Translate(object->mPosition, position)
			|| !Elapsed(object->mEpoch, tick, mSecondsPerTick, dt)
			|| !TranslateProduct(position, object->mVelocity, dt, position)
			|| !TranslateProduct(position, object->mVelocity, ElapsedRemainder(object->mEpoch, tick, mSecondsPerTick, dt), position)) return false;
		out = position;
		return true;
	}

	// Caller supplies the COMPLETE interacting set, including projectiles and
	// pending contacts, at a quiescent boundary. Failure changes nothing. Success
	// bumps a revision: physics/query/history owners must rebuild their caches
	// before continuing. This class does not own or silently mutate Jolt caches.
	bool TransformInteractingSet(const std::vector<ObjectID> &ids, FrameID destination, Tick tick)
	{
		auto target = FindFrame(destination);
		if (!IsValid() || target == nullptr || ids.empty() || mRepresentationVersion == UINT64_MAX) return false;
		std::set<ObjectID> unique;
		std::vector<ObjectState> staged;
		for (ObjectID id : ids)
		{
			auto object = FindObject(id);
			if (object == nullptr || !unique.insert(id).second) return false;
			ObjectState next;
			if (!object->Evaluate(tick, mSecondsPerTick, next)) return false;
			auto source = FindFrame(object->mFrame);
			GlobalPosition indexed, absolute;
			Vector delta;
			if (!source->OriginAt(tick, mSecondsPerTick, absolute)
				|| !target->OriginAt(tick, mSecondsPerTick, absolute)
				|| !source->IndexOrigin(*target, tick, mSecondsPerTick, indexed)
				|| !indexed.Difference(target->mOrigin, cMaxLocalRange, delta)) return false;
			for (size_t i = 0; i < 4; ++i)
			{
				double relative, relative_tail, velocity, velocity_tail, tail;
				if (!CheckedSum(next.mPosition[i], delta[i], next.mPosition[i])
					|| !SplitSum(source->mVelocity[i], -target->mVelocity[i], relative, relative_tail)
					|| !SplitSum(relative, next.mVelocity[i], velocity, velocity_tail)
					|| !CheckedSum(relative_tail, velocity_tail, tail)
					|| !CheckedSum(velocity, tail, next.mVelocity[i])) return false;
			}
			next.mFrame = destination;
			if (!next.IsValid() || !absolute.Translate(next.mPosition, absolute)) return false;
			staged.push_back(next);
		}
		for (const ObjectState &object : staged) mObjects.find(object.mID)->second = object;
		++mRepresentationVersion;
		return true;
	}

	// Read-only plan, sorted by stable IDs. The index frame is deterministic and
	// removes common velocity before the sweep; it is unrelated to any camera.
	bool FindInteractions(Tick tick, Tick interval, ActivationResult &out) const
	{
		if (!IsValid() || tick > UINT64_MAX - interval) return false;
		ActivationResult result;
		result.mTick = tick;
		result.mInterval = interval;
		result.mResidentCount = uint64_t(mObjects.size());
		if (mObjects.empty()) { out = std::move(result); return true; }
		const Frame &reference = mFrames.begin()->second;
		std::vector<Detail::SpatialIndex::Item> items;
		items.reserve(mObjects.size());
		for (const auto &entry : mObjects)
		{
			GlobalBounds bounds;
			if (!BuildSweptBounds(entry.second, *FindFrame(entry.second.mFrame), reference, tick, interval, mSecondsPerTick, bounds)) return false;
			items.push_back({entry.first, bounds});
		}
		Detail::SpatialIndex(std::move(items)).Find(result);
		std::map<ObjectID, ObjectID> parents;
		auto root = [&parents](ObjectID id)
		{
			ObjectID r = id;
			while (parents[r] != r) r = parents[r];
			while (parents[id] != id) { ObjectID next = parents[id]; parents[id] = r; id = next; }
			return r;
		};
		for (const auto &pair : result.mPairs)
		{
			parents.emplace(pair.first, pair.first);
			parents.emplace(pair.second, pair.second);
			ObjectID a = root(pair.first), b = root(pair.second);
			parents[std::max(a, b)] = std::min(a, b);
		}
		std::map<ObjectID, std::vector<ObjectID>> components;
		for (const auto &entry : parents) components[root(entry.first)].push_back(entry.first);
		for (auto &entry : components) result.mGroups.push_back({FindObject(entry.first)->mFrame, std::move(entry.second)});
		out = std::move(result);
		return true;
	}

	// Commit promotion and coherent frame ownership BEFORE the caller performs
	// collision/weapon work. A connected component cannot be independently solved
	// by two regions. Any unrepresentable component aborts the entire operation.
	bool PrepareInteractions(Tick tick, Tick interval, ActivationResult &out)
	{
		ActivationResult result;
		if (!FindInteractions(tick, interval, result)) return false;
		World staged = *this;
		for (const auto &group : result.mGroups)
		{
			if (!staged.TransformInteractingSet(group.mObjects, group.mFrame, tick)) return false;
			for (ObjectID id : group.mObjects)
			{
				ObjectState &object = staged.mObjects.find(id)->second;
				if (object.mTier != Tier::Detailed)
				{
					object.mTier = Tier::Detailed;
					object.mLastPromotion = tick;
					result.mPromoted.push_back(id);
				}
			}
		}
		std::sort(result.mPromoted.begin(), result.mPromoted.end());
		*this = std::move(staged);
		out = std::move(result);
		return true;
	}

	// Hysteresis is an integer residence duration. All blockers and a fresh
	// conservative future query must be clear; merely being off camera is never
	// sufficient. A failed demotion leaves epoch, metadata and tier untouched.
	bool TryDemote(ObjectID id, Tier target, Tick tick, Tick minimum_residence, Tick lookahead)
	{
		const ObjectState *object = FindObject(id);
		if (object == nullptr || (target != Tier::Trajectory && target != Tier::Distant)
			|| uint8_t(target) <= uint8_t(object->mTier) || object->mDemotionBlockers != 0
			|| tick < object->mEpoch || tick < object->mLastPromotion || tick - object->mLastPromotion < minimum_residence
			|| mRepresentationVersion == UINT64_MAX) return false;
		ActivationResult plan;
		if (!FindInteractions(tick, lookahead, plan)) return false;
		for (const auto &pair : plan.mPairs)
			if (pair.first == id || pair.second == id) return false;
		ObjectState next;
		if (!object->Evaluate(tick, mSecondsPerTick, next)) return false;
		next.mTier = target;
		mObjects.find(id)->second = next;
		++mRepresentationVersion;
		return true;
	}

private:
	double mSecondsPerTick;
	uint64_t mRepresentationVersion = 0;
	std::map<FrameID, Frame> mFrames;
	std::map<ObjectID, ObjectState> mObjects;
};
}
