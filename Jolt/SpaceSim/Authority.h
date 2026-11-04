// Engine-support contracts, not a transport or a second physics integrator.
#pragma once
#include <Jolt/Jolt.h>
#include <Jolt/Math/Rotor.h>
#include <cstring>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <functional>
#include <memory>
#include <tuple>
#include <vector>

namespace JPH::SpaceSim
{
using Tick = uint64;
using ID = uint64;

namespace Detail
{
// v1 prelude + four count fields, and exact fixed element wire sizes.
inline constexpr size_t cFixedBytes = 192, cEntityBytes = 276, cCommandBytes = 72, cStreamBytes = 32, cEventBytes = 108;
}

struct Limits
{
	size_t entities = 1024;
	size_t commands = 1024;
	size_t streams = 1024;
	size_t events = 1024;
	size_t history = 32;
	size_t allocation_bytes = 8 * 1024 * 1024;
	size_t bytes = 4 * 1024 * 1024;
	Tick future_ticks = 120, replay_ticks = 120;
};

struct Command
{
	ID client = 0, entity = 0;
	uint64 sequence = 0;
	Tick tick = 0;
	std::array<float, 4> translation {};
	std::array<float, 6> rotation {};
};

struct StreamState
{
	ID client = 0, entity = 0;
	uint64 accepted = 0, applied = 0;
};

struct EntityState
{
	ID id = 0, owner = 0, frame = 1, topology_version = 0, geometry = 1;
	std::array<double, 4> position {};
	std::array<float, 8> orientation { 1, 0, 0, 0, 0, 0, 0, 0 };
	std::array<float, 4> velocity {}, center_of_mass {};
	std::array<float, 6> angular_velocity {}, inertia { 1, 1, 1, 1, 1, 1 };
	std::array<float, 8> inertia_orientation { 1, 0, 0, 0, 0, 0, 0, 0 };
	std::array<float, 10> actuators {};
	float mass = 1, fuel = 0, energy = 0, heat = 0, health = 1;
};

// One nonrotating local frame per coordinator. origin/velocity are global,
// evaluated at epoch (an authoritative tick); entity positions/velocities are local.
struct FrameState
{
	ID id = 1;
	Tick epoch = 0;
	std::array<double, 4> origin {}, velocity {};
};

template <class T, size_t N>
inline bool Finite(const std::array<T, N> &inValues)
{
	for (T v : inValues)
		if (!std::isfinite(v)) return false;
	return true;
}

inline bool ValidCommand(const Command &inCommand)
{
	if (inCommand.client == 0 || inCommand.entity == 0 || inCommand.sequence == 0
		|| !Finite(inCommand.translation) || !Finite(inCommand.rotation)) return false;
	for (float v : inCommand.translation) if (std::abs(v) > 1.0f) return false;
	for (float v : inCommand.rotation) if (std::abs(v) > 1.0f) return false;
	return true;
}

enum class ResponseOwner : uint8 { Solver, Query };

// fraction is an ordering key in millionths of a tick, not a request to substep.
// (tick, id) is the impact identity shared by all producers; conflicting owners
// for that identity are rejected before either response is delivered.
// Query response is applied at the tick-end boundary; true TOI integration is deferred.
struct ImpactEvent
{
	Tick tick = 0;
	uint64 fraction = 1000000;
	ID id = 0, source = 0, target = 0, component = 0, frame = 1;
	ResponseOwner owner = ResponseOwner::Solver;
	std::array<float, 4> impulse {};
	std::array<float, 6> angular_impulse {};
	float energy = 0;
};

// A causal application checkpoint for the fields below, NOT a PhysicsSystem
// rollback image: solver/contact caches, shape assets and external system state
// are not serialized. geometry/model_version identify caller-owned configuration.
struct Checkpoint
{
	Tick tick = 0;
	uint64 random_seed = 1, backlog = 0;
	uint64 physics_hz = 120, resource_period = 1, snapshot_period = 6, model_version = 1;
	FrameState frame;
	std::vector<EntityState> entities;
	std::vector<Command> commands;
	std::vector<StreamState> streams;
	std::vector<ImpactEvent> events;
};

inline bool ValidImpact(const ImpactEvent &e)
{
	return e.tick != 0 && e.id != 0 && e.source != 0 && e.target != 0 && e.frame != 0
		&& e.fraction <= 1000000 && (e.owner == ResponseOwner::Solver || e.owner == ResponseOwner::Query)
		&& Finite(e.impulse) && Finite(e.angular_impulse) && std::isfinite(e.energy) && e.energy >= 0;
}

inline bool ImpactLess(const ImpactEvent &a, const ImpactEvent &b)
{
	return std::tie(a.tick, a.fraction, a.source, a.target, a.component, a.id)
		< std::tie(b.tick, b.fraction, b.source, b.target, b.component, b.id);
}

class ImpactBatch
{
public:
	ImpactBatch(const std::vector<ImpactEvent> &inEvents, const Limits &inLimits)
	{
		if (inEvents.size() > inLimits.events) return;
		for (size_t i = 0; i < inEvents.size(); ++i)
		{
			if (!ValidImpact(inEvents[i])) return;
			for (size_t j = 0; j < i; ++j)
				if (inEvents[i].tick == inEvents[j].tick && inEvents[i].id == inEvents[j].id) return;
		}
		mEvents = inEvents;
		std::sort(mEvents.begin(), mEvents.end(), ImpactLess);
		mValid = true;
	}
	ImpactBatch(const ImpactBatch &) = delete;
	ImpactBatch &operator = (const ImpactBatch &) = delete;

	// Callback failure consumes the batch as well: the enclosing coordinator must
	// fail-stop, never retry a partially delivered external impulse or observation.
	template <class ApplyQuery, class Observe>
	bool Process(ApplyQuery &&inApply, Observe &&inObserve)
	{
		if (!mValid || mConsumed) return false;
		mConsumed = true;
		for (const ImpactEvent &e : mEvents)
		{
			if (e.owner == ResponseOwner::Query && !inApply(e)) return false;
			if (!inObserve(e)) return false;
		}
		return true;
	}
private:
	std::vector<ImpactEvent> mEvents;
	bool mValid = false, mConsumed = false;
};

inline bool FitsStorage(const Checkpoint &s, const Limits &l, bool inAddCommand = false, bool inAddStream = false)
{
	size_t budget = l.allocation_bytes;
	auto account = [&](size_t count, size_t size) { if (count > budget / size) return false; budget -= count * size; return true; };
	if (!account(s.entities.size(), sizeof(EntityState)) || !account(s.commands.size(), sizeof(Command))
		|| !account(s.streams.size(), sizeof(StreamState)) || !account(s.events.size(), sizeof(ImpactEvent))
		|| !account(inAddCommand ? 1 : 0, sizeof(Command)) || !account(inAddStream ? 1 : 0, sizeof(StreamState))) return false;
	budget = l.bytes;
	return account(1, Detail::cFixedBytes) && account(s.entities.size(), Detail::cEntityBytes)
		&& account(s.commands.size(), Detail::cCommandBytes) && account(s.streams.size(), Detail::cStreamBytes)
		&& account(s.events.size(), Detail::cEventBytes) && account(inAddCommand ? 1 : 0, Detail::cCommandBytes)
		&& account(inAddStream ? 1 : 0, Detail::cStreamBytes);
}

enum class CommandResult { Accepted, Invalid, Authority, Sequence, Late, Future, Full, Stopped };

// Late policy: reject every command at or before the last committed tick.
// A gap is allowed (packet loss); an older reordered sequence is rejected. Rejections
// never consume a sequence number. Acknowledgements distinguish accepted from applied.
inline CommandResult SubmitCommand(Checkpoint &ioState, const Limits &inLimits, const Command &inCommand)
{
	if (!ValidCommand(inCommand)) return CommandResult::Invalid;
	auto entity = std::find_if(ioState.entities.begin(), ioState.entities.end(), [&](const EntityState &e) { return e.id == inCommand.entity; });
	if (entity == ioState.entities.end() || entity->owner != inCommand.client) return CommandResult::Authority;
	if (inCommand.tick <= ioState.tick) return CommandResult::Late;
	if (inCommand.tick - ioState.tick > inLimits.future_ticks) return CommandResult::Future;
	auto stream = std::find_if(ioState.streams.begin(), ioState.streams.end(), [&](const StreamState &s) { return s.client == inCommand.client && s.entity == inCommand.entity; });
	if (stream != ioState.streams.end() && inCommand.sequence <= stream->accepted) return CommandResult::Sequence;
	if (ioState.commands.size() >= inLimits.commands || (stream == ioState.streams.end() && ioState.streams.size() >= inLimits.streams)
		|| !FitsStorage(ioState, inLimits, true, stream == ioState.streams.end())) return CommandResult::Full;
	// Sequence order must also follow intended time, otherwise a later ack could
	// acknowledge an earlier command that has not yet been executed.
	for (const Command &c : ioState.commands)
		if (c.client == inCommand.client && c.entity == inCommand.entity && c.tick > inCommand.tick) return CommandResult::Sequence;
	ioState.commands.push_back(inCommand);
	if (stream == ioState.streams.end()) ioState.streams.push_back({ inCommand.client, inCommand.entity, inCommand.sequence, 0 });
	else stream->accepted = inCommand.sequence;
	return CommandResult::Accepted;
}

inline bool ValidRotor(const std::array<float, 8> &r)
{
	if (!Finite(r)) return false;
	for (float v : r) if (std::abs(v) > 1.001f) return false;
	return Rotor(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]).IsNormalized(1.0e-4f);
}

inline bool ValidCheckpoint(const Checkpoint &s, const Limits &l)
{
	if (s.physics_hz == 0 || s.physics_hz > 1000000 || s.resource_period == 0 || s.snapshot_period == 0
		|| s.model_version == 0 || s.frame.id == 0 || s.frame.epoch > s.tick
		|| !Finite(s.frame.origin) || !Finite(s.frame.velocity)
		|| s.entities.size() > l.entities || s.commands.size() > l.commands || s.streams.size() > l.streams || s.events.size() > l.events) return false;
	if (!FitsStorage(s, l)) return false;
	for (size_t i = 0; i < s.entities.size(); ++i)
	{
		const EntityState &e = s.entities[i];
		if (e.id == 0 || e.geometry == 0 || e.frame != s.frame.id || !Finite(e.position) || !Finite(e.velocity)
			|| !Finite(e.center_of_mass) || !Finite(e.angular_velocity) || !Finite(e.inertia) || !Finite(e.actuators)
			|| !ValidRotor(e.orientation) || !ValidRotor(e.inertia_orientation)
			|| !std::isfinite(e.mass) || e.mass <= 0 || !std::isfinite(e.fuel) || e.fuel < 0
			|| !std::isfinite(e.energy) || e.energy < 0 || !std::isfinite(e.heat) || e.heat < 0
			|| !std::isfinite(e.health) || e.health < 0 || e.health > 1) return false;
		for (float v : e.inertia) if (v <= 0) return false;
		for (size_t j = 0; j < i; ++j) if (s.entities[j].id == e.id) return false;
	}
	for (size_t i = 0; i < s.streams.size(); ++i)
	{
		const StreamState &stream = s.streams[i];
		if (stream.client == 0 || stream.entity == 0 || stream.applied > stream.accepted) return false;
		for (size_t j = 0; j < i; ++j)
			if (s.streams[j].client == stream.client && s.streams[j].entity == stream.entity) return false;
		uint64 latest = stream.applied;
		for (const Command &c : s.commands)
			if (c.client == stream.client && c.entity == stream.entity) latest = std::max(latest, c.sequence);
		if (latest != stream.accepted) return false;
	}
	for (size_t i = 0; i < s.commands.size(); ++i)
	{
		const Command &c = s.commands[i];
		if (!ValidCommand(c) || c.tick <= s.tick || c.tick - s.tick > l.future_ticks) return false;
		auto e = std::find_if(s.entities.begin(), s.entities.end(), [&](const EntityState &v) { return v.id == c.entity && v.owner == c.client; });
		auto stream = std::find_if(s.streams.begin(), s.streams.end(), [&](const StreamState &v) { return v.entity == c.entity && v.client == c.client; });
		if (e == s.entities.end() || stream == s.streams.end() || c.sequence <= stream->applied || c.sequence > stream->accepted) return false;
		for (size_t j = 0; j < i; ++j)
		{
			const Command &prev = s.commands[j];
			if (prev.client == c.client && prev.entity == c.entity
				&& (prev.sequence == c.sequence || (prev.sequence < c.sequence && prev.tick > c.tick) || (prev.sequence > c.sequence && prev.tick < c.tick))) return false;
		}
	}
	for (size_t i = 0; i < s.events.size(); ++i)
	{
		const ImpactEvent &e = s.events[i];
		if (!ValidImpact(e) || e.tick <= s.tick || e.tick - s.tick > l.future_ticks || e.frame != s.frame.id) return false;
		if (std::none_of(s.entities.begin(), s.entities.end(), [&](const EntityState &v) { return v.id == e.target; })) return false;
		for (size_t j = 0; j < i; ++j) if (s.events[j].tick == e.tick && s.events[j].id == e.id) return false;
	}
	return true;
}

namespace Detail
{
// Fixed little-endian IEEE754 wire format. No native struct dumps, padding or
// unchecked length-to-allocation conversion. Byte limit and live allocation
// limit are independent. Caller supplies immutable model/shape assets separately.
struct Writer
{
	std::vector<uint8> bytes;
	size_t limit;
	bool ok = true;
	void Bits(uint64 value, size_t size)
	{
		if (!ok || bytes.size() > limit || size > limit - bytes.size()) { ok = false; return; }
		for (size_t i = 0; i < size; ++i) bytes.push_back(uint8(value >> (i * 8)));
	}
	void Value(uint64 v) { Bits(v, 8); }
	void Value(float v) { uint32 bits; std::memcpy(&bits, &v, 4); Bits(bits, 4); }
	void Value(double v) { uint64 bits; std::memcpy(&bits, &v, 8); Bits(bits, 8); }
	template <class T, size_t N> void Value(const std::array<T, N> &v) { for (T x : v) Value(x); }
	template <class... T> void operator () (const T &...v) { (Value(v), ...); }
	void Owner(ResponseOwner v) { Value(uint64(v)); }
	template <class T, class F> void Vector(const std::vector<T> &v, size_t, size_t, F fields)
	{
		Value(uint64(v.size()));
		for (const T &item : v) { if (!ok) break; fields(*this, item); }
	}
};

struct Reader
{
	const std::vector<uint8> &bytes;
	size_t budget, offset = 0;
	bool ok = true;
	uint64 Bits(size_t size)
	{
		if (!ok || size > bytes.size() - offset) { ok = false; return 0; }
		uint64 result = 0;
		for (size_t i = 0; i < size; ++i) result |= uint64(bytes[offset++]) << (8 * i);
		return result;
	}
	void Value(uint64 &v) { v = Bits(8); }
	void Value(float &v) { uint32 bits = uint32(Bits(4)); std::memcpy(&v, &bits, 4); }
	void Value(double &v) { uint64 bits = Bits(8); std::memcpy(&v, &bits, 8); }
	template <class T, size_t N> void Value(std::array<T, N> &v) { for (T &x : v) Value(x); }
	template <class... T> void operator () (T &...v) { (Value(v), ...); }
	void Owner(ResponseOwner &v) { uint64 value = Bits(8); if (value > 1) ok = false; else v = ResponseOwner(value); }
	template <class T, class F> void Vector(std::vector<T> &v, size_t limit, size_t wire_size, F fields)
	{
		uint64 count = Bits(8);
		if (!ok || count > limit || count > (bytes.size() - offset) / wire_size || count > budget / sizeof(T)) { ok = false; return; }
		budget -= size_t(count) * sizeof(T);
		v.resize(size_t(count));
		for (T &item : v) { if (!ok) break; fields(*this, item); }
	}
};

template <class A, class T> void EntityFields(A &a, T &v)
{
	a(v.id, v.owner, v.frame, v.topology_version, v.geometry, v.position, v.orientation,
		v.velocity, v.center_of_mass, v.angular_velocity, v.inertia, v.inertia_orientation, v.actuators,
		v.mass, v.fuel, v.energy, v.heat, v.health);
}
template <class A, class T> void CommandFields(A &a, T &v) { a(v.client, v.entity, v.sequence, v.tick, v.translation, v.rotation); }
template <class A, class T> void StreamFields(A &a, T &v) { a(v.client, v.entity, v.accepted, v.applied); }
template <class A, class T> void EventFields(A &a, T &v)
{
	a(v.tick, v.fraction, v.id, v.source, v.target, v.component, v.frame);
	a.Owner(v.owner);
	a(v.impulse, v.angular_impulse, v.energy);
}
template <class A, class T> void StateFields(A &a, T &s, const Limits &l)
{
	a(s.tick, s.random_seed, s.backlog, s.physics_hz, s.resource_period, s.snapshot_period, s.model_version,
		s.frame.id, s.frame.epoch, s.frame.origin, s.frame.velocity);
	a.Vector(s.entities, l.entities, cEntityBytes, [](auto &ar, auto &v) { EntityFields(ar, v); });
	a.Vector(s.commands, l.commands, cCommandBytes, [](auto &ar, auto &v) { CommandFields(ar, v); });
	a.Vector(s.streams, l.streams, cStreamBytes, [](auto &ar, auto &v) { StreamFields(ar, v); });
	a.Vector(s.events, l.events, cEventBytes, [](auto &ar, auto &v) { EventFields(ar, v); });
}
}

inline bool EncodeCheckpoint(const Checkpoint &inState, const Limits &inLimits, std::vector<uint8> &outBytes)
{
	static_assert(sizeof(float) == 4 && sizeof(double) == 8 && std::numeric_limits<float>::is_iec559 && std::numeric_limits<double>::is_iec559);
	if (!ValidCheckpoint(inState, inLimits)) return false;
	Detail::Writer writer { {}, inLimits.bytes };
	writer(uint64(0x5350414345415554ULL), uint64(1), uint64(0));
	Detail::StateFields(writer, inState, inLimits);
	if (!writer.ok) return false;
	uint64 length = uint64(writer.bytes.size());
	for (size_t i = 0; i < 8; ++i) writer.bytes[16 + i] = uint8(length >> (8 * i));
	outBytes = std::move(writer.bytes);
	return true;
}

inline bool DecodeCheckpoint(const std::vector<uint8> &inBytes, const Limits &inLimits, Checkpoint &outState)
{
	if (inBytes.size() > inLimits.bytes) return false;
	Detail::Reader reader { inBytes, inLimits.allocation_bytes };
	uint64 magic, version, length;
	reader(magic, version, length);
	if (!reader.ok || magic != 0x5350414345415554ULL || version != 1 || length != inBytes.size()) return false;
	Checkpoint draft;
	Detail::StateFields(reader, draft, inLimits);
	if (!reader.ok || reader.offset != inBytes.size() || !ValidCheckpoint(draft, inLimits)) return false;
	outState = std::move(draft);
	return true;
}

enum class Stage { Boundary, Region, Resources, Forces, Damage, Validate };
enum class Error { None, InvalidState, MissingCallback, Reservation, StageFailed, Physics, Impact, ClockOverflow, TickOverflow, Reentrant };

struct TickContext
{
	Tick tick;
	float seconds;
	bool resources_due, snapshot_due;
	const std::vector<Command> &commands;
};

struct Callbacks
{
	// Reserve body/fragment/ledger capacity BEFORE any boundary transaction.
	// If reservation fails, no stage or PhysicsSystem update runs. release(false)
	// must release reservations; it is not a rollback of external physics state.
	std::function<bool(const Checkpoint &, const TickContext &)> reserve;
	std::function<void(bool)> release;
	// Trusted adapters may mutate entities, frame, seed and pending events only.
	// Resource constraints run EVERY tick; resources_due schedules slower models.
	std::function<bool(Stage, Checkpoint &, const TickContext &)> stage;
	// Call PhysicsSystem::Update here exactly once, copy body state to the draft,
	// return its error bits (0 = success). There is no integration in this header.
	std::function<uint32(float, Checkpoint &)> physics_update;
	std::function<bool(Checkpoint &, const ImpactEvent &)> apply_query;
	// Solver-owned events NEVER call apply_query. This observer has no mutable
	// world parameter; damage/topology writes belong to the later Damage stage.
	std::function<bool(const ImpactEvent &)> observe;
	// Notification after the single commit. Must not throw or mutate coordinator.
	std::function<void(const Checkpoint &, bool)> publish;
};

struct AdvanceResult
{
	uint64 steps = 0, backlog_ticks = 0;
	bool overloaded = false;
	Error error = Error::None;
};

// Single-thread owner. Hand immutable shared snapshots to other threads using
// caller synchronization. Internally retained history is count- and per-image-
// bounded; callers retaining old handles own that additional memory. The view
// includes geometry IDs/poses/frame metadata, not a historical collision world.
// Query historical geometry through caller-owned immutable assets, never live bodies.
//
// Error policy is fail-stop: a failed Update may already have changed live Jolt
// bodies/caches. No snapshot/ack is published and Update is NEVER retried. Rebuild
// the external world and create a new coordinator explicitly; this is NOT rollback.
// Callbacks return errors, not exceptions. STL allocation exhaustion follows the
// embedding process policy (Jolt normally builds without C++ exceptions).
class Coordinator
{
public:
	using Snapshot = std::shared_ptr<const Checkpoint>;
	Coordinator(const Checkpoint &inInitial, const Limits &inLimits) : mLimits(inLimits)
	{
		std::vector<uint8> check;
		if (inLimits.history == 0 || inLimits.history > inLimits.allocation_bytes / sizeof(Snapshot)
			|| !EncodeCheckpoint(inInitial, inLimits, check)) { mError = Error::InvalidState; return; }
		mHistory.reserve(inLimits.history);
		mState = inInitial;
		mPublished = std::make_shared<const Checkpoint>(mState);
		mHistory.push_back(mPublished);
	}
	Coordinator(const Coordinator &) = delete;
	Coordinator &operator = (const Coordinator &) = delete;

	const Checkpoint &State() const { return mState; }
	Snapshot Published() const { return mPublished; }
	Snapshot Historical(Tick inTick) const
	{
		for (const Snapshot &s : mHistory) if (s->tick == inTick) return s;
		return {};
	}
	template <class Query> bool QueryHistory(Tick inTick, ID inFrame, Query &&inQuery) const
	{
		Snapshot s = Historical(inTick);
		if (!s || s->frame.id != inFrame) return false;
		inQuery(*s);
		return true;
	}
	bool Stopped() const { return mError != Error::None; }
	Error Failure() const { return mError; }
	CommandResult Submit(const Command &inCommand)
	{
		if (Stopped() || mRunning) return CommandResult::Stopped;
		return SubmitCommand(mState, mLimits, inCommand);
	}

	// Time debt is exact integer nanoseconds * physics_hz. One step costs 1e9
	// units, so 120 Hz never truncates to 8,333,333 ns. max_steps only limits work
	// per call: leftover debt is kept, visible, and can be drained with elapsed=0.
	AdvanceResult Advance(uint64 inElapsedNanoseconds, uint64 inMaxSteps, const Callbacks &inCallbacks)
	{
		AdvanceResult result;
		if (mRunning) mError = Error::Reentrant;
		if (!Stopped())
		{
			if (inElapsedNanoseconds > (std::numeric_limits<uint64>::max() - mState.backlog) / mState.physics_hz) mError = Error::ClockOverflow;
			else
			{
				mState.backlog += inElapsedNanoseconds * mState.physics_hz;
				mRunning = true;
				while (mState.backlog >= cTickCost && result.steps < inMaxSteps && !Stopped())
				{
					if (!Step(inCallbacks)) break;
					++result.steps;
				}
				mRunning = false;
			}
		}
		result.error = mError;
		result.backlog_ticks = mState.backlog / cTickCost;
		result.overloaded = result.backlog_ticks != 0;
		return result;
	}

private:
	bool Fail(Error inError) { if (!Stopped()) mError = inError; return false; }
	bool Step(const Callbacks &cb)
	{
		if (mState.tick == std::numeric_limits<Tick>::max()) return Fail(Error::TickOverflow);
		if (!cb.reserve || !cb.physics_update) return Fail(Error::MissingCallback);
		std::vector<Command> commands;
		for (const Command &c : mState.commands) if (c.tick == mState.tick + 1) commands.push_back(c);
		std::sort(commands.begin(), commands.end(), [](const Command &a, const Command &b) { return std::tie(a.client, a.entity, a.sequence) < std::tie(b.client, b.entity, b.sequence); });
		TickContext context { mState.tick + 1, 1.0f / float(mState.physics_hz), (mState.tick + 1) % mState.resource_period == 0, (mState.tick + 1) % mState.snapshot_period == 0, commands };
		// Allocate the unpublished image before external side effects.
		auto draft = std::make_shared<Checkpoint>(mState);
		bool success = false;
		struct Release
		{
			const Callbacks &cb;
			bool &success;
			~Release() { if (cb.release) cb.release(success); }
		} release { cb, success };
		if (!cb.reserve(mState, context) || Stopped()) return Fail(Error::Reservation);
		auto stage = [&](Stage s) { return !Stopped() && (!cb.stage || cb.stage(s, *draft, context)); };
		for (Stage s : { Stage::Boundary, Stage::Region, Stage::Resources, Stage::Forces })
			if (!stage(s)) return Fail(Error::StageFailed);
		if (!ValidCheckpoint(*draft, mLimits)) return Fail(Error::InvalidState);
		if (cb.physics_update(context.seconds, *draft) != 0 || Stopped()) return Fail(Error::Physics);
		// Validate the entire batch before delivering any event. A query response
		// here takes effect at the end boundary, not at the recorded fractional TOI.
		if (!ValidCheckpoint(*draft, mLimits)) return Fail(Error::InvalidState);
		std::vector<ImpactEvent> impacts;
		for (const ImpactEvent &e : draft->events) if (e.tick == context.tick) impacts.push_back(e);
		ImpactBatch batch(impacts, mLimits);
		// Remove ONLY the captured batch before callbacks. A callback emitting a
		// new same-tick event is rejected by final validation, never silently lost.
		// Such chains require a real bounded event/substep loop, not another Update.
		draft->events.erase(std::remove_if(draft->events.begin(), draft->events.end(), [&](const ImpactEvent &e) { return e.tick == context.tick; }), draft->events.end());
		if (!batch.Process([&](const ImpactEvent &e) { return cb.apply_query && cb.apply_query(*draft, e); },
			[&](const ImpactEvent &e) { return !cb.observe || cb.observe(e); }) || Stopped()) return Fail(Error::Impact);
		if (!stage(Stage::Damage) || !stage(Stage::Validate)) return Fail(Error::StageFailed);
		if (draft->tick != mState.tick || draft->backlog != mState.backlog || draft->physics_hz != mState.physics_hz
			|| draft->resource_period != mState.resource_period || draft->snapshot_period != mState.snapshot_period) return Fail(Error::InvalidState);
		for (const Command &c : commands)
			for (StreamState &s : draft->streams)
				if (s.client == c.client && s.entity == c.entity) s.applied = std::max(s.applied, c.sequence);
		draft->commands.erase(std::remove_if(draft->commands.begin(), draft->commands.end(), [&](const Command &c) { return c.tick == context.tick; }), draft->commands.end());
		draft->tick = context.tick;
		draft->backlog -= cTickCost;
		std::sort(draft->entities.begin(), draft->entities.end(), [](const EntityState &a, const EntityState &b) { return a.id < b.id; });
		std::sort(draft->streams.begin(), draft->streams.end(), [](const StreamState &a, const StreamState &b) { return std::tie(a.client, a.entity) < std::tie(b.client, b.entity); });
		std::sort(draft->commands.begin(), draft->commands.end(), [](const Command &a, const Command &b) { return std::tie(a.tick, a.client, a.entity, a.sequence) < std::tie(b.tick, b.client, b.entity, b.sequence); });
		std::sort(draft->events.begin(), draft->events.end(), ImpactLess);
		std::vector<uint8> check;
		if (!EncodeCheckpoint(*draft, mLimits, check) || Stopped()) return Fail(Error::InvalidState);
		// Everything observable changes at this one publication boundary.
		mState = *draft;
		mPublished = std::move(draft);
		if (mHistory.size() == mLimits.history) mHistory.erase(mHistory.begin());
		mHistory.push_back(mPublished);
		success = true;
		if (cb.publish) cb.publish(*mPublished, context.snapshot_due);
		return true;
	}
	static constexpr uint64 cTickCost = 1000000000;
	Limits mLimits;
	Checkpoint mState;
	Snapshot mPublished;
	std::vector<Snapshot> mHistory;
	Error mError = Error::None;
	bool mRunning = false;
};

enum class ReconcileResult { Applied, Stale, Invalid, TooFar, ReplayFailed };

// Local-entity prediction only. The callback supplies the matching model and
// advances one full fixed tick, including ticks with no input. This does not
// rewind the coupled world, resources of other ships or solver caches. Pending
// inputs at/before the server cutoff are discarded under the server late policy.
class PredictionBuffer
{
public:
	PredictionBuffer(ID inClient, ID inEntity, uint64 inModelVersion, const Limits &inLimits)
		: mClient(inClient), mEntity(inEntity), mModelVersion(inModelVersion), mLimits(inLimits) { }
	bool Queue(const Command &c)
	{
		if (!ValidCommand(c) || c.client != mClient || c.entity != mEntity || c.sequence <= mLastSequence
			|| c.tick <= mServerTick || c.tick < mLastIntendedTick || c.tick - mServerTick > mLimits.future_ticks
			|| mPending.size() >= mLimits.commands || mPending.size() >= mLimits.allocation_bytes / sizeof(Command)) return false;
		mPending.push_back(c);
		mLastSequence = c.sequence;
		mLastIntendedTick = c.tick;
		return true;
	}
	const std::vector<Command> &Pending() const { return mPending; }
	const EntityState &Predicted() const { return mPredicted; }
	Tick PredictedTick() const { return mPredictedTick; }
	size_t LateDiscarded() const { return mLateDiscarded; }

	template <class Replay> ReconcileResult Reconcile(const Checkpoint &server, Tick inTargetTick, Replay &&inReplay)
	{
		if (mReady && server.tick <= mServerTick) return ReconcileResult::Stale;
		if (!ValidCheckpoint(server, mLimits) || server.model_version != mModelVersion) return ReconcileResult::Invalid;
		if (inTargetTick < server.tick || inTargetTick - server.tick > mLimits.replay_ticks) return ReconcileResult::TooFar;
		auto entity = std::find_if(server.entities.begin(), server.entities.end(), [&](const EntityState &e) { return e.id == mEntity && e.owner == mClient; });
		if (entity == server.entities.end()) return ReconcileResult::Invalid;
		uint64 applied = 0, accepted = 0;
		for (const StreamState &s : server.streams)
			if (s.entity == mEntity && s.client == mClient) { applied = s.applied; accepted = s.accepted; }
		std::vector<Command> pending;
		size_t late = 0;
		for (const Command &c : mPending)
			if (c.sequence > applied)
			{
				if (c.tick <= server.tick) ++late;
				else pending.push_back(c);
			}
		EntityState predicted = *entity;
		// Count up by a bounded delta, not tick <= target, which wraps at UINT64_MAX.
		for (Tick delta = 0; delta < inTargetTick - server.tick; ++delta)
		{
			Tick tick = server.tick + delta + 1;
			std::vector<Command> commands;
			for (const Command &c : pending) if (c.tick == tick) commands.push_back(c);
			if (!inReplay(predicted, tick, 1.0f / float(server.physics_hz), commands)) return ReconcileResult::ReplayFailed;
		}
		Checkpoint probe;
		probe.tick = inTargetTick;
		probe.frame = server.frame;
		probe.entities.push_back(predicted);
		if (predicted.id != mEntity || predicted.owner != mClient || !ValidCheckpoint(probe, mLimits)) return ReconcileResult::ReplayFailed;
		mPredicted = predicted;
		mPredictedTick = inTargetTick;
		mPending = std::move(pending);
		mServerTick = server.tick;
		mLastSequence = std::max(mLastSequence, accepted);
		mLateDiscarded = late;
		mReady = true;
		return ReconcileResult::Applied;
	}
private:
	ID mClient, mEntity;
	uint64 mModelVersion;
	Limits mLimits;
	std::vector<Command> mPending;
	EntityState mPredicted;
	Tick mServerTick = 0, mPredictedTick = 0, mLastIntendedTick = 0;
	uint64 mLastSequence = 0;
	size_t mLateDiscarded = 0;
	bool mReady = false;
};
}
