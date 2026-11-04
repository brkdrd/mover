#include <UnitTestFramework.h>
#include <Jolt/SpaceSim/Authority.h>
#include <PhysicsTestContext.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <limits>

using namespace JPH;
using namespace JPH::SpaceSim;

#define AUTH_REQUIRE(...) do { bool auth_ok = (__VA_ARGS__); CHECK(auth_ok); if (!auth_ok) return; } while (false)

TEST_SUITE("SpaceSimAuthority")
{
	TEST_CASE("PhysicsSystem adapter advances a real body exactly once using four plus six axis requests")
	{
		PhysicsTestContext physics;
		physics.ZeroGravity();
		BodyCreationSettings settings(new BoxShape(Vec4::sReplicate(0.5f), 0.0f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
		settings.mLinearDamping = 0; settings.mAngularDamping = 0;
		settings.mOverrideMassProperties = EOverrideMassProperties::CalculateInertia;
		settings.mMassPropertiesOverride.mMass = 1;
		Body &body = physics.CreateBody(settings, EActivation::Activate);
		Checkpoint initial;
		EntityState ship; ship.id = 7; ship.owner = 3; ship.fuel = 1;
		for (uint i = 0; i < 6; ++i) ship.inertia[i] = 1.0f / body.GetMotionProperties()->GetInverseInertiaDiagonal().GetVec8()[i];
		initial.entities.push_back(ship);
		Coordinator c(initial, Limits {});
		Command cmd; cmd.client = 3; cmd.entity = 7; cmd.sequence = 1; cmd.tick = 1;
		cmd.translation.fill(0.25f); cmd.rotation.fill(0.125f);
		AUTH_REQUIRE(c.Submit(cmd) == CommandResult::Accepted);
		uint32 updates = 0, reservations = 0, released = 0;
		Callbacks cb;
		cb.reserve = [&](const Checkpoint &, const TickContext &) { ++reservations; return true; };
		cb.release = [&](bool committed) { CHECK(committed); ++released; };
		cb.stage = [&](Stage stage, Checkpoint &s, const TickContext &ctx)
		{
			if (stage == Stage::Forces)
				for (const Command &request : ctx.commands)
				{
					const auto &f = request.translation; const auto &t = request.rotation;
					body.AddForce(Vec4(f[0], f[1], f[2], f[3]));
					body.AddTorque(Bivec(t[0], t[1], t[2], t[3], t[4], t[5]));
					s.entities[0].fuel -= ctx.seconds;
				}
			return true;
		};
		cb.physics_update = [&](float dt, Checkpoint &s)
		{
			++updates;
			uint32 error = uint32(physics.GetSystem()->Update(dt, 1, physics.GetTempAllocator(), physics.GetJobSystem()));
			EntityState &e = s.entities[0];
			for (uint i = 0; i < 4; ++i) { e.position[i] = body.GetPosition()[i]; e.velocity[i] = body.GetLinearVelocity()[i]; }
			for (uint i = 0; i < 6; ++i) e.angular_velocity[i] = body.GetAngularVelocity().GetVec8()[i];
			for (uint i = 0; i < 8; ++i) { e.orientation[i] = body.GetRotation().GetVec8()[i]; e.inertia_orientation[i] = body.GetMotionProperties()->GetInertiaRotation().GetVec8()[i]; }
			return error;
		};
		AUTH_REQUIRE(c.Advance(10000000, 1, cb).steps == 1);
		CHECK(updates == 1); CHECK(reservations == 1); CHECK(released == 1);
		for (uint i = 0; i < 4; ++i)
		{
			CHECK_APPROX_EQUAL(c.Published()->entities[0].velocity[i], 0.25f / 120.0f);
			CHECK(c.Published()->entities[0].position[i] == double(body.GetPosition()[i]));
		}
		for (float w : c.Published()->entities[0].angular_velocity) CHECK(w > 0);
		CHECK_APPROX_EQUAL(c.Published()->entities[0].fuel, 1.0f - 1.0f / 120.0f);
		// No second integration or duplicate force on the next tick.
		AUTH_REQUIRE(c.Advance(10000000, 1, cb).steps == 1);
		CHECK(updates == 2);
		for (uint i = 0; i < 4; ++i) CHECK_APPROX_EQUAL(c.Published()->entities[0].velocity[i], 0.25f / 120.0f);
	}

	TEST_CASE("Seeded pending events and commands resume from a checkpoint in deterministic order")
	{
		std::vector<uint8> reference;
		for (bool reversed : { false, true })
		{
			Checkpoint initial; EntityState ship; ship.id = 7; ship.owner = 3; ship.fuel = 10; initial.entities.push_back(ship);
			for (uint64 i = 1; i <= 3; ++i)
			{
				Command cmd; cmd.client = 3; cmd.entity = 7; cmd.sequence = i; cmd.tick = i; cmd.translation[3] = 0.25f;
				AUTH_REQUIRE(SubmitCommand(initial, Limits {}, cmd) == CommandResult::Accepted);
			}
			for (uint64 i : { 9ULL, 8ULL, 7ULL })
			{
				ImpactEvent impact; impact.tick = i == 7 ? 3 : 1; impact.id = i; impact.source = 7; impact.target = 7;
				impact.owner = i == 8 ? ResponseOwner::Solver : ResponseOwner::Query;
				impact.impulse[3] = 0.5f;
				initial.events.push_back(impact);
			}
			if (reversed) std::reverse(initial.events.begin(), initial.events.end());
			std::vector<ID> order;
			uint32 applied = 0;
			Callbacks cb;
			cb.reserve = [](const Checkpoint &, const TickContext &) { return true; };
			cb.stage = [](Stage stage, Checkpoint &s, const TickContext &ctx)
			{
				if (stage == Stage::Boundary)
				{
					s.random_seed = s.random_seed * 6364136223846793005ULL + 1;
					for (const Command &cmd : ctx.commands) s.entities[0].actuators[3] = cmd.translation[3];
				}
				if (stage == Stage::Resources) s.entities[0].fuel -= 0.25f;
				return true;
			};
			cb.physics_update = [](float dt, Checkpoint &s) { s.entities[0].position[3] += s.entities[0].velocity[3] * dt; return uint32(0); };
			cb.apply_query = [&](Checkpoint &s, const ImpactEvent &e) { ++applied; s.entities[0].velocity[3] += e.impulse[3]; return true; };
			cb.observe = [&](const ImpactEvent &e) { order.push_back(e.id); return true; };
			Coordinator uninterrupted(initial, Limits {});
			AUTH_REQUIRE(uninterrupted.Advance(25000000, 1, cb).steps == 1);
			CHECK(order == std::vector<ID> { 8, 9 });
			CHECK(applied == 1);
			std::vector<uint8> checkpoint;
			AUTH_REQUIRE(EncodeCheckpoint(uninterrupted.State(), Limits {}, checkpoint));
			Checkpoint restored;
			AUTH_REQUIRE(DecodeCheckpoint(checkpoint, Limits {}, restored));
			Coordinator resumed(restored, Limits {});
			AUTH_REQUIRE(uninterrupted.Advance(0, 8, cb).steps == 2);
			std::vector<uint8> expected;
			AUTH_REQUIRE(EncodeCheckpoint(uninterrupted.State(), Limits {}, expected));
			AUTH_REQUIRE(resumed.Advance(0, 8, cb).steps == 2);
			std::vector<uint8> actual;
			AUTH_REQUIRE(EncodeCheckpoint(resumed.State(), Limits {}, actual));
			CHECK(actual == expected);
			CHECK(resumed.Published()->entities[0].fuel == 9.25f);
			CHECK(resumed.Published()->entities[0].velocity[3] == 1.0f);
			CHECK(resumed.Published()->streams[0].applied == 3);
			if (reversed) CHECK(reference == actual); else reference = actual;
		}
	}

	TEST_CASE("Impact callbacks cannot silently erase newly emitted same tick events")
	{
		Checkpoint state; EntityState ship; ship.id = 7; ship.owner = 3; state.entities.push_back(ship);
		ImpactEvent impact; impact.tick = 1; impact.id = 1; impact.source = 7; impact.target = 7; impact.owner = ResponseOwner::Query;
		state.events.push_back(impact);
		Coordinator c(state, Limits {});
		Callbacks cb;
		cb.reserve = [](const Checkpoint &, const TickContext &) { return true; };
		cb.physics_update = [](float, Checkpoint &) { return uint32(0); };
		cb.apply_query = [&](Checkpoint &s, const ImpactEvent &) { ImpactEvent chained = impact; chained.id = 2; s.events.push_back(chained); return true; };
		CHECK(c.Advance(10000000, 1, cb).error == Error::InvalidState);
		CHECK(c.Published()->tick == 0);
		CHECK(c.State().events.size() == 1);
	}

	TEST_CASE("Impact callback failure cannot retry a response or publish partial damage")
	{
		Checkpoint state; EntityState ship; ship.id = 7; ship.owner = 3; state.entities.push_back(ship);
		for (ID id : { 1ULL, 2ULL }) { ImpactEvent e; e.tick = 1; e.id = id; e.source = 7; e.target = 7; e.owner = ResponseOwner::Query; state.events.push_back(e); }
		Coordinator c(state, Limits {});
		uint32 applied = 0, released = 0;
		Callbacks cb;
		cb.reserve = [](const Checkpoint &, const TickContext &) { return true; };
		cb.release = [&](bool committed) { CHECK_FALSE(committed); ++released; };
		cb.physics_update = [](float, Checkpoint &) { return uint32(0); };
		cb.apply_query = [&](Checkpoint &s, const ImpactEvent &) { ++applied; s.entities[0].velocity[3] += 1; return true; };
		cb.observe = [](const ImpactEvent &) { return false; };
		CHECK(c.Advance(10000000, 1, cb).error == Error::Impact);
		CHECK(c.Published()->tick == 0);
		CHECK(c.Published()->entities[0].velocity[3] == 0);
		CHECK(c.State().events.size() == 2);
		CHECK(c.Advance(0, 1, cb).error == Error::Impact);
		CHECK(applied == 1); CHECK(released == 1);
	}

	TEST_CASE("Clock partitioning is exact and overload or integer overflow never silently drops time")
	{
		Checkpoint initial;
		Callbacks cb; cb.reserve = [](const Checkpoint &, const TickContext &) { return true; };
		cb.physics_update = [](float, Checkpoint &) { return uint32(0); };
		Coordinator whole(initial, Limits {}), split(initial, Limits {});
		AUTH_REQUIRE(whole.Advance(1000000000, 120, cb).steps == 120);
		for (uint32 i = 0; i < 1000; ++i) CHECK(split.Advance(1000000, 1, cb).error == Error::None);
		CHECK(split.State().tick == 120);
		CHECK(split.State().backlog == 0);
		std::vector<uint8> a, b;
		AUTH_REQUIRE(EncodeCheckpoint(whole.State(), Limits {}, a));
		AUTH_REQUIRE(EncodeCheckpoint(split.State(), Limits {}, b));
		CHECK(a == b);
		Coordinator overflow(initial, Limits {});
		CHECK(overflow.Advance(std::numeric_limits<uint64>::max(), 1, cb).error == Error::ClockOverflow);
		CHECK(overflow.State().backlog == 0);
		initial.tick = std::numeric_limits<uint64>::max();
		Coordinator tick_overflow(initial, Limits {});
		CHECK(tick_overflow.Advance(10000000, 1, cb).error == Error::TickOverflow);
		CHECK(tick_overflow.Published()->tick == initial.tick);
		CHECK(tick_overflow.State().backlog == 1200000000ULL);
		Coordinator missing(Checkpoint {}, Limits {});
		CHECK(missing.Advance(10000000, 1, Callbacks {}).error == Error::MissingCallback);
	}


	TEST_CASE("Untrusted snapshot counts and nonfinite bit patterns are rejected before allocation or arithmetic")
	{
		Checkpoint state; EntityState ship; ship.id = 7; ship.owner = 3; state.entities.push_back(ship);
		std::vector<uint8> bytes;
		AUTH_REQUIRE(EncodeCheckpoint(state, Limits {}, bytes));
		Checkpoint output; output.tick = 999;
		// v1: fixed prelude is 160 bytes, followed by entity count and fields.
		auto bad = bytes;
		std::fill(bad.begin() + 160, bad.begin() + 168, uint8(255));
		CHECK_FALSE(DecodeCheckpoint(bad, Limits {}, output));
		for (uint32 bits : { 0x7f800000U, 0x7fc00000U, 0x7f800001U })
		{
			bad = bytes;
			for (size_t i = 0; i < 4; ++i) bad[424 + i] = uint8(bits >> (8 * i));
			CHECK_FALSE(DecodeCheckpoint(bad, Limits {}, output));
			Command cmd; cmd.client = 3; cmd.entity = 7; cmd.sequence = 1; cmd.tick = 1;
			std::memcpy(&cmd.rotation[5], &bits, 4);
			CHECK(SubmitCommand(state, Limits {}, cmd) == CommandResult::Invalid);
		}
		bad = bytes;
		uint64 signaling = 0x7ff0000000000001ULL;
		for (size_t i = 0; i < 8; ++i) bad[208 + i] = uint8(signaling >> (8 * i));
		CHECK_FALSE(DecodeCheckpoint(bad, Limits {}, output));
		// A finite but non-Spin(4) orientation is invalid as well.
		bad = bytes;
		std::fill(bad.begin() + 240, bad.begin() + 272, uint8(0));
		CHECK_FALSE(DecodeCheckpoint(bad, Limits {}, output));
		CHECK(output.tick == 999);
	}

	TEST_CASE("Command admission cannot accept an uncheckpointable future queue")
	{
		Checkpoint s; EntityState e; e.id = 7; e.owner = 3; s.entities.push_back(e);
		std::vector<uint8> encoded;
		AUTH_REQUIRE(EncodeCheckpoint(s, Limits {}, encoded));
		Limits l; l.bytes = encoded.size();
		Command cmd; cmd.client = 3; cmd.entity = 7; cmd.sequence = 1; cmd.tick = 2;
		CHECK(SubmitCommand(s, l, cmd) == CommandResult::Full);
		CHECK(s.commands.empty()); CHECK(s.streams.empty());
	}

	TEST_CASE("Command admission respects the total causal allocation budget")
	{
		Checkpoint s; EntityState e; e.id = 7; e.owner = 3; s.entities.push_back(e);
		Limits l; l.allocation_bytes = sizeof(EntityState) + sizeof(Command) + sizeof(StreamState) - 1;
		Command cmd; cmd.client = 3; cmd.entity = 7; cmd.sequence = 1; cmd.tick = 1;
		CHECK(SubmitCommand(s, l, cmd) == CommandResult::Full);
		CHECK(s.commands.empty());
		CHECK(s.streams.empty());
		++l.allocation_bytes;
		CHECK(SubmitCommand(s, l, cmd) == CommandResult::Accepted);
	}

	TEST_CASE("Bounded reconciliation handles loss reorder duplicate snapshots and replay failure")
	{
		Checkpoint initial; EntityState ship; ship.id = 7; ship.owner = 3; initial.entities.push_back(ship);
		Limits limits; limits.commands = 6; limits.replay_ticks = 6;
		Coordinator server(initial, limits);
		PredictionBuffer client(3, 7, 1, limits);
		std::vector<Command> packets;
		for (uint64 i = 1; i <= 6; ++i)
		{
			Command cmd; cmd.client = 3; cmd.entity = 7; cmd.sequence = i; cmd.tick = i; cmd.translation[3] = 0.25f;
			packets.push_back(cmd);
			CHECK(client.Queue(cmd));
		}
		CHECK_FALSE(client.Queue(packets.back()));
		// Bounded channel: command 4 lost; 2 arrives after 3; 5 duplicated.
		CHECK(server.Submit(packets[0]) == CommandResult::Accepted);
		CHECK(server.Submit(packets[2]) == CommandResult::Accepted);
		CHECK(server.Submit(packets[1]) == CommandResult::Sequence);
		CHECK(server.Submit(packets[4]) == CommandResult::Accepted);
		CHECK(server.Submit(packets[4]) == CommandResult::Sequence);
		CHECK(server.Submit(packets[5]) == CommandResult::Accepted);
		Callbacks cb;
		cb.reserve = [](const Checkpoint &, const TickContext &) { return true; };
		cb.stage = [](Stage stage, Checkpoint &s, const TickContext &ctx)
		{
			if (stage == Stage::Boundary) for (const Command &cmd : ctx.commands) s.entities[0].velocity[3] = cmd.translation[3];
			return true;
		};
		cb.physics_update = [](float dt, Checkpoint &s) { s.entities[0].position[3] += s.entities[0].velocity[3] * dt; return uint32(0); };
		AUTH_REQUIRE(server.Advance(50000000, 4, cb).steps == 4);
		auto old = server.Historical(2);
		auto snapshot = server.Published();
		auto replay = [](EntityState &e, Tick, float dt, const std::vector<Command> &commands)
		{
			for (const Command &cmd : commands) e.velocity[3] = cmd.translation[3];
			e.position[3] += e.velocity[3] * dt;
			return true;
		};
		CHECK(client.Reconcile(*snapshot, 6, replay) == ReconcileResult::Applied);
		CHECK(client.Pending().size() == 2);
		CHECK(client.LateDiscarded() == 1);
		CHECK(client.Reconcile(*old, 6, replay) == ReconcileResult::Stale);
		CHECK(client.Reconcile(*snapshot, 6, replay) == ReconcileResult::Stale);
		AUTH_REQUIRE(server.Advance(0, 2, cb).steps == 2);
		CHECK(client.Predicted().position[3] == server.Published()->entities[0].position[3]);
		EntityState before = client.Predicted();
		CHECK(client.Reconcile(*server.Published(), 13, replay) == ReconcileResult::TooFar);
		CHECK(client.Reconcile(*server.Published(), 7, [](EntityState &e, Tick, float, const std::vector<Command> &) { e.position[3] = 999; return false; }) == ReconcileResult::ReplayFailed);
		CHECK(client.Predicted().position == before.position);
		CHECK(client.Pending().size() == 2);
		CHECK(client.Reconcile(*server.Published(), 6, replay) == ReconcileResult::Applied);
		CHECK(client.Pending().empty());
		CHECK(server.Submit(packets[3]) == CommandResult::Late);
		uint32 queries = 0;
		CHECK(server.QueryHistory(5, 1, [&](const Checkpoint &s) { ++queries; CHECK(s.tick == 5); CHECK(s.entities[0].geometry == 1); }));
		CHECK_FALSE(server.QueryHistory(5, 99, [&](const Checkpoint &) { ++queries; }));
		CHECK_FALSE(server.QueryHistory(99, 1, [&](const Checkpoint &) { ++queries; }));
		CHECK(queries == 1);
	}

	TEST_CASE("Coordinator preserves fixed time backlog ordered stages and atomic publication")
	{
		Checkpoint initial;
		initial.resource_period = 2; initial.snapshot_period = 3;
		EntityState ship; ship.id = 7; ship.owner = 3;
		initial.entities.push_back(ship);
		Limits limits; limits.history = 2;
		Coordinator coordinator(initial, limits);
		std::vector<Stage> stages;
		uint32 updates = 0, resource_updates = 0, network_updates = 0, reservations = 0;
		Callbacks callbacks;
		callbacks.reserve = [&](const Checkpoint &, const TickContext &) { ++reservations; return true; };
		callbacks.stage = [&](Stage stage, Checkpoint &draft, const TickContext &ctx)
		{
			stages.push_back(stage);
			CHECK(coordinator.Published()->tick == ctx.tick - 1);
			if (stage == Stage::Resources && ctx.resources_due) ++resource_updates;
			if (stage == Stage::Boundary)
				for (const Command &c : ctx.commands) draft.entities[0].actuators[9] = c.rotation[5];
			return true;
		};
		callbacks.physics_update = [&](float dt, Checkpoint &draft) { CHECK(dt == 1.0f / 120.0f); ++updates; draft.entities[0].position[3] += 1; return uint32(0); };
		Command cmd; cmd.entity = 7; cmd.client = 3; cmd.sequence = 1; cmd.tick = 1; cmd.rotation[5] = 0.5f;
		AUTH_REQUIRE(coordinator.Submit(cmd) == CommandResult::Accepted);
		auto first = coordinator.Advance(25000000, 1, callbacks);
		CHECK(first.steps == 1);
		CHECK(first.overloaded);
		CHECK(coordinator.State().backlog == 2000000000ULL);
		CHECK(coordinator.Published()->tick == 1);
		CHECK(coordinator.Published()->entities[0].actuators[9] == 0.5f);
		CHECK(coordinator.Published()->streams[0].applied == 1);
		CHECK(coordinator.Published()->commands.empty());
		CHECK(stages == std::vector<Stage> { Stage::Boundary, Stage::Region, Stage::Resources, Stage::Forces, Stage::Damage, Stage::Validate });
		auto held = coordinator.Historical(1);
		AUTH_REQUIRE(held != nullptr);
		callbacks.publish = [&](const Checkpoint &s, bool network_due) { CHECK(s.tick > 1); if (network_due) ++network_updates; };
		auto rest = coordinator.Advance(0, 8, callbacks);
		CHECK(rest.steps == 2);
		CHECK_FALSE(rest.overloaded);
		CHECK(updates == 3);
		CHECK(reservations == 3);
		CHECK(resource_updates == 1);
		CHECK(network_updates == 1);
		CHECK(coordinator.State().backlog == 0);
		CHECK(coordinator.Published()->tick == 3);
		CHECK(coordinator.Historical(1) == nullptr);
		CHECK(coordinator.Historical(2) != nullptr);
		CHECK(held->entities[0].position[3] == 1);
		CHECK(coordinator.Published()->entities[0].position[3] == 3);
		CHECK(coordinator.Submit(cmd) == CommandResult::Late);
	}

	TEST_CASE("Coordinator failure stops permanently without half publishing or acknowledging")
	{
		for (int fail = 0; fail < 5; ++fail)
		{
			Checkpoint initial; EntityState ship; ship.id = 7; ship.owner = 3; initial.entities.push_back(ship);
			Coordinator c(initial, Limits {});
			Command cmd; cmd.client = 3; cmd.entity = 7; cmd.tick = 1; cmd.sequence = 1;
			AUTH_REQUIRE(c.Submit(cmd) == CommandResult::Accepted);
			auto before = c.Published();
			uint32 updates = 0, published = 0;
			Callbacks cb;
			cb.reserve = [&](const Checkpoint &, const TickContext &) { return fail != 0; };
			cb.stage = [&](Stage stage, Checkpoint &s, const TickContext &)
			{
				s.entities[0].topology_version = 99;
				if (fail == 4 && stage == Stage::Validate) s.entities[0].fuel = -1;
				return !(fail == 1 && stage == Stage::Forces) && !(fail == 3 && stage == Stage::Damage);
			};
			cb.physics_update = [&](float, Checkpoint &) { ++updates; return uint32(fail == 2 ? 1 : 0); };
			cb.publish = [&](const Checkpoint &, bool) { ++published; };
			auto result = c.Advance(10000000, 1, cb);
			CHECK(result.error != Error::None);
			CHECK(c.Stopped());
			CHECK(c.Published() == before);
			CHECK(c.State().tick == 0);
			CHECK(c.State().entities[0].topology_version == 0);
			CHECK(c.State().streams[0].applied == 0);
			CHECK(c.State().commands.size() == 1);
			CHECK(c.State().backlog == 1200000000ULL);
			CHECK(published == 0);
			CHECK(c.Advance(0, 1, cb).steps == 0);
			CHECK(updates == uint32(fail < 2 ? 0 : 1));
			CHECK(c.Submit(cmd) == CommandResult::Stopped);
		}
	}

	TEST_CASE("Checkpoint codec preserves causal state and rejects malformed input atomically")
	{
		Checkpoint state;
		state.tick = 11; state.random_seed = 123; state.backlog = 456;
		state.physics_hz = 120; state.resource_period = 4; state.snapshot_period = 6;
		state.frame.id = 5; state.frame.epoch = 11;
		state.frame.origin[3] = 123456789.25;
		state.frame.velocity[0] = 42;
		EntityState ship;
		ship.id = 7; ship.owner = 3; ship.frame = 5; ship.topology_version = 9; ship.geometry = 10;
		ship.position[3] = 1.25; ship.velocity[3] = 2; ship.angular_velocity[5] = 3;
		ship.mass = 8; ship.inertia[5] = 4; ship.fuel = 5; ship.energy = 6; ship.heat = 7;
		ship.health = 0.75f; ship.actuators[9] = 0.5f; ship.center_of_mass[3] = 0.25f;
		state.entities.push_back(ship);
		Command cmd; cmd.entity = 7; cmd.client = 3; cmd.sequence = 2; cmd.tick = 12;
		AUTH_REQUIRE(SubmitCommand(state, Limits {}, cmd) == CommandResult::Accepted);
		state.streams[0].applied = 1;
		ImpactEvent impact; impact.id = 1; impact.tick = 13; impact.source = 7; impact.target = 7; impact.frame = 5;
		state.events.push_back(impact);
		std::vector<uint8> bytes;
		AUTH_REQUIRE(EncodeCheckpoint(state, Limits {}, bytes));
		Checkpoint decoded;
		AUTH_REQUIRE(DecodeCheckpoint(bytes, Limits {}, decoded));
		std::vector<uint8> again;
		AUTH_REQUIRE(EncodeCheckpoint(decoded, Limits {}, again));
		CHECK(bytes == again);
		CHECK(decoded.entities[0].actuators[9] == 0.5f);
		CHECK(decoded.frame.origin[3] == state.frame.origin[3]);
		CHECK(decoded.streams[0].applied == 1);
		CHECK(decoded.events[0].tick == 13);
		// Every truncation, trailing data and unsupported version must leave output intact.
		decoded.tick = 999;
		for (size_t size = 0; size < bytes.size(); ++size)
		{
			std::vector<uint8> bad(bytes.begin(), bytes.begin() + ptrdiff_t(size));
			CHECK_FALSE(DecodeCheckpoint(bad, Limits {}, decoded));
			CHECK(decoded.tick == 999);
		}
		auto bad = bytes; bad.push_back(0);
		CHECK_FALSE(DecodeCheckpoint(bad, Limits {}, decoded));
		bad = bytes; bad[8] = 99;
		CHECK_FALSE(DecodeCheckpoint(bad, Limits {}, decoded));
		bad = bytes; bad[0] ^= 1;
		CHECK_FALSE(DecodeCheckpoint(bad, Limits {}, decoded));
		Limits small; small.entities = 0;
		CHECK_FALSE(DecodeCheckpoint(bytes, small, decoded));
		small = Limits {}; small.allocation_bytes = 1;
		CHECK_FALSE(DecodeCheckpoint(bytes, small, decoded));
		small = Limits {}; small.bytes = bytes.size() - 1;
		CHECK_FALSE(DecodeCheckpoint(bytes, small, decoded));
		CHECK_FALSE(EncodeCheckpoint(state, small, again));
		CHECK(again == bytes);
		state.entities[0].mass = -1;
		CHECK_FALSE(EncodeCheckpoint(state, Limits {}, again));
		CHECK(decoded.tick == 999);
	}

	TEST_CASE("Impact ordering uses tick fraction source target component then stable ID")
	{
		ImpactEvent a; a.tick = 1; a.fraction = 2; a.source = 2; a.target = 7; a.component = 1; a.id = 90;
		ImpactEvent b = a; b.fraction = 1; b.id = 99;
		ImpactEvent c = a; c.source = 1; c.id = 98;
		ImpactEvent d = a; d.target = 6; d.id = 97;
		ImpactEvent e = a; e.component = 0; e.id = 96;
		ImpactEvent f = a; f.id = 1;
		ImpactEvent g = a; g.tick = 2; g.fraction = 0; g.id = 100;
		std::vector<ImpactEvent> input { a, b, c, d, e, f, g };
		for (size_t permutation = 0; permutation < input.size(); ++permutation)
		{
			std::rotate(input.begin(), input.begin() + 1, input.end());
			std::vector<ID> order;
			ImpactBatch batch(input, Limits {});
			CHECK(batch.Process([](const ImpactEvent &) { CHECK(false); return false; }, [&](const ImpactEvent &event) { order.push_back(event.id); return true; }));
			CHECK(order == std::vector<ID> { 99, 98, 97, 96, 1, 90, 100 });
		}
	}

	TEST_CASE("Impact batches order deterministically and have one response owner")
	{
		ImpactEvent query;
		query.tick = 1; query.id = 9; query.source = 2; query.target = 7;
		query.owner = ResponseOwner::Query;
		query.impulse[3] = 2;
		ImpactEvent solver = query;
		solver.id = 8; solver.owner = ResponseOwner::Solver;
		std::vector<ID> observed;
		float impulse = 0;
		ImpactBatch batch({ query, solver }, Limits {});
		CHECK(batch.Process([&](const ImpactEvent &e) { impulse += e.impulse[3]; return true; },
			[&](const ImpactEvent &e) { observed.push_back(e.id); return true; }));
		CHECK(observed == std::vector<ID> { 8, 9 });
		CHECK(impulse == 2);
		CHECK_FALSE(batch.Process([](const ImpactEvent &) { return true; }, [](const ImpactEvent &) { return true; }));
		ImpactBatch duplicate({ query, query }, Limits {});
		CHECK_FALSE(duplicate.Process([&](const ImpactEvent &) { ++impulse; return true; }, [](const ImpactEvent &) { return true; }));
		CHECK(impulse == 2);
		Limits small; small.events = 1;
		ImpactBatch overflow({ query, solver }, small);
		CHECK_FALSE(overflow.Process([](const ImpactEvent &) { return true; }, [](const ImpactEvent &) { return true; }));
		query.energy = -1;
		ImpactBatch invalid({ query }, Limits {});
		CHECK_FALSE(invalid.Process([](const ImpactEvent &) { return true; }, [](const ImpactEvent &) { return true; }));
	}

	TEST_CASE("Commands validate authority timing sequence bounds and all ten axes")
	{
		Checkpoint state;
		EntityState ship;
		ship.id = 7;
		ship.owner = 3;
		state.entities.push_back(ship);
		Limits limits;
		limits.commands = 2;
		limits.future_ticks = 4;
		Command cmd;
		cmd.client = 3;
		cmd.entity = 7;
		cmd.sequence = 1;
		cmd.tick = 1;
		CHECK(SubmitCommand(state, limits, cmd) == CommandResult::Accepted);
		CHECK(SubmitCommand(state, limits, cmd) == CommandResult::Sequence);
		cmd.sequence = 2;
		cmd.client = 8;
		CHECK(SubmitCommand(state, limits, cmd) == CommandResult::Authority);
		cmd.client = 3;
		cmd.tick = 0;
		CHECK(SubmitCommand(state, limits, cmd) == CommandResult::Late);
		cmd.tick = 5;
		CHECK(SubmitCommand(state, limits, cmd) == CommandResult::Future);
		cmd.tick = 2;
		for (size_t axis = 0; axis < 10; ++axis)
		{
			Command invalid = cmd;
			if (axis < 4) invalid.translation[axis] = std::numeric_limits<float>::infinity();
			else invalid.rotation[axis - 4] = std::numeric_limits<float>::quiet_NaN();
			CHECK(SubmitCommand(state, limits, invalid) == CommandResult::Invalid);
		}
		cmd.rotation[5] = 1.01f;
		CHECK(SubmitCommand(state, limits, cmd) == CommandResult::Invalid);
		cmd.rotation[5] = -1.0f;
		CHECK(SubmitCommand(state, limits, cmd) == CommandResult::Accepted);
		cmd.sequence = 3;
		CHECK(SubmitCommand(state, limits, cmd) == CommandResult::Full);
		CHECK(state.streams[0].accepted == 2);
		CHECK(state.streams[0].applied == 0);
		CHECK(state.commands.size() == 2);
	}
}
