// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#include <Jolt/SpaceSim/World.h>
#include <limits>
#include <random>
using namespace JPH::SpaceSim;

TEST_SUITE("SpaceSimWorldTests")
{
	TEST_CASE("SectorArithmeticPreservesMillimetersAtInt64Extremes")
	{
		GlobalPosition p;
		p.mSector = {INT64_MAX - 1, INT64_MIN + 1, INT64_MAX - 1, INT64_MIN + 1};
		p.mOffset = {0.125, 0.25, 0.5, 0.75};
		GlobalPosition q;
		CHECK(p.Translate({0.001, -0.002, 0.003, -0.004}, q));
		Vector d;
		CHECK(q.Difference(p, 1.0, d));
		for (size_t i = 0; i < 4; ++i)
			CHECK(d[i] == doctest::Approx(Vector{0.001, -0.002, 0.003, -0.004}[i]).epsilon(1.0e-12));
		GlobalPosition edge;
		edge.mSector[0] = INT64_MAX;
		edge.mOffset[0] = GlobalPosition::cSectorSize - 0.5;
		q = p;
		CHECK_FALSE(edge.Translate({1, 0, 0, 0}, q));
		CHECK(q.mSector == p.mSector);
		CHECK_FALSE(p.Translate({std::numeric_limits<double>::infinity(), 0, 0, 0}, q));
		edge.mSector[0] = INT64_MIN;
		edge.mOffset[0] = 0;
		CHECK_FALSE(edge.Translate({-1, 0, 0, 0}, q));
		CHECK_FALSE(edge.Difference(p, 65536, d));
	}

	TEST_CASE("FrameSetTransferPreservesTinyMotionAndMetadataAtHugeCommonVelocity")
	{
		World world(0.5);
		Frame a;
		a.mID = 1;
		a.mOrigin.mSector = {INT64_MAX - (INT64_C(1) << 50), INT64_MIN + (INT64_C(1) << 50), INT64_MAX - (INT64_C(1) << 50), INT64_MIN + (INT64_C(1) << 50)};
		a.mVelocity = {1.0e18, -1.0e18, 1.0e18, -1.0e18};
		Frame b = a;
		b.mID = 2;
		CHECK(a.mOrigin.Translate({10, -20, 30, -40}, b.mOrigin));
		CHECK(world.AddFrame(a));
		CHECK(world.AddFrame(b));
		ObjectState first;
		first.mID = 100;
		first.mFrame = 1;
		first.mPosition = {1, 2, 3, 4};
		first.mVelocity = {0.001, -0.002, 0.003, -0.004};
		first.mHealth = 0.73;
		first.mTopologyVersion = UINT64_MAX;
		first.mPopulationSeed = 987654321;
		first.mOrientation = Rotor::sExp(Bivec(0, 0, 0.2f, 0, 0, 0));
		first.mAngularVelocity = Bivec(0.1f, 0, 0, 0, 0, 0.2f);
		ObjectState second = first;
		second.mID = 101;
		second.mPosition[3] += 0.0005;
		CHECK(world.AddObject(first));
		CHECK(world.AddObject(second));
		CHECK(world.TransformInteractingSet({101, 100}, 2, 20));
		auto result = world.FindObject(100);
		if (!result) { CHECK(false); return; }
		CHECK(result->mFrame == 2);
		CHECK(result->mEpoch == 20);
		for (size_t i = 0; i < 4; ++i)
		{
			CHECK(result->mVelocity[i] == first.mVelocity[i]);
			CHECK(result->mPosition[i] == doctest::Approx(first.mPosition[i] + first.mVelocity[i] * 10 - Vector{10, -20, 30, -40}[i]).epsilon(1.0e-12));
		}
		CHECK(result->mHealth == first.mHealth);
		CHECK(result->mTopologyVersion == first.mTopologyVersion);
		CHECK(result->mPopulationSeed == first.mPopulationSeed);
		CHECK(result->mOrientation == first.mOrientation);
		CHECK(result->mAngularVelocity == first.mAngularVelocity);
		CHECK(world.FindObject(101)->mPosition[3] - result->mPosition[3] == doctest::Approx(0.0005).epsilon(1.0e-10));
		CHECK(world.GetRepresentationVersion() == 1);
		CHECK(world.TransformInteractingSet({100, 101}, 1, 20));
		CHECK(world.FindObject(100)->mPosition[3] == doctest::Approx(3.96).epsilon(1.0e-12));
	}
	TEST_CASE("FrameTransferCannotDiscardTinyVelocityThatCreatesALongHorizonEncounter")
	{
		for (size_t axis = 0; axis < 4; ++axis)
			for (double boost : {1.0e18, -1.0e18})
			{
				World world(1);
				Frame stationary; stationary.mID = 1;
				Frame moving = stationary; moving.mID = 2; moving.mVelocity[axis] = boost;
				CHECK(world.AddFrame(stationary)); CHECK(world.AddFrame(moving));
				ObjectState approaching; approaching.mID = 2; approaching.mFrame = 1;
				approaching.mVelocity[axis] = 1.0e-10;
				approaching.mHasSourceEpoch = true; approaching.mSourceEpoch = 7;
				ObjectState target = approaching; target.mID = 1;
				target.mPosition[axis] = 0.1; target.mVelocity[axis] = 0;
				CHECK(world.AddObject(target)); CHECK(world.AddObject(approaching));
				ActivationResult before, after;
				CHECK(world.FindInteractions(0, 1000000000, before));
				CHECK(before.mPairs == std::vector<InteractionPair>{{1, 2}});
				GlobalPosition expected;
				CHECK(world.EvaluateGlobal(2, 1000000000, expected));
				// Stage the exactly representable target first, then reject the
				// approaching body's lost speed. Neither body may be committed.
				CHECK_FALSE(world.TransformInteractingSet({1, 2}, 2, 0));
				CHECK(world.FindObject(1)->mFrame == 1);
				CHECK(world.FindObject(2)->mFrame == 1);
				CHECK(world.FindObject(2)->mVelocity == approaching.mVelocity);
				CHECK(world.FindObject(2)->mEpoch == 0);
				CHECK(world.FindObject(2)->GetSourceEpoch() == 7);
				CHECK(world.GetRepresentationVersion() == 0);
				CHECK(world.TransformInteractingSet({1, 2}, 1, 0));
				CHECK(world.FindInteractions(0, 1000000000, after));
				CHECK(after.mPairs == before.mPairs);
				GlobalPosition actual;
				Vector difference;
				CHECK(world.EvaluateGlobal(2, 1000000000, actual));
				CHECK(actual.Difference(expected, 1, difference));
				CHECK(difference == Vector{});
			}
	}

	TEST_CASE("RepresentableBoostRoundTripPreservesLongHorizonCandidatesAndGlobalDrift")
	{
		for (double boost : {0x1p20, -0x1p20})
		{
			World world(1);
			Frame source; source.mID = 1;
			Frame destination = source; destination.mID = 2; destination.mVelocity[3] = boost;
			CHECK(world.AddFrame(source)); CHECK(world.AddFrame(destination));
			ObjectState a; a.mID = 1; a.mFrame = 1; a.mVelocity[3] = 0x1p-30;
			ObjectState b = a; b.mID = 2; b.mVelocity[3] = 0; b.mPosition[3] = 0.5;
			CHECK(world.AddObject(a)); CHECK(world.AddObject(b));
			ActivationResult before, during, after;
			CHECK(world.FindInteractions(0, 1000000000, before));
			CHECK(before.mPairs == std::vector<InteractionPair>{{1, 2}});
			GlobalPosition expected, actual;
			CHECK(world.EvaluateGlobal(1, 1000000000, expected));
			CHECK(world.TransformInteractingSet({1, 2}, 2, 0));
			CHECK(world.FindInteractions(0, 1000000000, during));
			CHECK(during.mPairs == before.mPairs);
			CHECK(world.EvaluateGlobal(1, 1000000000, actual));
			Vector difference;
			CHECK(actual.Difference(expected, 1, difference));
			CHECK(std::abs(difference[3]) <= cTransferError);
			CHECK(world.TransformInteractingSet({1, 2}, 1, 0));
			CHECK(world.FindObject(1)->mVelocity == a.mVelocity);
			CHECK(world.FindObject(1)->GetSourceEpoch() == a.GetSourceEpoch());
			CHECK(world.FindInteractions(0, 1000000000, after));
			CHECK(after.mPairs == before.mPairs);
			CHECK(world.EvaluateGlobal(1, 1000000000, actual));
			CHECK(actual.mSector == expected.mSector);
			CHECK(actual.mOffset == expected.mOffset);
		}
	}

	TEST_CASE("PreparationPropagatesUnrepresentableVelocityWithoutPromotingOrPublishing")
	{
		World world(1);
		Frame source; source.mID = 1;
		Frame destination = source; destination.mID = 2; destination.mVelocity[0] = 1.0e18;
		CHECK(world.AddFrame(source)); CHECK(world.AddFrame(destination));
		ObjectState a; a.mID = 1; a.mFrame = 2; a.mTier = Tier::Distant;
		ObjectState b = a; b.mID = 2; b.mFrame = 1; b.mVelocity[0] = 1.0e-10;
		CHECK(world.AddObject(a)); CHECK(world.AddObject(b));
		ActivationResult plan;
		CHECK(world.FindInteractions(0, 0, plan));
		CHECK(plan.mPairs == std::vector<InteractionPair>{{1, 2}});
		CHECK(plan.mGroups.size() == 1);
		if (plan.mGroups.empty()) return;
		CHECK(plan.mGroups[0].mFrame == 2);
		ActivationResult result; result.mTick = 123; result.mPairs = {{7, 9}};
		CHECK_FALSE(world.PrepareInteractions(0, 0, result));
		CHECK(result.mTick == 123);
		CHECK(result.mPairs == std::vector<InteractionPair>{{7, 9}});
		CHECK(world.FindObject(1)->mTier == Tier::Distant);
		CHECK(world.FindObject(2)->mTier == Tier::Distant);
		CHECK(world.FindObject(2)->mFrame == 1);
		CHECK(world.FindObject(2)->mVelocity == b.mVelocity);
		CHECK(world.GetRepresentationVersion() == 0);
	}

	TEST_CASE("VelocityTransferCanCancelIntermediateSubtractionResidualExactly")
	{
		World world(1);
		Frame source; source.mID = 1; source.mVelocity[0] = 1.0e-10;
		Frame destination; destination.mID = 2; destination.mVelocity[0] = 1.0e18;
		CHECK(world.AddFrame(source)); CHECK(world.AddFrame(destination));
		ObjectState object; object.mID = 1; object.mFrame = 1; object.mVelocity[0] = -source.mVelocity[0];
		CHECK(world.AddObject(object));
		CHECK(world.TransformInteractingSet({1}, 2, 0));
		CHECK(world.FindObject(1)->mVelocity[0] == -destination.mVelocity[0]);
		CHECK(world.TransformInteractingSet({1}, 1, 0));
		CHECK(world.FindObject(1)->mVelocity == object.mVelocity);
	}

	TEST_CASE("CrossFrameProjectilePromotesWholeInteractionBeforeSolving")
	{
		for (double common : {0.0, 1.0e18})
		{
			World world(1.0);
			Frame ship_frame;
			ship_frame.mID = 1;
			if (common != 0)
				ship_frame.mOrigin.mSector = {INT64_C(1) << 60, -(INT64_C(1) << 60), 0, INT64_C(1) << 60};
			ship_frame.mVelocity[0] = common;
			Frame projectile_frame = ship_frame;
			projectile_frame.mID = 2;
			CHECK(ship_frame.mOrigin.Translate({0, 0, 0, 100}, projectile_frame.mOrigin));
			CHECK(world.AddFrame(ship_frame));
			CHECK(world.AddFrame(projectile_frame));
			ObjectState ship;
			ship.mID = 99;
			ship.mFrame = 1;
			ship.mTier = Tier::Distant;
			ship.mRotationRadius = 2;
			ship.mHealth = 0.33;
			ship.mTopologyVersion = 81;
			ObjectState projectile = ship;
			projectile.mID = 7;
			projectile.mFrame = 2;
			projectile.mVelocity[3] = -200;
			projectile.mRotationRadius = 0.01;
			CHECK(world.AddObject(ship));
			CHECK(world.AddObject(projectile));
			ActivationResult result;
			CHECK(world.PrepareInteractions(0, 1, result));
			CHECK(result.mPairs == std::vector<InteractionPair>{{7, 99}});
			CHECK(result.mPromoted == std::vector<ObjectID>{7, 99});
			CHECK(result.mGroups.size() == 1);
			CHECK(world.FindObject(7)->mFrame == world.FindObject(99)->mFrame);
			CHECK(world.FindObject(99)->mTier == Tier::Detailed);
			CHECK(world.FindObject(99)->mHealth == 0.33);
			CHECK(world.FindObject(99)->mTopologyVersion == 81);
		}
	}

	TEST_CASE("AccelerationRotationWeaponReachAndSchedulingAllExpandBounds")
	{
		for (int mode = 0; mode < 4; ++mode)
		{
			World world(1);
			Frame frame;
			frame.mID = 1;
			CHECK(world.AddFrame(frame));
			ObjectState a;
			a.mID = 1;
			a.mFrame = 1;
			ObjectState b = a;
			b.mID = 2;
			b.mPosition[3] = 10;
			if (mode == 0) a.mAccelerationBound = 20;
			if (mode == 1) a.mRotationRadius = 10;
			if (mode == 2) a.mWeaponReach = 10;
			if (mode == 3) { a.mVelocity[3] = 1; a.mSchedulingInterval = 10; }
			CHECK(world.AddObject(a));
			CHECK(world.AddObject(b));
			ActivationResult result;
			CHECK(world.PrepareInteractions(0, 1, result));
			CHECK(result.mPairs == std::vector<InteractionPair>{{1, 2}});
		}
	}

	TEST_CASE("SparseFourDimensionalIndexIsDeterministicWithoutAllPairs")
	{
		World first(1), reversed(1);
		Frame frame;
		frame.mID = 1;
		frame.mVelocity[0] = 1.0e18;
		CHECK(first.AddFrame(frame));
		CHECK(reversed.AddFrame(frame));
		for (uint64_t i = 0; i < 2048; ++i)
		{
			ObjectState object;
			object.mID = i + 1;
			object.mFrame = 1;
			object.mPosition[3] = double(i) * 16;
			object.mRotationRadius = 1;
			CHECK(first.AddObject(object));
			object.mID = 2048 - i;
			object.mPosition[3] = double(2047 - i) * 16;
			CHECK(reversed.AddObject(object));
		}
		ActivationResult a, b;
		CHECK(first.FindInteractions(3, 1, a));
		CHECK(reversed.FindInteractions(3, 1, b));
		CHECK(a.mPairs.empty());
		CHECK(a.mPairs == b.mPairs);
		CHECK(a.mNodeVisits < 2048 * 64);
		CHECK(a.mNodeVisits == b.mNodeVisits);
		CHECK(a.mResidentCount == 2048);
	}

	TEST_CASE("FramePropagationRetainsProductResidualRatherThanRoundingHugeDisplacement")
	{
		Frame frame;
		frame.mID = 1;
		frame.mVelocity = {1.0e18 + 128, 0, 0, 0};
		GlobalPosition expected, actual;
		CHECK(frame.mOrigin.Translate(frame.mVelocity, expected));
		CHECK(expected.Translate({0.5 * frame.mVelocity[0], 0, 0, 0}, expected));
		CHECK(frame.OriginAt(3, 0.5, actual));
		Vector delta;
		CHECK(actual.Difference(expected, 100, delta));
		CHECK(delta[0] == 0);
	}

	TEST_CASE("ResamplingCannotResetAccelerationEnvelopeOrAttitudeEpoch")
	{
		Frame frame;
		frame.mID = 1;
		ObjectState object;
		object.mID = 1; object.mFrame = 1;
		object.mAccelerationBound = 1;
		ObjectState advanced;
		CHECK(object.Evaluate(20, 1, advanced));
		GlobalBounds original, resampled;
		CHECK(BuildSweptBounds(object, frame, frame, 20, 0, 1, original));
		CHECK(BuildSweptBounds(advanced, frame, frame, 20, 0, 1, resampled));
		GlobalPosition extreme;
		CHECK(frame.mOrigin.Translate({200, 0, 0, 0}, extreme));
		CHECK(original.Contains(extreme));
		CHECK(resampled.Contains(extreme));
	}

	TEST_CASE("FloatPhysicsBoundaryRequiresExplicitRangeAndRoundingBudget")
	{
		ObjectState object;
		object.mID = 1;
		object.mFrame = 1;
		object.mPosition = {0.001, -2, 3, -4};
		object.mVelocity = {0.0001, 0, 0, 0};
		PhysicsState result;
		PhysicsLimits limits;
		CHECK(TryPhysicsState(object, limits, result));
		CHECK(result.mPosition.GetW() == -4);
		CHECK(result.mVelocity.GetX() == doctest::Approx(0.0001));
		object.mPosition[0] = 65536;
		CHECK_FALSE(TryPhysicsState(object, limits, result));
		limits.mMaxPosition = 65536;
		object.mPosition[0] = 65535.0001;
		limits.mPositionError = 1.0e-8;
		CHECK_FALSE(TryPhysicsState(object, limits, result));
		object.mPosition[0] = 0;
		object.mVelocity[3] = 1.0e18 + 0.001;
		CHECK_FALSE(TryPhysicsState(object, limits, result));
		limits.mVelocityError = std::numeric_limits<double>::infinity();
		CHECK_FALSE(TryPhysicsState(object, limits, result));
	}

	TEST_CASE("DemotionUsesHysteresisBlockersAndFutureInteractions")
	{
		for (uint32_t blocker : {0u, uint32_t(ActiveContact), uint32_t(PendingImpulse), uint32_t(WeaponInteraction), uint32_t(TopologyTransaction)})
		{
			World world(1);
			Frame frame;
			frame.mID = 1;
			CHECK(world.AddFrame(frame));
			ObjectState object;
			object.mID = 1;
			object.mFrame = 1;
			object.mTier = Tier::Detailed;
			object.mLastPromotion = 10;
			object.mDemotionBlockers = blocker;
			object.mVelocity[3] = 0.001;
			object.mHealth = 0.42;
			object.mTopologyVersion = 17;
			CHECK(world.AddObject(object));
			CHECK_FALSE(world.TryDemote(1, Tier::Distant, 14, 5, 1));
			CHECK(world.TryDemote(1, Tier::Distant, 15, 5, 1) == (blocker == 0));
			if (blocker == 0)
			{
				CHECK(world.FindObject(1)->mPosition[3] == doctest::Approx(0.015));
				CHECK(world.FindObject(1)->mHealth == 0.42);
				CHECK(world.FindObject(1)->mTopologyVersion == 17);
				ObjectState predicted;
				CHECK(world.FindObject(1)->Evaluate(1000, 1, predicted));
				CHECK(predicted.mPosition[3] == doctest::Approx(1));
			}
		}
		World world(1);
		Frame frame;
		frame.mID = 1;
		CHECK(world.AddFrame(frame));
		ObjectState a;
		a.mID = 1; a.mFrame = 1; a.mTier = Tier::Detailed;
		ObjectState b = a;
		b.mID = 2; b.mPosition[3] = 100; b.mVelocity[3] = -100;
		CHECK(world.AddObject(a)); CHECK(world.AddObject(b));
		CHECK_FALSE(world.TryDemote(1, Tier::Trajectory, 0, 0, 1));
	}

	TEST_CASE("DetailedOwnerCanPublishStateWithoutBypassingFrameAndTierOwnership")
	{
		World world(1);
		Frame frame; frame.mID = 1;
		CHECK(world.AddFrame(frame));
		ObjectState object; object.mID = 1; object.mFrame = 1; object.mTier = Tier::Detailed;
		object.mDemotionBlockers = ActiveContact;
		CHECK(world.AddObject(object));
		ObjectState updated = object;
		updated.mEpoch = 10; updated.mHealth = 0.5; updated.mTopologyVersion = 4;
		updated.mPosition[3] = 7; updated.mDemotionBlockers = 0;
		CHECK(world.CommitDetailedState(updated));
		CHECK(world.FindObject(1)->mHealth == 0.5);
		CHECK(world.FindObject(1)->GetSourceEpoch() == 10);
		CHECK_FALSE(world.CommitDetailedState(object));
		updated.mFrame = 2;
		CHECK_FALSE(world.CommitDetailedState(updated));
		updated.mFrame = 1; updated.mTier = Tier::Distant;
		CHECK_FALSE(world.CommitDetailedState(updated));
		CHECK(world.TryDemote(1, Tier::Distant, 10, 5, 1));
		CHECK(world.FindObject(1)->mPosition[3] == 7);

	}

	TEST_CASE("DistantFreeDriftCrossesSectorsWithoutFreezingOrFlattening")
	{
		World world(0.5);
		Frame frame; frame.mID = 1;
		frame.mOrigin.mSector = {INT64_C(1) << 60, -(INT64_C(1) << 60), INT64_C(1) << 60, -(INT64_C(1) << 60)};
		CHECK(world.AddFrame(frame));
		ObjectState object; object.mID = 1; object.mFrame = 1; object.mTier = Tier::Distant;
		object.mPosition = {0.125, 0.25, 0.5, 0.75};
		object.mVelocity = {1048576, -1048576, 2097152, -2097152};
		CHECK(world.AddObject(object));
		GlobalPosition result;
		CHECK(world.EvaluateGlobal(1, 2000000, result));
		CHECK(result.mSector == std::array<int64_t, 4>{frame.mOrigin.mSector[0] + 1000000, frame.mOrigin.mSector[1] - 1000000,
			frame.mOrigin.mSector[2] + 2000000, frame.mOrigin.mSector[3] - 2000000});
		CHECK(result.mOffset == object.mPosition);
		CHECK(world.FindObject(1)->mEpoch == 0); // On-demand trajectory queries do not accumulate step error.

	}

	TEST_CASE("IndexContainsEveryAnalyticalFourDimensionalSphereEncounter")
	{
		World world(0.25);
		Frame frame; frame.mID = 1;
		frame.mOrigin.mSector = {INT64_C(1) << 60, -(INT64_C(1) << 60), INT64_C(1) << 60, -(INT64_C(1) << 60)};
		frame.mVelocity = {1.0e18, -1.0e18, 1.0e18, -1.0e18};
		CHECK(world.AddFrame(frame));
		std::mt19937 random(0x714be3u);
		std::vector<ObjectState> objects;
		for (ObjectID id = 1; id <= 192; ++id)
		{
			ObjectState object; object.mID = id; object.mFrame = 1;
			for (size_t axis = 0; axis < 4; ++axis)
			{
				object.mPosition[axis] = double(int(random() % 41) - 20);
				object.mVelocity[axis] = double(int(random() % 21) - 10);
			}
			object.mRotationRadius = 2;
			objects.push_back(object);
			CHECK(world.AddObject(object));
		}
		ActivationResult result;
		CHECK(world.FindInteractions(0, 8, result));
		CHECK(std::is_sorted(result.mPairs.begin(), result.mPairs.end()));
		CHECK(std::adjacent_find(result.mPairs.begin(), result.mPairs.end()) == result.mPairs.end());
		uint64_t encounters = 0;
		for (size_t a = 0; a < objects.size(); ++a)
			for (size_t b = a + 1; b < objects.size(); ++b)
			{
				double vv = 0, pv = 0;
				Vector p, v;
				for (size_t axis = 0; axis < 4; ++axis)
				{
					p[axis] = objects[a].mPosition[axis] - objects[b].mPosition[axis];
					v[axis] = objects[a].mVelocity[axis] - objects[b].mVelocity[axis];
					vv += v[axis] * v[axis]; pv += p[axis] * v[axis];
				}
				double t = vv == 0? 0 : std::clamp(-pv / vv, 0.0, 2.0);
				double separation_sq = 0;
				for (size_t axis = 0; axis < 4; ++axis) { double d = p[axis] + t * v[axis]; separation_sq += d * d; }
				if (separation_sq <= 16)
				{
					++encounters;
					CHECK(std::binary_search(result.mPairs.begin(), result.mPairs.end(), InteractionPair{objects[a].mID, objects[b].mID}));
				}
			}
		CHECK(encounters > 10);
	}

	TEST_CASE("TransitiveInteractionComponentsShareOneFrameAndFailuresAreAtomic")
	{
		World world(1);
		for (FrameID id : {1u, 2u, 3u})
		{
			Frame frame; frame.mID = id;
			CHECK(world.AddFrame(frame));
			ObjectState object; object.mID = id; object.mFrame = id;
			object.mPosition[3] = double(id) * 3;
			object.mRotationRadius = 2;
			CHECK(world.AddObject(object));
		}
		ActivationResult result;
		CHECK(world.PrepareInteractions(0, 1, result));
		CHECK(result.mPairs == std::vector<InteractionPair>{{1, 2}, {2, 3}});
		CHECK(result.mGroups.size() == 1);
		CHECK(result.mGroups[0].mObjects == std::vector<ObjectID>{1, 2, 3});
		for (ObjectID id : {1u, 2u, 3u})
		{
			CHECK(world.FindObject(id)->mFrame == 1);
			CHECK(world.FindObject(id)->mTier == Tier::Detailed);
		}
		World impossible(1);
		Frame a; a.mID = 1;
		Frame b = a; b.mID = 2; b.mOrigin.mSector[3] = 1;
		CHECK(impossible.AddFrame(a)); CHECK(impossible.AddFrame(b));
		ObjectState one; one.mID = 1; one.mFrame = 1; one.mWeaponReach = GlobalPosition::cSectorSize;
		ObjectState two = one; two.mID = 2; two.mFrame = 2;
		CHECK(impossible.AddObject(one)); CHECK(impossible.AddObject(two));
		result.mTick = 987;
		CHECK_FALSE(impossible.PrepareInteractions(0, 1, result));
		CHECK(result.mTick == 987);
		CHECK(impossible.FindObject(1)->mTier == Tier::Trajectory);
		CHECK(impossible.FindObject(2)->mFrame == 2);
		CHECK(impossible.GetRepresentationVersion() == 0);
	}

	TEST_CASE("SectorSeparatedPopulationCannotAliasAtSignedExtremes")
	{
		World world(1);
		for (FrameID id : {1u, 2u})
		{
			Frame frame; frame.mID = id;
			frame.mOrigin.mSector[3] = id == 1? INT64_MIN + 1 : INT64_MAX - 1;
			CHECK(world.AddFrame(frame));
			ObjectState object; object.mID = id; object.mFrame = id; object.mRotationRadius = 1;
			CHECK(world.AddObject(object));
		}
		ActivationResult result;
		CHECK(world.FindInteractions(0, 1, result));
		CHECK(result.mPairs.empty());
	}

	TEST_CASE("DifferentFrameEpochsAndVelocityDeltasRequireRepresentableLocalMotion")
	{
		for (double local_velocity : {0.001, 0x1p-10})
		{
			World world(0.5);
			Frame a; a.mID = 1; a.mVelocity[0] = 1.0e18 + 128;
			Frame b = a; b.mID = 2; b.mEpoch = 3;
			CHECK(a.OriginAt(3, 0.5, b.mOrigin));
			b.mVelocity[0] += 128;
			CHECK(world.AddFrame(a)); CHECK(world.AddFrame(b));
			ObjectState object; object.mID = 1; object.mFrame = 1;
			object.mVelocity[0] = local_velocity; object.mPosition[3] = 0.00025;
			CHECK(world.AddObject(object));
			if (local_velocity == 0.001)
			{
				// -128 + binary64(0.001) is not exactly representable. Its tiny
				// speed residual cannot be charged to a fixed position budget.
				CHECK_FALSE(world.TransformInteractingSet({1}, 2, 5));
				CHECK(world.FindObject(1)->mFrame == 1);
				CHECK(world.FindObject(1)->mEpoch == 0);
				CHECK(world.FindObject(1)->mVelocity == object.mVelocity);
				CHECK(world.GetRepresentationVersion() == 0);
				continue;
			}
			CHECK(world.TransformInteractingSet({1}, 2, 5));
			const auto *result = world.FindObject(1);
			CHECK(result->mPosition[0] == std::fma(local_velocity, 2.5, -128.0));
			CHECK(result->mVelocity[0] == -128 + local_velocity);
			CHECK(result->mPosition[3] == 0.00025);
			CHECK(result->GetSourceEpoch() == 0);
		}
	}

	TEST_CASE("FrameVelocitySubtractionResidualCannotAmplifyIntoLostPosition")
	{
		Frame source; source.mID = 1; source.mVelocity[0] = 1.0e-10;
		Frame reference; reference.mID = 2; reference.mVelocity[0] = 1.0e18;
		GlobalPosition expected, actual;
		CHECK(TranslateProduct(source.mOrigin, {-1.0e18, 0, 0, 0}, 1000000, expected));
		CHECK(expected.Translate({0.0001, 0, 0, 0}, expected));
		CHECK(source.IndexOrigin(reference, 1000000, 1, actual));
		Vector difference;
		CHECK(actual.Difference(expected, 1, difference));
		CHECK(std::abs(difference[0]) < 1.0e-9);
	}

	TEST_CASE("GlobalOverflowOfLocalBodyCannotHideBehindValidFrameOrigin")
	{
		World world(1);
		Frame frame; frame.mID = 1;
		frame.mOrigin.mSector[0] = INT64_MAX;
		frame.mOrigin.mOffset[0] = GlobalPosition::cSectorSize - 10;
		CHECK(world.AddFrame(frame));
		ObjectState object; object.mID = 1; object.mFrame = 1;
		object.mPosition[0] = 100;
		CHECK_FALSE(world.AddObject(object));
		object.mID = 2; object.mPosition[0] = 0; object.mVelocity[0] = 100;
		CHECK(world.AddObject(object));
		CHECK_FALSE(world.TransformInteractingSet({2}, 1, 1));
		CHECK(world.FindObject(2)->mEpoch == 0);
		CHECK(world.GetRepresentationVersion() == 0);
	}

	TEST_CASE("SubnormalAccelerationCannotDisappearBeforeHugeTimeSquaring")
	{
		Frame frame;
		frame.mID = 1;
		ObjectState object;
		object.mID = 1; object.mFrame = 1;
		object.mAccelerationBound = std::numeric_limits<double>::denorm_min();
		GlobalBounds bounds;
		CHECK_FALSE(BuildSweptBounds(object, frame, frame, 0, 1, 1.0e180, bounds));
	}

	TEST_CASE("AtomicFrameFailurePreservesEveryObjectAndRevision")
	{
		World world(1);
		Frame a; a.mID = 1; a.mVelocity[0] = 1.0e18;
		Frame b = a; b.mID = 2;
		Frame far = a; far.mID = 3; far.mOrigin.mSector[3] = INT64_MAX;
		Frame stationary; stationary.mID = 4;
		CHECK(world.AddFrame(a)); CHECK(world.AddFrame(b));
		CHECK(world.AddFrame(far)); CHECK(world.AddFrame(stationary));
		ObjectState first; first.mID = 1; first.mFrame = 1;
		first.mVelocity[0] = 0.001;
		ObjectState second = first; second.mID = 2; second.mFrame = 3;
		CHECK(world.AddObject(first)); CHECK(world.AddObject(second));
		CHECK_FALSE(world.TransformInteractingSet({1, 2}, 2, 1));
		CHECK(world.FindObject(1)->mEpoch == 0);
		CHECK(world.FindObject(1)->mFrame == 1);
		CHECK(world.FindObject(1)->mPosition == first.mPosition);
		CHECK(world.FindObject(2)->mFrame == 3);
		CHECK(world.GetRepresentationVersion() == 0);
		CHECK_FALSE(world.TransformInteractingSet({1}, 4, 0)); // Would lose 0.001 in 1e18.
		CHECK_FALSE(world.TransformInteractingSet({1, 1}, 2, 0));
		CHECK_FALSE(world.TransformInteractingSet({1, 99}, 2, 0));
		CHECK_FALSE(world.TransformInteractingSet({}, 2, 0));
		CHECK(world.GetRepresentationVersion() == 0);
	}

	TEST_CASE("ValidationRejectsNonfiniteOverflowAndMalformedIdentityWithoutPartialOutput")
	{
		World invalid(std::numeric_limits<double>::quiet_NaN());
		Frame frame; frame.mID = 1;
		CHECK_FALSE(invalid.AddFrame(frame));
		World world(1);
		CHECK(world.AddFrame(frame));
		CHECK_FALSE(world.AddFrame(frame));
		Frame bad_frame = frame;
		bad_frame.mID = 2; bad_frame.mVelocity[3] = std::numeric_limits<double>::infinity();
		CHECK_FALSE(world.AddFrame(bad_frame));
		ObjectState object; object.mID = UINT64_MAX; object.mFrame = 1;
		CHECK(world.AddObject(object));
		CHECK_FALSE(world.AddObject(object));
		ObjectState bad = object;
		bad.mID = 2; bad.mPosition[3] = std::numeric_limits<double>::quiet_NaN();
		CHECK_FALSE(world.AddObject(bad));
		bad = object; bad.mID = 2; bad.mVelocity[3] = std::numeric_limits<double>::infinity();
		CHECK_FALSE(world.AddObject(bad));
		bad = object; bad.mID = 2; bad.mAccelerationBound = -1;
		CHECK_FALSE(world.AddObject(bad));
		bad = object; bad.mID = 2; bad.mRotationRadius = std::numeric_limits<double>::infinity();
		CHECK_FALSE(world.AddObject(bad));
		bad = object; bad.mID = 2; bad.mFrame = 33;
		CHECK_FALSE(world.AddObject(bad));
		ActivationResult output; output.mTick = 123;
		CHECK_FALSE(world.FindInteractions(UINT64_MAX, 1, output));
		CHECK(output.mTick == 123);
		CHECK_FALSE(world.FindInteractions(UINT64_MAX, 0, output)); // Interval exceeds exact integer-double range.
		CHECK(output.mTick == 123);
		ObjectState saved = object;
		object.mVelocity[0] = std::numeric_limits<double>::max();
		CHECK_FALSE(object.Evaluate(2, 1, saved));
		CHECK(saved.mVelocity[0] == 0);
		GlobalPosition edge;
		edge.mSector[0] = INT64_MAX; edge.mOffset[0] = GlobalPosition::cSectorSize - 1;
		frame.mOrigin = edge; frame.mVelocity[0] = 2;
		GlobalPosition unchanged;
		CHECK_FALSE(frame.OriginAt(1, 1, unchanged));
		CHECK(unchanged.mSector[0] == 0);
	}

	TEST_CASE("SchedulingHorizonsAreCheckedDirectlyAtArbitraryEpochs")
	{
		constexpr Tick exact_limit = UINT64_C(1) << 53;
		for (Tick tick : {Tick(0), Tick(37), UINT64_MAX - exact_limit - 1})
			for (Tick horizon : {exact_limit - 1, exact_limit, exact_limit + 1})
				for (bool scheduled : {false, true})
				{
					World world(1);
					// Both endpoints are individually within 2^53 of the epochs,
					// even when their separation is 2^53 + 1. This must not bypass
					// the direct requested/effective horizon contract.
					Frame frame; frame.mID = 1; frame.mEpoch = tick + 1;
					ObjectState object; object.mID = 1; object.mFrame = 1;
					object.mEpoch = frame.mEpoch; object.mTier = Tier::Detailed;
					if (scheduled) object.mSchedulingInterval = horizon;
					Tick interval = scheduled? 1 : horizon;
					CHECK(world.AddFrame(frame)); CHECK(world.AddObject(object));
					bool allowed = horizon <= exact_limit;
					GlobalBounds bounds;
					bounds.mMin.mOffset[0] = 123;
					CHECK(BuildSweptBounds(object, frame, frame, tick, interval, 1, bounds) == allowed);
					if (!allowed) CHECK(bounds.mMin.mOffset[0] == 123);
					ActivationResult result; result.mTick = 123; result.mPairs = {{7, 9}};
					CHECK(world.FindInteractions(tick, interval, result) == allowed);
					if (!allowed)
					{
						CHECK(result.mTick == 123);
						CHECK(result.mPairs == std::vector<InteractionPair>{{7, 9}});
						CHECK_FALSE(world.PrepareInteractions(tick, interval, result));
						CHECK(result.mTick == 123);
						CHECK_FALSE(world.TryDemote(1, Tier::Distant, tick + 1, 0, interval));
						CHECK(world.FindObject(1)->mTier == Tier::Detailed);
						CHECK(world.FindObject(1)->mEpoch == tick + 1);
						CHECK(world.GetRepresentationVersion() == 0);
					}
				}
	}

	TEST_CASE("EmptyWorldStillRejectsRequestedHorizonsBeyondExactTickRange")
	{
		constexpr Tick exact_limit = UINT64_C(1) << 53;
		for (Tick tick : {Tick(0), Tick(37), UINT64_MAX - exact_limit - 1})
			for (bool add_frame : {false, true})
			{
				World world(1);
				if (add_frame)
				{
					Frame frame; frame.mID = 1; frame.mEpoch = tick + 1;
					CHECK(world.AddFrame(frame));
				}
				ActivationResult result;
				CHECK(world.FindInteractions(tick, exact_limit - 1, result));
				CHECK(world.FindInteractions(tick, exact_limit, result));
				result.mTick = 123; result.mPairs = {{7, 9}};
				CHECK_FALSE(world.FindInteractions(tick, exact_limit + 1, result));
				CHECK(result.mTick == 123);
				CHECK(result.mPairs == std::vector<InteractionPair>{{7, 9}});
				CHECK_FALSE(world.PrepareInteractions(tick, exact_limit + 1, result));
				CHECK(result.mTick == 123);
				CHECK(world.GetRepresentationVersion() == 0);
			}
	}

	TEST_CASE("TicksNearUint64MaximumUseSmallDifferencesNotAbsoluteDoubleTime")
	{
		Frame frame; frame.mID = 1; frame.mEpoch = UINT64_MAX - 10;
		ObjectState object; object.mID = 1; object.mFrame = 1; object.mEpoch = frame.mEpoch;
		object.mVelocity = {0.001, -0.002, 0.003, -0.004};
		ObjectState result;
		CHECK(object.Evaluate(UINT64_MAX - 1, 0.5, result));
		CHECK(result.mPosition[0] == std::fma(object.mVelocity[0], 4.5, 0.0));
		GlobalBounds bounds;
		CHECK(BuildSweptBounds(object, frame, frame, UINT64_MAX - 1, 1, 0.5, bounds));
		CHECK_FALSE(BuildSweptBounds(object, frame, frame, UINT64_MAX - 1, 2, 0.5, bounds));
	}

	TEST_CASE("SweptBoundsCoverSampledFourAxisAccelerationRotationAndOldEpoch")
	{
		Frame frame; frame.mID = 1;
		frame.mOrigin.mSector = {INT64_MAX - 20, INT64_MIN + 20, INT64_MAX - 20, INT64_MIN + 20};
		frame.mOrigin.mOffset.fill(GlobalPosition::cSectorSize - 0.5);
		frame.mVelocity = {7, -11, 13, -17};
		ObjectState object; object.mID = 1; object.mFrame = 1;
		object.mPosition = {-1, 2, -3, 4}; object.mVelocity = {3, -5, 7, -9};
		object.mAccelerationBound = 4; object.mRotationRadius = 5; object.mWeaponReach = 2;
		object.mSchedulingInterval = 16;
		GlobalBounds bounds;
		CHECK(BuildSweptBounds(object, frame, frame, 4, 1, 0.25, bounds));
		// Independent samples of x = x0 + vt + 1/2 at^2, with |a|=4,
		// plus rotating hull/reach support in every positive/negative axis.
		for (int sample = 0; sample <= 64; ++sample)
		{
			double t = 1 + double(sample) / 16;
			for (size_t axis = 0; axis < 4; ++axis)
				for (double sign : {-1.0, 1.0})
				{
					Vector location;
					for (size_t i = 0; i < 4; ++i) location[i] = object.mPosition[i] + object.mVelocity[i] * t;
					location[axis] += sign * (2 * t * t + 7);
					GlobalPosition point;
					CHECK(frame.mOrigin.Translate(location, point));
					CHECK(bounds.Contains(point));
				}
		}
	}
}
