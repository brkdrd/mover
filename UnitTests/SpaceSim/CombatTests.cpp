// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include <Jolt/Physics/Body/Body.h>
#include <Jolt/SpaceSim/Combat.h>

using namespace JPH::SpaceSim;

TEST_SUITE("CombatTests")
{
	struct Scene
	{
		CombatBox box;
		CombatCompound target;
		CombatProjectile projectile;
		CombatQuery query;
		Scene() { target.mBodyID = 10; box.mComponentID = 20; target.mBoxes = &box; target.mBoxCount = 1; }
		CombatHit Sweep() const { return SweepProjectile(projectile, &target, 1, query); }
	};

	TEST_CASE("LongRotatingHullCrossesBetweenSeparatedEndpointPoses")
	{
		Scene s;
		s.box.mHalfExtent = Vec4(10, 0.05f, 0.05f, 0.05f);
		s.target.mMotion.mAngularVelocity = Bivec(100 * JPH_PI, 0, 0, 0, 0, 0);
		s.projectile.mMotion.mPosition = Vec4(0, -8, 0, 0);
		s.projectile.mRadius = 0.1;
		s.query.mDuration = 0.01;
		auto hit = s.Sweep();
		CHECK(hit.mStatus == CombatQueryStatus::Hit);
		CHECK(hit.mTime > 0.004);
		CHECK(hit.mTime < 0.005);
		CHECK(hit.mNormal.Length() == doctest::Approx(1).epsilon(0.0001));
		// Repeat the same geometry in a plane involving W, not a 3D projection.
		s.target.mMotion.mAngularVelocity = Bivec(0, 0, 100 * JPH_PI, 0, 0, 0);
		s.projectile.mMotion.mPosition = Vec4(0, 0, 0, -8);
		auto w_hit = s.Sweep();
		CHECK(w_hit.mStatus == CombatQueryStatus::Hit);
		CHECK(std::abs(w_hit.mTime - hit.mTime) < 1.0e-7);
	}

	TEST_CASE("FastTargetRelativeTranslationAndFourthDimension")
	{
		Scene s;
		s.target.mMotion.mPosition = Vec4(0, 0, 0, -10);
		s.target.mMotion.mVelocity = Vec4(0, 0, 0, 20000);
		s.query.mDuration = 0.001;
		s.projectile.mRadius = 0.25;
		auto hit = s.Sweep();
		CHECK(hit.mStatus == CombatQueryStatus::Hit);
		CHECK(std::abs(hit.mTime - 0.0004625) < 1.0e-7);
		CHECK(hit.mNormal.GetW() == doctest::Approx(1));
		s.projectile.mMotion.mPosition = Vec4(0, 0, 2, 0);
		CHECK(s.Sweep().mStatus == CombatQueryStatus::NoHit);
	}

	TEST_CASE("FiniteSphereGrazesBoxCornerWithoutInflatedBoxFalsePositive")
	{
		Scene s;
		s.projectile.mRadius = 0.25;
		s.projectile.mMotion.mPosition = Vec4(-2, 0.75f, 0, 0);
		s.projectile.mMotion.mVelocity = Vec4(400, 0, 0, 0);
		s.query.mDuration = 0.01;
		CHECK(s.Sweep().mStatus == CombatQueryStatus::Hit);
		s.projectile.mMotion.mPosition = Vec4(-2, 0.7f, 0.7f, 0);
		CHECK(s.Sweep().mStatus == CombatQueryStatus::NoHit);
	}

	TEST_CASE("RejectsIncoherentFramesEpochsAndInvalidGeometry")
	{
		Scene s;
		s.target.mMotion.mFrameID = 1;
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
		CHECK(s.Sweep().mFailure == CombatQueryFailure::FrameEpochMismatch);
		s.target.mMotion.mFrameID = 0;
		s.projectile.mMotion.mEpoch = 1;
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
		s.projectile.mMotion.mEpoch = 0;
		s.box.mHalfExtent = Vec4(-1, 1, 1, 1);
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
		s.box.mHalfExtent = Vec4::sReplicate(0.5f);
		s.box.mLocalRotation = Rotor(1, 0, 0, 0, 0, 0, 0, 1);
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
		s.box.mLocalRotation = Rotor::sIdentity();
		s.projectile.mRadius = std::numeric_limits<double>::quiet_NaN();
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
		s.projectile.mRadius = 0.1;
		s.query.mDuration = -1;
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
		s.query.mDuration = 0.01;
		s.projectile.mMotion.mVelocity = Vec4(1.0e9f, 0, 0, 0);
		CHECK(s.Sweep().mStatus == CombatQueryStatus::Indeterminate);
		CHECK(s.Sweep().mFailure == CombatQueryFailure::OutsideEnvelope);
	}

	TEST_CASE("StableComponentAndBodyIDsBreakSimultaneousTiesNotArrayOrder")
	{
		Scene s;
		CombatBox boxes[2] = {s.box, s.box};
		boxes[0].mComponentID = 200; boxes[1].mComponentID = 100;
		s.target.mBoxes = boxes; s.target.mBoxCount = 2;
		CHECK(s.Sweep().mComponentID == 100);
		std::swap(boxes[0], boxes[1]);
		CHECK(s.Sweep().mComponentID == 100);
		CombatCompound bodies[2] = {s.target, s.target};
		bodies[0].mBodyID = 200; bodies[1].mBodyID = 100;
		CHECK(SweepProjectile(s.projectile,bodies,2,s.query).mBodyID == 100);
		std::swap(boxes[0],boxes[1]);
		s.projectile.mMotion.mPosition = Vec4(-2, 0, 0, 0);
		s.projectile.mMotion.mVelocity = Vec4(400, 0, 0, 0);
		CHECK(s.Sweep().mComponentID == 100);
		boxes[0].mComponentID = boxes[1].mComponentID;
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
	}

	TEST_CASE("CurvedConstantAccelerationCrossesAndReturnsRatherThanEndChord")
	{
		Scene s;
		s.projectile.mMotion.mPosition = Vec4(-2, 0, 0, 0);
		s.projectile.mMotion.mVelocity = Vec4(800, 0, 0, 0);
		s.projectile.mMotion.mAcceleration = Vec4(-160000, 0, 0, 0);
		s.query.mDuration = 0.01;
		auto hit = s.Sweep();
		CHECK(hit.mStatus == CombatQueryStatus::Hit);
		CHECK(hit.mTime > 0.002);
		CHECK(hit.mTime < 0.003);
		// Same relative curve produced by target acceleration instead.
		s.target.mMotion.mAcceleration = -s.projectile.mMotion.mAcceleration;
		s.projectile.mMotion.mAcceleration = Vec4::sZero();
		CHECK(std::abs(s.Sweep().mTime-hit.mTime) < 1.0e-7);
	}

	TEST_CASE("BudgetExhaustionIsNeverSilentMissOrProvisionalHit")
	{
		Scene s;
		s.projectile.mMotion.mPosition = Vec4(-2, 0, 0, 0);
		s.projectile.mMotion.mVelocity = Vec4(400, 0, 0, 0);
		s.query.mMaxDistanceEvaluations = 1;
		auto hit = s.Sweep();
		CHECK(hit.mStatus == CombatQueryStatus::Indeterminate);
		CHECK(hit.mFailure == CombatQueryFailure::BudgetExceeded);
		CHECK(hit.mTelemetry.mDistanceEvaluations == 1);
		CHECK(hit.mTelemetry.mBudgetFailures == 1);
		// A hit on one component cannot prove earliest hit when another is unresolved.
		CombatBox boxes[2] = {s.box,s.box};
		boxes[1].mComponentID = 21;
		s.target.mBoxes = boxes; s.target.mBoxCount = 2;
		s.projectile.mMotion.mPosition = Vec4::sZero();
		hit = s.Sweep();
		CHECK(hit.mStatus == CombatQueryStatus::Indeterminate);
		CHECK(hit.mBodyID == 0);
		CHECK(hit.mComponentID == 0);
		CHECK(hit.mNormal == Vec4::sZero());
		CHECK(hit.mNumericalAllowance == 0);
	}

	TEST_CASE("SweptBroadPhaseRejectsWithoutSpendingNarrowPhaseBudget")
	{
		Scene s;
		s.target.mMotion.mPosition = Vec4(0, 0, 0, 100);
		s.query.mMaxDistanceEvaluations = 0;
		auto miss = s.Sweep();
		CHECK(miss.mStatus == CombatQueryStatus::NoHit);
		CHECK(miss.mTelemetry.mBroadPhaseRejected == 1);
		CHECK(miss.mTelemetry.mDistanceEvaluations == 0);
		CHECK(miss.mTelemetry.mCandidateBoxes == 0);
	}

	TEST_CASE("CompoundOffsetAndLocalOrientationAreAuthoritative")
	{
		Scene s;
		s.box.mLocalCenter = Vec4(8, 0, 0, 0);
		s.box.mLocalRotation = Rotor::sRotation(Vec4::sAxisZ(),Vec4::sAxisW(),0.7f);
		s.box.mHalfExtent = Vec4(0.05f, 0.05f, 0.1f, 0.05f);
		s.target.mMotion.mAngularVelocity = Bivec(100*JPH_PI, 0, 0, 0, 0, 0);
		s.projectile.mMotion.mPosition = Vec4(0, -8, 0, 0);
		s.query.mDuration = 0.01;
		auto hit = s.Sweep();
		CHECK(hit.mStatus == CombatQueryStatus::Hit);
		CHECK(hit.mTelemetry.mCandidateBoxes == 1);
		CHECK(hit.mTime > 0.004);
		CHECK(hit.mTime < 0.005);
		// Static local rotation turns the long X dimension into W.
		s.target.mMotion.mAngularVelocity = Bivec::sZero();
		s.box.mLocalCenter = Vec4::sZero();
		s.box.mHalfExtent = Vec4(4, 0.1f, 0.1f, 0.1f);
		s.box.mLocalRotation = Rotor::sRotation(Vec4::sAxisX(),Vec4::sAxisW(),JPH_PI/2);
		s.projectile.mMotion.mPosition = Vec4(0, 0, 0, -3);
		CHECK(s.Sweep().mStatus == CombatQueryStatus::Hit);
	}

	TEST_CASE("LaunchInheritsRotatingPointVelocityAndConservesEjectionMomentum")
	{
		CombatMotion ship;
		ship.mVelocity = Vec4(4, 5, 6, 7);
		ship.mPosition = Vec4(20, 30, 40, 50);
		ship.mAngularVelocity = Bivec(1, 2, 3, 4, 5, 6);
		Vec4 offset(1, 2, 3, 4), delta(10, -20, 30, -40);
		const float mass = 2;
		auto launch = PrepareCombatLaunch(ship,offset,delta,mass);
		CHECK(launch.mValid);
		Vec4 inherited = ship.mVelocity + ship.mAngularVelocity.Contract(offset);
		CHECK_APPROX_EQUAL(launch.mPosition,ship.mPosition+offset,1.0e-5f);
		CHECK_APPROX_EQUAL(launch.mProjectileVelocity,inherited+delta,1.0e-5f);
		Vec4 momentum_change = mass*(launch.mProjectileVelocity-inherited);
		CHECK_APPROX_EQUAL(momentum_change+launch.mRecoilImpulse,Vec4::sZero(),1.0e-5f);
		Bivec angular_change = Bivec::sWedge(momentum_change,offset);
		CHECK((angular_change+launch.mRecoilAngularImpulse).Length() < 1.0e-5f);
		CHECK(!PrepareCombatLaunch(ship,offset,delta,-1).mValid);
	}

	TEST_CASE("SolverOwnedDamageCannotRequestASecondImpulse")
	{
		Vec4 incoming(10, 0, 0, -4), outgoing(2, 0, 0, -1), offset(0, 2, 0, 3);
		auto solver = PrepareCombatImpact(CombatResponseOwner::PhysicsSolver,2,incoming,outgoing,offset);
		CHECK(solver.mValid);
		CHECK(!solver.mApplyImpulse);
		CHECK(solver.mTargetImpulse == Vec4::sZero());
		CHECK(solver.mTargetAngularImpulse == Bivec::sZero());
		auto event = PrepareCombatImpact(CombatResponseOwner::QueryEvent,2,incoming,outgoing,offset);
		CHECK(event.mValid);
		CHECK(event.mApplyImpulse);
		CHECK_APPROX_EQUAL(event.mTargetImpulse,2*(incoming-outgoing),1.0e-5f);
		CHECK((event.mTargetAngularImpulse-Bivec::sWedge(event.mTargetImpulse,offset)).Length() < 1.0e-5f);
		CHECK(event.mProjectileEnergyLoss == doctest::Approx(incoming.LengthSq()-outgoing.LengthSq()));
		CHECK(!PrepareCombatImpact(CombatResponseOwner::QueryEvent,0,incoming,outgoing,offset).mValid);
	}

	TEST_CASE("PrecisionBudgetRejectsUnreliableWitnessRatherThanPretendingHit")
	{
		Scene s;
		auto hit = s.Sweep();
		CHECK(hit.mStatus == CombatQueryStatus::Hit);
		CHECK(hit.mNumericalAllowance > 0);
		CHECK(hit.mNumericalAllowance <= s.query.mMaxNumericalAllowance);
		s.query.mMaxNumericalAllowance = 1.0e-9;
		CHECK(s.Sweep().mStatus == CombatQueryStatus::Indeterminate);
		CHECK(s.Sweep().mFailure == CombatQueryFailure::PrecisionExceeded);
		s.query.mMaxNumericalAllowance = 0.001;
		s.target.mMotion.mPosition = Vec4(500000, 0, 0, 0);
		s.projectile.mMotion.mPosition = s.target.mMotion.mPosition;
		CHECK(s.Sweep().mStatus == CombatQueryStatus::Indeterminate);
	}

	TEST_CASE("MovingIdenticalComponentTiesRemainStableAcrossRounding")
	{
		Scene s;
		CombatBox boxes[2] = {s.box,s.box};
		s.target.mBoxes = boxes; s.target.mBoxCount = 2;
		s.query.mDistanceTolerance = 1.0e-6;
		for (int k=1; k<=128; ++k)
		{
			CAPTURE(k);
			s.projectile.mMotion.mPosition = Vec4(-0.501f-float(k)*0.017f, 0, 0, 0);
			s.projectile.mMotion.mVelocity = Vec4(321, 0, 0, 0);
			boxes[0].mComponentID = 200; boxes[1].mComponentID = 100;
			auto a = s.Sweep();
			std::swap(boxes[0],boxes[1]);
			auto b = s.Sweep();
			CHECK(a.mStatus == CombatQueryStatus::Hit);
			CHECK(a.mComponentID == 100);
			CHECK(b.mComponentID == a.mComponentID);
			CHECK(b.mTime == a.mTime);
		}
	}

	TEST_CASE("FourDimensionalRoundedCornerOracleAcrossGeneralRotors")
	{
		Scene s;
		s.target.mMotion.mPosition = Vec4(100, -50, 25, 75);
		s.box.mLocalCenter = Vec4(1, -2, 3, -4);
		for (int k=0; k<64; ++k)
		{
			CAPTURE(k);
			float f = float(k)*0.03125f;
			s.target.mMotion.mRotation = Rotor::sExp(Bivec(f, -0.7f*f, 0.3f*f, 0.2f*f, 0.8f*f, -0.5f*f));
			s.box.mLocalRotation = Rotor::sExp(Bivec(-0.3f*f, 0.1f*f, f, -0.8f*f, 0.5f*f, 0.2f*f));
			Rotor rotation = s.target.mMotion.mRotation*s.box.mLocalRotation;
			s.projectile.mMotion.mPosition = s.target.mMotion.mPosition+s.target.mMotion.mRotation*s.box.mLocalCenter+rotation*Vec4::sReplicate(0.6f);
			// Four 0.1 excess coordinates give Euclidean distance 0.2, independently of orientation.
			s.projectile.mRadius = 0.21;
			auto hit = s.Sweep();
			CHECK(hit.mStatus == CombatQueryStatus::Hit);
			CHECK(std::abs(hit.mSeparation+0.01) < 5.0e-5);
			CHECK(hit.mNormal.Dot(rotation*Vec4::sReplicate(0.5f)) > 0.999f);
			s.projectile.mRadius = 0.19;
			CHECK(s.Sweep().mStatus == CombatQueryStatus::NoHit);
		}
	}

	TEST_CASE("IsolatedTangentCornerContactIsNotSkipped")
	{
		Scene s;
		s.projectile.mRadius = 0.25;
		float tangent = 0.5f+0.25f/std::sqrt(2.0f);
		s.projectile.mMotion.mPosition = Vec4(tangent-1,tangent+1,0,0);
		s.projectile.mMotion.mVelocity = Vec4(200,-200,0,0);
		s.query.mDuration = 0.01;
		auto hit = s.Sweep();
		CHECK(hit.mStatus == CombatQueryStatus::Hit);
		CHECK(hit.mTime > 0.004);
		CHECK(hit.mTime <= 0.005);
		s.projectile.mMotion.mPosition += Vec4(0.01f,0.01f,0,0);
		CHECK(s.Sweep().mStatus == CombatQueryStatus::NoHit);
	}

	TEST_CASE("FiniteMassQueryResponseGivesAsteroidRecoilUsingCoreBodyImpulse")
	{
		PhysicsTestContext c;
		c.ZeroGravity();
		Body &target = c.CreateBox(RVec4::sZero(),Rotor::sIdentity(),EMotionType::Dynamic,EMotionQuality::Discrete,Layers::MOVING,Vec4::sReplicate(1));
		float mass = 1/target.GetMotionProperties()->GetInverseMass();
		const float projectile_mass = 2;
		Vec4 incoming(0,0,0,100);
		// Explicit elastic 1D model, NOT a general ramming/rotational CCD claim.
		Vec4 outgoing = incoming*((projectile_mass-mass)/(projectile_mass+mass));
		auto plan = PrepareCombatImpact(CombatResponseOwner::QueryEvent,projectile_mass,incoming,outgoing,Vec4::sZero());
		CHECK(plan.mValid);
		CHECK(plan.mApplyImpulse);
		target.AddImpulse(plan.mTargetImpulse);
		CHECK(target.GetLinearVelocity().GetW() > 0);
		CHECK_APPROX_EQUAL(mass*target.GetLinearVelocity()+projectile_mass*outgoing,projectile_mass*incoming,0.001f);
		float energy_before = 0.5f*projectile_mass*incoming.LengthSq();
		float energy_after = 0.5f*projectile_mass*outgoing.LengthSq()+0.5f*mass*target.GetLinearVelocity().LengthSq();
		CHECK(std::abs(energy_after-energy_before) < 0.01f);
	}

	TEST_CASE("EmptyQueriesNullStorageAndNonfiniteLanesAreExplicit")
	{
		Scene s;
		CHECK(SweepProjectile(s.projectile,nullptr,0,s.query).mStatus == CombatQueryStatus::NoHit);
		CHECK(SweepProjectile(s.projectile,nullptr,1,s.query).mStatus == CombatQueryStatus::InvalidInput);
		s.target.mBoxes = nullptr;
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
		s.target.mBoxes = &s.box;
		for (uint axis=0; axis<4; ++axis)
		{
			s.projectile.mMotion.mPosition = Vec4::sZero();
			s.projectile.mMotion.mPosition[axis] = std::numeric_limits<float>::infinity();
			CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
		}
		s.projectile.mMotion.mPosition = Vec4::sZero();
		s.target.mMotion.mAngularVelocity = Bivec(Vec8(0,0,0,0,0,0,1,0));
		CHECK(s.Sweep().mStatus == CombatQueryStatus::InvalidInput);
	}

	TEST_CASE("EarliestComponentWinsRatherThanLowestIdentity")
	{
		Scene s;
		CombatBox boxes[2] = {s.box,s.box};
		boxes[0].mComponentID = 1; boxes[0].mLocalCenter = Vec4(1,0,0,0);
		boxes[1].mComponentID = 200;
		s.target.mBoxes = boxes; s.target.mBoxCount = 2;
		s.projectile.mMotion.mPosition = Vec4(-2,0,0,0);
		s.projectile.mMotion.mVelocity = Vec4(400,0,0,0);
		CHECK(s.Sweep().mComponentID == 200);
		std::swap(boxes[0],boxes[1]);
		CHECK(s.Sweep().mComponentID == 200);
	}

	TEST_CASE("EndpointSeparatedProjectileCrossesAuthoritativeBox")
	{
		CombatBox box;
		box.mComponentID = 42;
		CombatCompound target;
		target.mBodyID = 7;
		target.mBoxes = &box;
		target.mBoxCount = 1;
		CombatProjectile projectile;
		projectile.mMotion.mPosition = Vec4(-10, 0, 0, 0);
		projectile.mMotion.mVelocity = Vec4(20000, 0, 0, 0);
		projectile.mRadius = 0.25;
		CombatQuery query;
		query.mDuration = 0.001;
		auto result = SweepProjectile(projectile, &target, 1, query);
		CHECK(result.mStatus == CombatQueryStatus::Hit);
		CHECK(result.mBodyID == 7);
		CHECK(result.mComponentID == 42);
		CHECK(result.mTime == doctest::Approx(0.0004625).epsilon(0.001));
		CHECK(result.mNormal.GetX() == doctest::Approx(-1));
		CHECK(result.mPosition.GetX() == doctest::Approx(-0.5).epsilon(0.001));
		CHECK(result.mTelemetry.mDistanceEvaluations > 0);
	}
}
