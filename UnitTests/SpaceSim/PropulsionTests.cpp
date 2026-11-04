// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"

#include <Jolt/SpaceSim/Propulsion.h>
using namespace JPH::SpaceSim;

TEST_SUITE("PropulsionTests")
{
	static Vec4 Axis(int inAxis)
	{
		return Vec4(inAxis == 0? 1.0f : 0.0f, inAxis == 1? 1.0f : 0.0f, inAxis == 2? 1.0f : 0.0f, inAxis == 3? 1.0f : 0.0f);
	}

	static Bivec Plane(int inPlane, float inValue)
	{
		return Bivec(inPlane == 0? inValue : 0, inPlane == 1? inValue : 0, inPlane == 2? inValue : 0,
			inPlane == 3? inValue : 0, inPlane == 4? inValue : 0, inPlane == 5? inValue : 0);
	}

	TEST_CASE("FailedDisabledAndDamagedActuatorsRespectCapacity")
	{
		Thruster t;
		t.mMaxForce = 10;
		t.mForce = 10;
		PropulsionRequest request;
		request.mForce = Vec4(20, 0, 0, 0);
		for (int mode = 0; mode < 3; ++mode)
		{
			t.mEnabled = mode != 0;
			t.mFailed = mode == 1;
			t.mHealth = mode == 2? 0.25 : 1;
			PropulsionResult result;
			CHECK(StepPropulsion(request, { t }, Vec4::sZero(), {}, 1, result));
			CHECK(result.mImpulse.GetX() == (mode == 2? 2.5f : 0.0f));
			CHECK(result.mActuators[0].mThruster.mForce == (mode == 2? 2.5 : 0));
		}
	}

	TEST_CASE("AllFourTranslationAndSixRotationAxesHaveBidirectionalAuthority")
	{
		std::vector<Thruster> rig;
		for (int axis = 0; axis < 4; ++axis)
			for (int sign : { -1, 1 })
			{
				Thruster t;
				t.mID = uint64(rig.size());
				t.mDirection = float(sign) * Axis(axis);
				t.mMaxForce = 10;
				rig.push_back(t);
			}
		for (int i = 0; i < 4; ++i)
			for (int j = i + 1; j < 4; ++j)
				for (int sign : { -1, 1 })
					for (int side : { -1, 1 })
					{
						Thruster t;
						t.mID = uint64(rig.size());
						t.mDirection = float(side) * Axis(i);
						t.mMount = float(side * sign) * Axis(j);
						t.mMaxForce = 10;
						rig.push_back(t);
					}
		for (int axis = 0; axis < 10; ++axis)
			for (int sign : { -1, 1 })
			{
				PropulsionRequest request;
				if (axis < 4) request.mForce = float(sign) * Axis(axis);
				else request.mTorque = Plane(axis - 4, float(sign));
				PropulsionResult result;
				CHECK(StepPropulsion(request, rig, Vec4::sZero(), {}, 1, result));
				CHECK(result.mResidualImpulse.Length() < 1.0e-5f);
				CHECK(result.mResidualAngularImpulse.Length() < 1.0e-5f);
				for (const auto &a : result.mActuators)
				{
					CHECK(a.mTargetForce >= 0);
					CHECK(a.mTargetForce <= 10);
				}
			}
	}

	TEST_CASE("CurrentCOMTorqueSignWeightedCompromiseAndStableIDs")
	{
		Thruster t;
		t.mID = uint64(1) << 48;
		t.mMount = Vec4(3, 3, 5, 7);
		t.mMaxForce = 10;
		PropulsionRequest request;
		request.mTorque = Bivec(1, 0, 0, 0, 0, 0);
		request.mWeights[4] = 3;
		PropulsionResult result;
		CHECK(StepPropulsion(request, { t }, Vec4(3, 2, 5, 7), {}, 2, result));
		CHECK(result.mActuators[0].mTargetForce == doctest::Approx(0.75));
		CHECK(result.mImpulse.GetX() == doctest::Approx(1.5));
		CHECK(result.mAngularImpulse.GetE12() == doctest::Approx(1.5));
		CHECK(result.mAngularImpulse.GetE13() == 0);
		Thruster other = t;
		other.mID = 2;
		PropulsionResult a, b;
		CHECK(StepPropulsion(request, { t, other }, Vec4(3, 2, 5, 7), {}, 2, a));
		CHECK(StepPropulsion(request, { other, t }, Vec4(3, 2, 5, 7), {}, 2, b));
		CHECK(a.mImpulse == b.mImpulse);
		CHECK(a.mAngularImpulse == b.mAngularImpulse);
		CHECK(a.mActuators[0].mThruster.mID == 2);
		CHECK(a.mActuators[1].mThruster.mID == t.mID);
		CHECK(a.mActuators[0].mTargetForce == b.mActuators[0].mTargetForce);
	}

	TEST_CASE("ResponseIntegratesExactlyIncludingSpoolDownAndTinySteps")
	{
		Thruster t;
		t.mMaxForce = 8;
		t.mResponseTime = 0.5;
		PropulsionRequest request;
		request.mForce = Vec4(8, 0, 0, 0);
		PropulsionResult result;
		CHECK(StepPropulsion(request, { t }, Vec4::sZero(), {}, 1, result));
		double reference = 8 * (1 - 0.5 * (1 - std::exp(-2.0)));
		CHECK(result.mActuators[0].mDeliveredImpulse == doctest::Approx(reference));
		CHECK(result.mActuators[0].mThruster.mForce == doctest::Approx(8 * (1 - std::exp(-2.0))));
		t.mForce = 8;
		request.mForce = Vec4::sZero();
		CHECK(StepPropulsion(request, { t }, Vec4::sZero(), {}, 1, result));
		CHECK(result.mActuators[0].mDeliveredImpulse == doctest::Approx(4 * (1 - std::exp(-2.0))));
		t.mForce = 0;
		request.mForce = Vec4(8, 0, 0, 0);
		CHECK(StepPropulsion(request, { t }, Vec4::sZero(), {}, 1.0e-12, result));
		CHECK(result.mActuators[0].mDeliveredImpulse > 0);
		CHECK(std::abs(result.mActuators[0].mDeliveredImpulse / 8.0e-24 - 1) < 1.0e-10);
	}

	TEST_CASE("ResponseTimeRefinementMatchesAnalyticalImpulse")
	{
		for (int count : { 1, 60, 120, 240 })
		{
			Thruster t;
			t.mMaxForce = 8;
			t.mResponseTime = 0.5;
			PropulsionRequest request;
			request.mForce = Vec4(8, 0, 0, 0);
			double total = 0;
			for (int i = 0; i < count; ++i)
			{
				PropulsionResult result;
				CHECK(StepPropulsion(request, { t }, Vec4::sZero(), {}, 1.0 / count, result));
				t = result.mActuators[0].mThruster;
				total += result.mActuators[0].mDeliveredImpulse;
			}
			CHECK(std::abs(total - 8 * (1 - 0.5 * (1 - std::exp(-2.0)))) < 1.0e-12);
		}
	}

	TEST_CASE("ResourcesChargeDeliveredImpulseAndHeatNotCommandedThrust")
	{
		Thruster t;
		t.mMaxForce = 8;
		t.mResponseTime = 0.5;
		t.mPropellantRate = 2;
		t.mPowerDraw = 4;
		t.mHeatRate = 6;
		t.mHeat = 3;
		PropulsionRequest request;
		request.mForce = Vec4(8, 0, 0, 0);
		PropulsionResources resources { 10, 20, 100 };
		PropulsionResult result;
		CHECK(StepPropulsion(request, { t }, Vec4::sZero(), resources, 1, result));
		double fractionSeconds = 1 - 0.5 * (1 - std::exp(-2.0));
		CHECK(result.mPropellantConsumed == doctest::Approx(2 * fractionSeconds));
		CHECK(result.mEnergyConsumed == doctest::Approx(4 * fractionSeconds));
		CHECK(result.mHeatGenerated == doctest::Approx(6 * fractionSeconds));
		CHECK(result.mActuators[0].mThruster.mHeat == doctest::Approx(3 + 6 * fractionSeconds));
		CHECK(result.mResources.mPropellant == doctest::Approx(10 - result.mPropellantConsumed));
		CHECK(result.mResources.mEnergy == doctest::Approx(20 - result.mEnergyConsumed));
		CHECK(t.mHeat == 3);
		CHECK(resources.mPropellant == 10);
	}

	TEST_CASE("MidIntervalFuelExhaustionStopsOnlyConsumersWithoutPhantomImpulse")
	{
		Thruster a;
		a.mID = 1;
		a.mMaxForce = 10;
		a.mPropellantRate = 2;
		a.mPowerDraw = 4;
		a.mHeatRate = 6;
		Thruster b;
		b.mID = 2;
		b.mDirection = Vec4::sAxisY();
		b.mMaxForce = 10;
		b.mPowerDraw = 1;
		b.mHeatRate = 2;
		PropulsionRequest request;
		request.mForce = Vec4(10, 10, 0, 0);
		PropulsionResult result;
		CHECK(StepPropulsion(request, { a, b }, Vec4::sZero(), { 0.5, 10, 100 }, 1, result));
		CHECK_APPROX_EQUAL(result.mImpulse, Vec4(2.5f, 10, 0, 0), 1.0e-6f);
		CHECK(result.mActuators[0].mThruster.mForce == 0);
		CHECK(result.mActuators[1].mThruster.mForce == 10);
		CHECK(result.mPropellantConsumed == doctest::Approx(0.5));
		CHECK(result.mEnergyConsumed == doctest::Approx(2));
		CHECK(result.mHeatGenerated == doctest::Approx(3.5));
		CHECK(result.mPropellantDepleted);
		CHECK_FALSE(result.mEnergyDepleted);
		CHECK(result.mResources.mPropellant >= 0);
	}

	TEST_CASE("SharedFuelAndEnergyAreCoupledAcrossDifferentResponseTimes")
	{
		Thruster a;
		a.mID = 1;
		a.mMaxForce = 10;
		a.mResponseTime = 0.5;
		a.mPropellantRate = 2;
		a.mPowerDraw = 4;
		Thruster b = a;
		b.mID = 2;
		b.mDirection = Vec4::sAxisY();
		b.mResponseTime = 0.25;
		b.mPropellantRate = 3;
		b.mPowerDraw = 2;
		PropulsionRequest request;
		request.mForce = Vec4(10, 10, 0, 0);
		// Choose a known depletion time: both resources expire at exactly t = 0.3.
		double ia = 10 * (0.3 - 0.5 * (1 - std::exp(-0.6)));
		double ib = 10 * (0.3 - 0.25 * (1 - std::exp(-1.2)));
		double fuel = ia * 0.2 + ib * 0.3;
		double energy = ia * 0.4 + ib * 0.2;
		for (int count : { 1, 60, 120, 240 })
		{
			std::vector<Thruster> actuators { a, b };
			PropulsionResources resources { fuel, energy, 100 };
			double totalA = 0, totalB = 0, usedFuel = 0, usedEnergy = 0;
			for (int i = 0; i < count; ++i)
			{
				PropulsionResult result;
				CHECK(StepPropulsion(request, actuators, Vec4::sZero(), resources, 1.0 / count, result));
				resources = result.mResources;
				for (size_t j = 0; j < actuators.size(); ++j) actuators[j] = result.mActuators[j].mThruster;
				totalA += result.mActuators[0].mDeliveredImpulse;
				totalB += result.mActuators[1].mDeliveredImpulse;
				usedFuel += result.mPropellantConsumed;
				usedEnergy += result.mEnergyConsumed;
				CHECK(resources.mPropellant >= 0);
				CHECK(resources.mEnergy >= 0);
			}
			CHECK(std::abs(totalA - ia) < 1.0e-12);
			CHECK(std::abs(totalB - ib) < 1.0e-12);
			CHECK(std::abs(usedFuel - fuel) < 1.0e-12);
			CHECK(std::abs(usedEnergy - energy) < 1.0e-12);
		}
	}

	TEST_CASE("SharedPowerLimitBoundsInstantaneousDeliveryAndDebitsEnergy")
	{
		Thruster a;
		a.mID = 1;
		a.mMaxForce = 10;
		a.mPowerDraw = 10;
		Thruster b = a;
		b.mID = 2;
		b.mDirection = Vec4::sAxisY();
		PropulsionRequest request;
		request.mForce = Vec4(10, 10, 0, 0);
		PropulsionResult result;
		CHECK(StepPropulsion(request, { a, b }, Vec4::sZero(), { 0, 100, 5 }, 2, result));
		CHECK_APPROX_EQUAL(result.mImpulse, Vec4(5, 5, 0, 0), 1.0e-6f);
		CHECK(result.mEnergyConsumed == doctest::Approx(10));
		CHECK(result.mPowerLimited);
		CHECK(result.mActuators[0].mThruster.mForce + result.mActuators[1].mThruster.mForce <= 5);
		// Falling power limits hard-cut preexisting actuator state, even during spool-down.
		a.mForce = 10;
		a.mResponseTime = 1;
		request.mForce = Vec4::sZero();
		CHECK(StepPropulsion(request, { a }, Vec4::sZero(), { 0, 100, 0 }, 1, result));
		CHECK(result.mImpulse == Vec4::sZero());
		CHECK(result.mEnergyConsumed == 0);
		CHECK(result.mActuators[0].mThruster.mForce == 0);
	}

	TEST_CASE("PowerLimitedResponseRefinesToSameAnalyticalSolution")
	{
		for (int count : { 1, 60, 120, 240 })
		{
			Thruster a;
			a.mID = 1;
			a.mMaxForce = 10;
			a.mPowerDraw = 10;
			a.mResponseTime = 0.5;
			Thruster b = a;
			b.mID = 2;
			b.mDirection = Vec4::sAxisY();
			PropulsionRequest request;
			request.mForce = Vec4(10, 10, 0, 0);
			PropulsionResources resources { 0, 100, 5 };
			double total = 0;
			for (int i = 0; i < count; ++i)
			{
				PropulsionResult result;
				CHECK(StepPropulsion(request, { a, b }, Vec4::sZero(), resources, 1.0 / count, result));
				a = result.mActuators[0].mThruster;
				b = result.mActuators[1].mThruster;
				resources = result.mResources;
				total += result.mActuators[0].mDeliveredImpulse;
			}
			CHECK(std::abs(total - 2.5 * (1 - 0.5 * (1 - std::exp(-2.0)))) < 1.0e-12);
		}
	}

	TEST_CASE("TinyResourceBudgetsStillDeliverTheirAccountedSubTickImpulse")
	{
		Thruster t;
		t.mMaxForce = 10;
		t.mPropellantRate = 2;
		PropulsionRequest request;
		request.mForce = Vec4(10, 0, 0, 0);
		PropulsionResult result;
		CHECK(StepPropulsion(request, { t }, Vec4::sZero(), { 1.0e-100, 0, 0 }, 1, result));
		CHECK(std::abs(result.mActuators[0].mDeliveredImpulse / 5.0e-100 - 1) < 1.0e-14);
		CHECK(std::abs(result.mPropellantConsumed / 1.0e-100 - 1) < 1.0e-14);
		CHECK(result.mResources.mPropellant >= 0);
		CHECK(result.mActuators[0].mThruster.mForce == 0);
	}

	TEST_CASE("SequentialFuelThenEnergyEventsKeepIndependentEnginesRunning")
	{
		Thruster fueled;
		fueled.mID = 1;
		fueled.mMaxForce = 10;
		fueled.mPropellantRate = 2;
		fueled.mPowerDraw = 4;
		Thruster electric;
		electric.mID = 2;
		electric.mMaxForce = 10;
		electric.mDirection = Vec4::sAxisY();
		electric.mPowerDraw = 2;
		Thruster free;
		free.mID = 3;
		free.mMaxForce = 10;
		free.mDirection = Vec4::sAxisW();
		PropulsionRequest request;
		request.mForce = Vec4(10, 10, 0, 10);
		PropulsionResult result;
		CHECK(StepPropulsion(request, { fueled, electric, free }, Vec4::sZero(), { 0.5, 2, 100 }, 1, result));
		CHECK_APPROX_EQUAL(result.mImpulse, Vec4(2.5f, 5, 0, 10), 1.0e-6f);
		CHECK(result.mPropellantConsumed == doctest::Approx(0.5));
		CHECK(result.mEnergyConsumed == doctest::Approx(2));
		CHECK(result.mPropellantDepleted);
		CHECK(result.mEnergyDepleted);
		CHECK(result.mActuators[0].mThruster.mForce == 0);
		CHECK(result.mActuators[1].mThruster.mForce == 0);
		CHECK(result.mActuators[2].mThruster.mForce == 10);
		// Energy-only shortage is not excused by fuel remaining.
		CHECK(StepPropulsion(request, { fueled, electric, free }, Vec4::sZero(), { 100, 1.5, 100 }, 1, result));
		CHECK_APPROX_EQUAL(result.mImpulse, Vec4(2.5f, 2.5f, 0, 10), 1.0e-6f);
		CHECK(result.mPropellantConsumed == doctest::Approx(0.5));
		CHECK_FALSE(result.mPropellantDepleted);
		CHECK(result.mEnergyDepleted);
	}

	TEST_CASE("ResourceLimitedImpulsesMatchIndependentScalarWedgeAndStablePermutation")
	{
		std::vector<Thruster> rig;
		for (int i = 0; i < 8; ++i)
		{
			Thruster t;
			t.mID = (uint64(1) << 50) + uint64(i);
			t.mMaxForce = 10;
			t.mDirection = Vec4(float(i - 3), float(2 - i), float(i % 3 + 1), float(i % 2 - 1)).Normalized();
			t.mMount = Vec4(float(i), float(i % 2), float(i % 3 - 1), float(i % 4));
			t.mResponseTime = 0.1 * (i + 1);
			t.mPropellantRate = 0.5 * (i % 2);
			t.mPowerDraw = i + 1;
			t.mHeatRate = 2;
			rig.push_back(t);
		}
		PropulsionRequest request;
		request.mForce = Vec4(30, -20, 10, 40);
		request.mTorque = Bivec(2, -3, 4, 5, -6, 7);
		Vec4 com(0.1f, -0.2f, 0.3f, -0.4f);
		PropulsionResult a, b;
		CHECK(StepPropulsion(request, rig, com, { 0.3, 0.4, 3 }, 1, a));
		std::reverse(rig.begin(), rig.end());
		CHECK(StepPropulsion(request, rig, com, { 0.3, 0.4, 3 }, 1, b));
		CHECK(a.mImpulse == b.mImpulse);
		CHECK(a.mAngularImpulse == b.mAngularImpulse);
		CHECK(a.mResources.mPropellant == b.mResources.mPropellant);
		CHECK(a.mResources.mEnergy == b.mResources.mEnergy);
		CHECK(a.mHeatGenerated == b.mHeatGenerated);
		double torque[6] {}, force[4] {}, fuel = 0, energy = 0, heat = 0;
		for (size_t index = 0; index < a.mActuators.size(); ++index)
		{
			const auto &actuator = a.mActuators[index];
			const auto &t = actuator.mThruster;
			CHECK(actuator.mDeliveredImpulse == b.mActuators[index].mDeliveredImpulse);
			double seconds = actuator.mDeliveredImpulse / t.mMaxForce;
			fuel += seconds * t.mPropellantRate;
			energy += seconds * t.mPowerDraw;
			heat += seconds * t.mHeatRate;
			int plane = 0;
			Vec4 r = t.mMount - com;
			for (int i = 0; i < 4; ++i)
			{
				force[i] += actuator.mDeliveredImpulse * t.mDirection[uint(i)];
				for (int j = i + 1; j < 4; ++j)
					torque[plane++] += actuator.mDeliveredImpulse * (double(t.mDirection[uint(i)]) * r[uint(j)] - double(t.mDirection[uint(j)]) * r[uint(i)]);
			}
		}
		for (uint i = 0; i < 4; ++i) CHECK(std::abs(a.mImpulse[i] - force[i]) < 1.0e-6);
		for (uint i = 0; i < 6; ++i) CHECK(std::abs(a.mAngularImpulse[i] - torque[i]) < 1.0e-6);
		CHECK(std::abs(a.mPropellantConsumed - fuel) < 1.0e-14);
		CHECK(std::abs(a.mEnergyConsumed - energy) < 1.0e-14);
		CHECK(std::abs(a.mHeatGenerated - heat) < 1.0e-14);
	}

	TEST_CASE("ResponseRoundoffCannotExceedActuatorCapacity")
	{
		Thruster t;
		t.mMaxForce = 10;
		t.mForce = 10;
		t.mResponseTime = 1;
		PropulsionRequest request;
		request.mForce = Vec4(10, 0, 0, 0);
		PropulsionResult result;
		CHECK(StepPropulsion(request, { t }, Vec4::sZero(), {}, 0.0004, result));
		CHECK(result.mActuators[0].mThruster.mForce <= 10);
		CHECK(result.mActuators[0].mDeliveredImpulse <= 10 * 0.0004);
	}

	TEST_CASE("InvalidInputAndNumericFailureLeaveOutputUnchanged")
	{
		Thruster valid;
		valid.mID = 1;
		valid.mMaxForce = 10;
		PropulsionRequest request;
		request.mForce = Vec4(10, 0, 0, 0);
		PropulsionResult output;
		output.mImpulse = Vec4(1, 2, 3, 4);
		output.mHeatGenerated = 123;
		output.mResources.mPropellant = 456;
		output.mActuators.push_back({ valid, 7, 8 });
		auto reject = [&](const PropulsionRequest &r, const std::vector<Thruster> &ts, Vec4Arg com, const PropulsionResources &resources, double dt)
		{
			CHECK_FALSE(StepPropulsion(r, ts, com, resources, dt, output));
			CHECK(output.mImpulse == Vec4(1, 2, 3, 4));
			CHECK(output.mHeatGenerated == 123);
			CHECK(output.mResources.mPropellant == 456);
			CHECK(output.mActuators.size() == 1);
			CHECK(output.mActuators[0].mTargetForce == 7);
		};
		const double nan = std::numeric_limits<double>::quiet_NaN();
		const double inf = std::numeric_limits<double>::infinity();
		for (double bad : { -1.0, nan, inf })
		{
			reject(request, { valid }, Vec4::sZero(), {}, bad);
			reject(request, { valid }, Vec4::sZero(), { bad, 0, 0 }, 1);
			reject(request, { valid }, Vec4::sZero(), { 0, bad, 0 }, 1);
			reject(request, { valid }, Vec4::sZero(), { 0, 0, bad }, 1);
			for (double Thruster::*field : { &Thruster::mMaxForce, &Thruster::mResponseTime, &Thruster::mPropellantRate,
				&Thruster::mPowerDraw, &Thruster::mHeatRate, &Thruster::mHealth, &Thruster::mForce, &Thruster::mHeat })
			{
				Thruster t = valid;
				t.*field = bad;
				reject(request, { t }, Vec4::sZero(), {}, 1);
			}
		}
		Thruster t = valid;
		t.mHealth = 1.1;
		reject(request, { t }, Vec4::sZero(), {}, 1);
		t = valid;
		t.mDirection = Vec4(2, 0, 0, 0);
		reject(request, { t }, Vec4::sZero(), {}, 1);
		t.mDirection = Vec4::sZero();
		reject(request, { t }, Vec4::sZero(), {}, 1);
		for (int axis = 0; axis < 4; ++axis)
		{
			Vec4 bad(axis == 0? float(nan) : 0, axis == 1? float(nan) : 0, axis == 2? float(nan) : 0, axis == 3? float(nan) : 0);
			reject(request, { valid }, bad, {}, 1);
			t = valid;
			t.mMount = bad;
			reject(request, { t }, Vec4::sZero(), {}, 1);
			t.mMount = Vec4::sZero();
			t.mDirection = bad;
			reject(request, { t }, Vec4::sZero(), {}, 1);
			PropulsionRequest r = request;
			r.mForce = bad;
			reject(r, { valid }, Vec4::sZero(), {}, 1);
		}
		for (int axis = 0; axis < 6; ++axis)
		{
			PropulsionRequest r = request;
			r.mTorque = Plane(axis, float(inf));
			reject(r, { valid }, Vec4::sZero(), {}, 1);
		}
		for (double bad : { -1.0, nan, inf })
		{
			PropulsionRequest r = request;
			r.mWeights[9] = bad;
			reject(r, { valid }, Vec4::sZero(), {}, 1);
		}
		reject(request, { valid, valid }, Vec4::sZero(), {}, 1);
		std::vector<Thruster> tooMany(4097, valid);
		for (size_t i = 0; i < tooMany.size(); ++i) tooMany[i].mID = uint64(i);
		reject(request, tooMany, Vec4::sZero(), {}, 1);
		reject(request, { valid }, Vec4::sZero(), {}, std::numeric_limits<double>::max());
	}

	TEST_CASE("EmptyZeroDurationAndUnreachableRequestsHaveNoImplicitForces")
	{
		PropulsionRequest request;
		request.mForce = Vec4(1, 2, 3, 4);
		request.mTorque = Bivec(1, 2, 3, 4, 5, 6);
		PropulsionResult result;
		CHECK(StepPropulsion(request, {}, Vec4::sZero(), {}, 2, result));
		CHECK(result.mImpulse == Vec4::sZero());
		CHECK(result.mAngularImpulse == Bivec::sZero());
		CHECK(result.mResidualImpulse == 2 * request.mForce);
		CHECK(result.mResidualAngularImpulse == 2 * request.mTorque);
		Thruster t;
		t.mMaxForce = 10;
		t.mForce = 3;
		t.mResponseTime = 1;
		CHECK(StepPropulsion({}, { t }, Vec4::sZero(), {}, 0, result));
		CHECK(result.mImpulse == Vec4::sZero());
		CHECK(result.mActuators[0].mThruster.mForce == 3);
		CHECK(result.mEnergyConsumed == 0);
	}

	TEST_CASE("TranslationSaturatesAndReportsResidual")
	{
		Thruster t;
		t.mID = 17;
		t.mMaxForce = 4;
		PropulsionRequest request;
		request.mForce = Vec4(10, 0, 0, 0);
		PropulsionResult result;
		CHECK(StepPropulsion(request, { t }, Vec4::sZero(), {}, 0.5, result));
		CHECK_APPROX_EQUAL(result.mImpulse, Vec4(2, 0, 0, 0), 1.0e-6f);
		CHECK_APPROX_EQUAL(result.mResidualImpulse, Vec4(3, 0, 0, 0), 1.0e-6f);
		CHECK(result.mActuators[0].mTargetForce == 4);
		CHECK(result.mActuators[0].mThruster.mForce == 4);
	}
}
