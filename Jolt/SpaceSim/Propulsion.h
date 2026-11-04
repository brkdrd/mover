// SPDX-License-Identifier: MIT
#pragma once

#include <Jolt/Jolt.h>
#include <Jolt/Math/Bivec.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

JPH_NAMESPACE_BEGIN
namespace SpaceSim {

/// SI-like consistent units. All geometry and requests are in the same ship-local frame.
struct Thruster
{
	uint64 mID = 0;
	Vec4 mMount = Vec4::sZero();
	Vec4 mDirection = Vec4::sAxisX(); ///< Unit direction of force on the ship, NOT exhaust.
	double mMaxForce = 0;
	double mResponseTime = 0;
	double mPropellantRate = 0; ///< Mass / second at nominal maximum force.
	double mPowerDraw = 0; ///< Energy / second at nominal maximum force.
	double mHeatRate = 0; ///< Heat / second at nominal maximum force.
	double mHealth = 1; ///< Linear capacity derating in [0, 1].
	bool mEnabled = true;
	bool mFailed = false;
	double mForce = 0; ///< Previous delivered force magnitude (actuator state).
	double mHeat = 0; ///< Accumulated heat; no implicit cooling or thermal damage.
};

struct PropulsionResources
{
	double mPropellant = 0;
	double mEnergy = 0;
	double mPowerLimit = std::numeric_limits<double>::max();
};

struct PropulsionRequest
{
	Vec4 mForce = Vec4::sZero();
	Bivec mTorque = Bivec::sZero();
	std::array<double, 10> mWeights { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
};

struct PropulsionActuatorResult
{
	Thruster mThruster;
	double mTargetForce = 0; ///< Boundary allocation after power limiting, before depletion. Compare with health * maximum for saturation.
	double mDeliveredImpulse = 0; ///< Integral of nonnegative delivered force, in double precision.
};

struct PropulsionResult
{
	Vec4 mImpulse = Vec4::sZero();
	Bivec mAngularImpulse = Bivec::sZero();
	Vec4 mResidualImpulse = Vec4::sZero();
	Bivec mResidualAngularImpulse = Bivec::sZero();
	Vec4 mAllocationResidualForce = Vec4::sZero();
	Bivec mAllocationResidualTorque = Bivec::sZero();
	PropulsionResources mResources;
	std::vector<PropulsionActuatorResult> mActuators; ///< Ascending stable ID order.
	double mPropellantConsumed = 0; ///< Resource mass-flux ledger ONLY; never changes body mass.
	double mEnergyConsumed = 0;
	double mHeatGenerated = 0;
	bool mPowerLimited = false;
	bool mPropellantDepleted = false;
	bool mEnergyDepleted = false;
};

namespace PropulsionDetail {
using Wrench = std::array<double, 10>;
inline bool Nonnegative(double inValue)
{
	return std::isfinite(inValue) && inValue >= 0;
}
// Deliberately bounded physical inputs, not a promise to support every binary64 exponent.
inline bool Physical(double inValue)
{
	return inValue == 0 || (inValue >= 1.0e-100 && inValue <= 1.0e100);
}
// Reject unsupported intermediate ranges instead of erasing a factor that a later
// multiply could recover. NaN propagates to the step's transactional validation.
inline double Product(double inA, double inB)
{
	if (!Nonnegative(inA) || !Nonnegative(inB)) return std::numeric_limits<double>::quiet_NaN();
	if (inA == 0 || inB == 0) return 0;
	double value = inA * inB;
	return std::isnormal(value)? value : std::numeric_limits<double>::quiet_NaN();
}
inline double Quotient(double inA, double inB)
{
	if (!Nonnegative(inA) || !std::isfinite(inB) || inB <= 0) return std::numeric_limits<double>::quiet_NaN();
	if (inA == 0) return 0;
	double value = inA / inB;
	return std::isnormal(value)? value : std::numeric_limits<double>::quiet_NaN();
}
inline double Charge(double inImpulse, const Thruster &inThruster, double inRate)
{
	if (!Nonnegative(inImpulse)) return std::numeric_limits<double>::quiet_NaN();
	if (inImpulse == 0 || inRate == 0) return 0;
	return Product(Quotient(inImpulse, inThruster.mMaxForce), inRate);
}
inline double PowerScale(double inBudget, double inPower)
{
	if (inBudget == 0) return 0;
	double scale = std::nextafter(Quotient(inBudget, inPower), 0.0);
	return std::isnormal(scale)? scale : std::numeric_limits<double>::quiet_NaN();
}
inline bool Finite(Vec4Arg inValue)
{
	for (uint i = 0; i < 4; ++i) if (!std::isfinite(inValue[i])) return false;
	return true;
}
inline bool Representable(const Wrench &inValue)
{
	for (double x : inValue)
		if (!std::isfinite(x) || std::abs(x) > std::numeric_limits<float>::max()) return false;
	return true;
}
inline bool Valid(const Thruster &inThruster)
{
	const auto &t = inThruster;
	if (!Finite(t.mMount) || !Finite(t.mDirection)) return false;
	double norm = 0;
	for (uint i = 0; i < 4; ++i) norm += double(t.mDirection[i]) * double(t.mDirection[i]);
	if (std::abs(norm - 1) > 1.0e-5) return false;
	for (double x : { t.mMaxForce, t.mResponseTime, t.mPropellantRate, t.mPowerDraw, t.mHeatRate, t.mHealth, t.mForce, t.mHeat })
		if (!Physical(x)) return false;
	return t.mHealth <= 1;
}
inline double Capacity(const Thruster &inThruster)
{
	return inThruster.mEnabled && !inThruster.mFailed? inThruster.mMaxForce * inThruster.mHealth : 0;
}

/// Analytic solution of f' = (target - f) / tau, with cancellation-free small-dt mean.
inline double IntegratedForce(const Thruster &inThruster, double inTarget, double inTime)
{
	if (inTime == 0) return 0;
	if (inThruster.mResponseTime == 0) return Product(inTarget, inTime);
	if (inThruster.mForce == 0 && inTarget == 0) return 0;
	double x = Quotient(inTime, inThruster.mResponseTime);
	double h = Quotient(-std::expm1(-x), x);
	double v = x < 1.0e-3? Product(x, 0.5 + x * (-1.0 / 6 + x * (1.0 / 24 + x * (-1.0 / 120 + x / 720)))) : 1 - h;
	double mean = Product(inThruster.mForce, h) + Product(inTarget, v);
	return Product(inTime, std::clamp(mean, std::min(inThruster.mForce, inTarget), std::max(inThruster.mForce, inTarget)));
}

inline double FinalForce(const Thruster &inThruster, double inTarget, double inTime)
{
	if (inThruster.mResponseTime == 0) return inTarget;
	if (inTime == 0 || (inThruster.mForce == 0 && inTarget == 0)) return inThruster.mForce;
	double x = Quotient(inTime, inThruster.mResponseTime);
	double decay = std::exp(-x);
	if (inThruster.mForce > 0 && !std::isnormal(decay)) return std::numeric_limits<double>::quiet_NaN();
	double force = Product(inThruster.mForce, decay) + Product(inTarget, -std::expm1(-x));
	return std::clamp(force, std::min(inThruster.mForce, inTarget), std::max(inThruster.mForce, inTarget));
}

inline Wrench Pack(Vec4Arg inForce, const Bivec &inTorque)
{
	return { inForce[0], inForce[1], inForce[2], inForce[3], inTorque[0], inTorque[1], inTorque[2], inTorque[3], inTorque[4], inTorque[5] };
}
inline Vec4 Force(const Wrench &inWrench)
{
	return Vec4(float(inWrench[0]), float(inWrench[1]), float(inWrench[2]), float(inWrench[3]));
}
inline Bivec Torque(const Wrench &inWrench)
{
	return Bivec(float(inWrench[4]), float(inWrench[5]), float(inWrench[6]), float(inWrench[7]), float(inWrench[8]), float(inWrench[9]));
}
}

/// A pure, bounded ship-system step. Minimize sum_j weight_j * (A * thrust - request)_j^2
/// with 0 <= thrust_i <= maximum_i * health_i, using 256 cyclic coordinate-descent sweeps
/// from zero. This is an approximate box least-squares solve, NOT a convergence guarantee or
/// a fuel-optimal allocator. Zero weights ignore axes. Unique uint64 IDs define all ordering.
/// At most 4096 thrusters; no hidden gravity or physics mutation. Saturation is observable in
/// per-actuator target/capacity; allocation residuals exclude response lag and later depletion,
/// whereas impulse residuals report the actual delivered interval error on all ten axes.
///
/// Model: each actuator follows df/dt = (target - f) / response_time (zero means instantaneous),
/// analytically integrated while its fixed boundary command remains available. A small-dt
/// series evaluates 1 - (1 - exp(-x)) / x through x^5 (absolute remainder <= x^6/5040 for
/// 0 <= x < 1e-3). Rates and heat are linear in delivered force / nominal maximum, including
/// spool-down. No idle consumption, cooling, recharge, automatic thermal failure or min throttle.
/// Shared instantaneous power reserves a conservative max(initial, target) envelope, clipping
/// excessive spool-down state at the boundary when necessary. Failed, disabled, derated and
/// starved actuators likewise hard-cut to safe capacity; this is not stored mechanical energy.
/// An outward-rounded post-scaling bound certifies power <= limit without an excess tolerance;
/// at most eight additional scale corrections are attempted before atomic numeric rejection.
/// Shared fuel/energy depletion splits the interval at the first event (at most two events),
/// using a feasible one-ULP double-time bracket in <= 64 evaluations per resource per phase.
/// Independent actuators continue, but commands are not reallocated until the next caller step.
/// No consumption is charged for the unflown remainder or for unused rounding-scale resources.
/// Depleted flags identify in-call cutoff events; a tiny unconsumed resource residue may remain.
///
/// Geometry, COM and requested wrench are held in ONE fixed ship-local frame over the interval;
/// torque is direction wedge (mount - COM), using the repository's F wedge r convention.
/// The physics owner handles orientation/time sampling and any COM changes by splitting calls.
/// Returned Vec4/Bivec impulses are floats; per-actuator impulses and resource ledgers are doubles.
/// mPropellantConsumed is a RESOURCE mass-flux ledger, not a body mass mutation. Per-thruster
/// flux is deliveredImpulse / maximum * propellantRate. The parent mass adapter must explicitly
/// model fuel location, exhaust momentum/energy and changing inertia; subtracting this ledger
/// from mass and applying an impulse is NOT by itself a variable-mass rocket integrator.
///
/// Inputs are never changed. Invalid/nonfinite inputs, duplicate IDs, negative quantities,
/// nonunit directions (squared-norm tolerance 1e-5), excess thruster count, numeric overflow or
/// unrepresentable float outputs return false with outResult unchanged. Allocation failures
/// follow the host STL/Jolt out-of-memory policy (not converted to false in no-exception builds).
/// Numeric envelope: each double-valued Thruster field, request weight and delta time is zero
/// or in [1e-100, 1e100] (health is also <= 1). Resource budgets/power limits remain arbitrary
/// finite nonnegative doubles, including the default unlimited power. This is NOT an all-exponent
/// binary64 solver; unsupported intermediate resource/response ranges also reject atomically.
/// In particular, nonzero products/quotients in resource charges, power scaling and analytic
/// response must be normal finite doubles; subnormal, underflowed or overflowed intermediates
/// (including a nonzero state's exponential decay factor) reject. Depletion search is confined
/// to durations >= 1e-120 seconds; an earlier event rejects, rather than granting free impulse.
/// Same-build ordering is deterministic; std::exp/expm1 do NOT promise cross-platform bit replay.
inline bool StepPropulsion(const PropulsionRequest &inRequest, const std::vector<Thruster> &inThrusters,
	Vec4Arg inCenterOfMass, const PropulsionResources &inResources, double inDeltaTime, PropulsionResult &outResult)
{
	using namespace PropulsionDetail;
	if (inThrusters.size() > 4096 || !Finite(inCenterOfMass) || !Representable(Pack(inRequest.mForce, inRequest.mTorque))
		|| !Physical(inDeltaTime) || !Nonnegative(inResources.mPropellant) || !Nonnegative(inResources.mEnergy) || !Nonnegative(inResources.mPowerLimit)) return false;
	for (double weight : inRequest.mWeights) if (!Physical(weight)) return false;
	for (const Thruster &t : inThrusters) if (!Valid(t)) return false;
	PropulsionResult result;
	result.mResources = inResources;
	for (const Thruster &t : inThrusters)
		result.mActuators.push_back({ t, 0, 0 });
	std::sort(result.mActuators.begin(), result.mActuators.end(), [](const auto &a, const auto &b) { return a.mThruster.mID < b.mThruster.mID; });
	for (size_t i = 1; i < result.mActuators.size(); ++i)
		if (result.mActuators[i - 1].mThruster.mID == result.mActuators[i].mThruster.mID) return false;
	std::vector<Wrench> columns;
	std::vector<double> capacities;
	for (auto &a : result.mActuators)
	{
		const Thruster &t = a.mThruster;
		columns.push_back(Pack(t.mDirection, Bivec::sWedge(t.mDirection, t.mMount - inCenterOfMass)));
		double capacity = Capacity(t);
		if ((inResources.mPropellant == 0 && t.mPropellantRate > 0) || (inResources.mEnergy == 0 && t.mPowerDraw > 0)) capacity = 0;
		capacities.push_back(capacity);
		a.mThruster.mForce = std::min(t.mForce, capacity);
	}
	Wrench residual = Pack(inRequest.mForce, inRequest.mTorque);
	for (int sweep = 0; sweep < 256; ++sweep)
		for (size_t i = 0; i < columns.size(); ++i)
		{
			double gradient = 0, diagonal = 0;
			for (size_t j = 0; j < 10; ++j)
			{
				gradient += columns[i][j] * inRequest.mWeights[j] * residual[j];
				diagonal += columns[i][j] * inRequest.mWeights[j] * columns[i][j];
			}
			if (!std::isfinite(gradient) || !std::isfinite(diagonal)) return false;
			if (diagonal == 0) continue;
			auto &a = result.mActuators[i];
			double next = std::clamp(a.mTargetForce + gradient / diagonal, 0.0, capacities[i]);
			double change = next - a.mTargetForce;
			a.mTargetForce = next;
			for (size_t j = 0; j < 10; ++j) residual[j] -= change * columns[i][j];
		}
	// Reserve power for the steady command first. Preserve existing spool state unless the
	// envelope sum(max(initial, target) * specific power) would exceed the shared limit.
	// Then hard-cut ONLY the excess above targets. Unlike scaling initial state on every tick,
	// this preserves the exact exponential response under a constant power-limited request.
	double targetPower = 0;
	for (const auto &a : result.mActuators)
		if (a.mThruster.mMaxForce > 0)
			targetPower += Charge(a.mTargetForce, a.mThruster, a.mThruster.mPowerDraw);
	if (!std::isfinite(targetPower)) return false;
	if (targetPower > inResources.mPowerLimit)
	{
		result.mPowerLimited = true;
		double scale = PowerScale(inResources.mPowerLimit, targetPower);
		for (auto &a : result.mActuators)
			if (a.mThruster.mPowerDraw > 0) a.mTargetForce = Product(a.mTargetForce, scale);
	}
	targetPower = 0;
	double excessPower = 0;
	for (const auto &a : result.mActuators)
		if (a.mThruster.mMaxForce > 0)
		{
			targetPower += Charge(a.mTargetForce, a.mThruster, a.mThruster.mPowerDraw);
			excessPower += Charge(std::max(0.0, a.mThruster.mForce - a.mTargetForce), a.mThruster, a.mThruster.mPowerDraw);
		}
	if (!std::isfinite(targetPower) || !std::isfinite(excessPower)) return false;
	double headroom = std::max(0.0, inResources.mPowerLimit - targetPower);
	if (excessPower > headroom)
	{
		result.mPowerLimited = true;
		double scale = PowerScale(headroom, excessPower);
		for (auto &a : result.mActuators)
			if (a.mThruster.mPowerDraw > 0 && a.mThruster.mForce > a.mTargetForce)
				a.mThruster.mForce = a.mTargetForce + Product(a.mThruster.mForce - a.mTargetForce, scale);
	}
	// One downward ULP on the scale does not bound the rounded products and a many-term
	// sum. Certify the actual post-scaling envelope with outward-rounded arithmetic.
	// Each nextafter bounds one correctly rounded nonnegative operation from above;
	// this remains conservative when long double has no extra precision (e.g. MSVC).
	auto powerEnvelope = [&]()
	{
		double upper = 0;
		const double infinity = std::numeric_limits<double>::infinity();
		for (const auto &a : result.mActuators)
		{
			if (!Nonnegative(a.mTargetForce) || !Nonnegative(a.mThruster.mForce)) return std::numeric_limits<double>::quiet_NaN();
			double force = std::max(a.mTargetForce, a.mThruster.mForce);
			if (force == 0 || a.mThruster.mPowerDraw == 0) continue;
			double fraction = std::nextafter(Quotient(force, a.mThruster.mMaxForce), infinity);
			double power = std::nextafter(Product(fraction, a.mThruster.mPowerDraw), infinity);
			upper = std::nextafter(upper + power, infinity);
		}
		return upper;
	};
	for (int correction = 0; ; ++correction)
	{
		double envelope = powerEnvelope();
		if (!Nonnegative(envelope)) return false;
		if (envelope <= inResources.mPowerLimit) break;
		if (correction == 8) return false;
		result.mPowerLimited = true;
		double scale = PowerScale(inResources.mPowerLimit, envelope);
		for (auto &a : result.mActuators)
			if (a.mThruster.mPowerDraw > 0)
			{
				a.mTargetForce = Product(a.mTargetForce, scale);
				a.mThruster.mForce = Product(a.mThruster.mForce, scale);
			}
	}
	Wrench desired = Pack(inRequest.mForce, inRequest.mTorque);
	residual = desired;
	for (size_t i = 0; i < columns.size(); ++i)
		for (size_t j = 0; j < 10; ++j) residual[j] -= columns[i][j] * result.mActuators[i].mTargetForce;
	if (!Representable(residual)) return false;
	result.mAllocationResidualForce = Force(residual);
	result.mAllocationResidualTorque = Torque(residual);

	std::vector<double> targets;
	for (const auto &a : result.mActuators) targets.push_back(a.mTargetForce);
	auto usage = [&](double time)
	{
		std::array<double, 2> cost {};
		for (size_t i = 0; i < targets.size(); ++i)
		{
			const auto &a = result.mActuators[i];
			if (a.mThruster.mMaxForce == 0) continue;
			double impulse = IntegratedForce(a.mThruster, targets[i], time);
			cost[0] += Charge(impulse, a.mThruster, a.mThruster.mPropellantRate);
			cost[1] += Charge(impulse, a.mThruster, a.mThruster.mPowerDraw);
		}
		return cost;
	};
	bool fuelBlocked = inResources.mPropellant == 0, energyBlocked = inResources.mEnergy == 0;
	double remaining = inDeltaTime;
	// At most two resource exhaustion events, followed by the resource-independent remainder.
	for (int phase = 0; phase < 3; ++phase)
	{
		for (size_t i = 0; i < targets.size(); ++i)
			if ((fuelBlocked && result.mActuators[i].mThruster.mPropellantRate > 0)
				|| (energyBlocked && result.mActuators[i].mThruster.mPowerDraw > 0))
			{
				targets[i] = 0;
				result.mActuators[i].mThruster.mForce = 0;
			}
		auto fullCost = usage(remaining);
		if (!Nonnegative(fullCost[0]) || !Nonnegative(fullCost[1])) return false;
		std::array<double, 2> budget { result.mResources.mPropellant, result.mResources.mEnergy };
		std::array<double, 2> eventTime { remaining, remaining };
		std::array<bool, 2> exhaust { false, false };
		for (size_t r = 0; r < 2; ++r)
			if (fullCost[r] > 0 && fullCost[r] >= budget[r])
			{
				exhaust[r] = true;
				// IEEE positive-double bit patterns are ordered. Bisect representable times,
				// not a fixed absolute time grid: supported tiny budgets get their impulse. At most
				// 64 evaluations leave adjacent doubles (one-ULP event bracket), low feasible.
				static_assert(sizeof(double) == sizeof(uint64) && std::numeric_limits<double>::is_iec559, "IEEE binary64 required");
				double minimumTime = std::min(remaining, 1.0e-120);
				auto minimumCost = usage(minimumTime);
				if (!Nonnegative(minimumCost[0]) || !Nonnegative(minimumCost[1]) || minimumCost[r] > budget[r]) return false;
				uint64 low, high;
				std::memcpy(&low, &minimumTime, sizeof(low));
				std::memcpy(&high, &remaining, sizeof(high));
				if (fullCost[r] == budget[r]) low = high;
				for (int iteration = 0; iteration < 64 && high - low > 1; ++iteration)
				{
					uint64 key = low + (high - low) / 2;
					double mid;
					std::memcpy(&mid, &key, sizeof(mid));
					auto middleCost = usage(mid);
					if (!Nonnegative(middleCost[0]) || !Nonnegative(middleCost[1])) return false;
					if (middleCost[r] <= budget[r]) low = key;
					else high = key;
				}
				std::memcpy(&eventTime[r], &low, sizeof(double));
			}
		double duration = std::min(eventTime[0], eventTime[1]);
		auto cost = usage(duration);
		if (!Nonnegative(cost[0]) || !Nonnegative(cost[1]) || cost[0] > budget[0] || cost[1] > budget[1]) return false;
		result.mResources.mPropellant -= cost[0];
		result.mResources.mEnergy -= cost[1];
		result.mPropellantConsumed += cost[0];
		result.mEnergyConsumed += cost[1];
		for (size_t i = 0; i < targets.size(); ++i)
		{
			auto &a = result.mActuators[i];
			double delivered = IntegratedForce(a.mThruster, targets[i], duration);
			a.mDeliveredImpulse += delivered;
			a.mThruster.mForce = FinalForce(a.mThruster, targets[i], duration);
			double heat = Charge(delivered, a.mThruster, a.mThruster.mHeatRate);
			// Validate before a subsequent resource hard-cut can overwrite a failed state.
			if (!Nonnegative(delivered) || !Nonnegative(a.mThruster.mForce) || !Nonnegative(heat)) return false;
			a.mThruster.mHeat += heat;
			result.mHeatGenerated += heat;
		}
		if (exhaust[0] && eventTime[0] <= eventTime[1]) fuelBlocked = result.mPropellantDepleted = true;
		if (exhaust[1] && eventTime[1] <= eventTime[0]) energyBlocked = result.mEnergyDepleted = true;
		remaining -= duration;
		if (!exhaust[0] && !exhaust[1]) break;
	}
	// Exhaustion is a hard shutdown even when it occurs precisely at the interval boundary.
	for (auto &a : result.mActuators)
		if ((fuelBlocked && a.mThruster.mPropellantRate > 0) || (energyBlocked && a.mThruster.mPowerDraw > 0)) a.mThruster.mForce = 0;
	Wrench impulse {}, missing {};
	for (size_t i = 0; i < columns.size(); ++i)
		for (size_t j = 0; j < 10; ++j) impulse[j] += columns[i][j] * result.mActuators[i].mDeliveredImpulse;
	for (size_t j = 0; j < 10; ++j) missing[j] = desired[j] * inDeltaTime - impulse[j];
	if (!Representable(impulse) || !Representable(missing) || !Nonnegative(result.mHeatGenerated)
		|| !Nonnegative(result.mPropellantConsumed) || !Nonnegative(result.mEnergyConsumed)) return false;
	for (const auto &a : result.mActuators)
		if (!Nonnegative(a.mDeliveredImpulse) || !Nonnegative(a.mThruster.mForce) || !Nonnegative(a.mThruster.mHeat)) return false;
	result.mImpulse = Force(impulse);
	result.mAngularImpulse = Torque(impulse);
	result.mResidualImpulse = Force(missing);
	result.mResidualAngularImpulse = Torque(missing);
	outResult = std::move(result);
	return true;
}

} // namespace SpaceSim
JPH_NAMESPACE_END
