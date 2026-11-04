// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include "Math/Cl4Reference.h"
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <array>
#include <cstdio>
#include <random>

// Independent double-precision SO(4) reference: no production inertia adjoints,
// exponential, commutator, orientation predictor, or time integrator.
namespace SpaceDynamicsReference
{
	constexpr int cI[6] = { 0, 0, 0, 1, 1, 2 }, cJ[6] = { 1, 2, 3, 2, 3, 3 };
	using Six = std::array<double, 6>;
	struct RefMatrix
	{
		double a[4][4] = {};
		RefMatrix operator + (const RefMatrix &b) const { RefMatrix r; for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) r.a[i][j] = a[i][j] + b.a[i][j]; return r; }
		RefMatrix operator * (double b) const { RefMatrix r; for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) r.a[i][j] = a[i][j] * b; return r; }
		RefMatrix operator * (const RefMatrix &b) const { RefMatrix r; for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) for (int k = 0; k < 4; ++k) r.a[i][j] += a[i][k] * b.a[k][j]; return r; }
		RefMatrix Transposed() const { RefMatrix r; for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) r.a[i][j] = a[j][i]; return r; }
	};
	static RefMatrix Rotation(RotorArg r)
	{
		RefMatrix result;
		auto q = Cl4Ref::FromRotor(r);
		for (int j = 0; j < 4; ++j)
		{
			Cl4Ref::MV axis; axis.mC[1 << j] = 1;
			auto v = Cl4Ref::Mul(Cl4Ref::Mul(q, axis), Cl4Ref::Reverse(q));
			for (int i = 0; i < 4; ++i) result.a[i][j] = v.mC[1 << i];
		}
		return result;
	}
	static Six Components(BivecArg b) { Six r; for (int k = 0; k < 6; ++k) r[k] = b[k]; return r; }
	static RefMatrix Skew(const Six &b) { RefMatrix r; for (int k = 0; k < 6; ++k) { r.a[cI[k]][cJ[k]] = b[k]; r.a[cJ[k]][cI[k]] = -b[k]; } return r; }
	static Six Components(const RefMatrix &b) { Six r; for (int k = 0; k < 6; ++k) r[k] = 0.5 * (b.a[cI[k]][cJ[k]] - b.a[cJ[k]][cI[k]]); return r; }
	static double Norm(const Six &a) { double r = 0; for (double v : a) r += v * v; return sqrt(r); }
	static double Error(const Six &a, const Six &b) { Six d; for (int k = 0; k < 6; ++k) d[k] = a[k] - b[k]; return Norm(d); }
	static double Error(const RefMatrix &a, const RefMatrix &b) { double r = 0; for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) r += Square(a.a[i][j] - b.a[i][j]); return sqrt(r); }
	static Six Inertia(const Body &b)
	{
		Six d; for (int k = 0; k < 6; ++k) d[k] = 1.0 / b.GetMotionProperties()->GetInverseInertiaDiagonal()[k]; return d;
	}
	static RefMatrix PrincipalRotation(const Body &b) { return Rotation(b.GetRotation()) * Rotation(b.GetMotionProperties()->GetInertiaRotation()); }
	static Six Momentum(const Body &b)
	{
		RefMatrix q = PrincipalRotation(b);
		Six l = Components(q.Transposed() * Skew(Components(b.GetAngularVelocity())) * q), d = Inertia(b);
		for (int k = 0; k < 6; ++k) l[k] *= d[k];
		return Components(q * Skew(l) * q.Transposed());
	}
	static double Energy(const Body &b) { Six w = Components(b.GetAngularVelocity()), l = Momentum(b); double r = 0; for (int k = 0; k < 6; ++k) r += 0.5 * w[k] * l[k]; return r; }
	static Six Velocity(const RefMatrix &q, const Six &world_l, const Six &d)
	{
		Six w = Components(q.Transposed() * Skew(world_l) * q);
		for (int k = 0; k < 6; ++k) w[k] /= d[k];
		return Components(q * Skew(w) * q.Transposed());
	}
	static RefMatrix Derivative(const RefMatrix &q, const Six &l, const Six &d) { return Skew(Velocity(q, l, d)) * q; }
	static RefMatrix Reference(RefMatrix q, const Six &l, const Six &d, double h, int n)
	{
		for (int step = 0; step < n; ++step)
		{
			RefMatrix k1 = Derivative(q, l, d), k2 = Derivative(q + k1 * (0.5 * h), l, d);
			RefMatrix k3 = Derivative(q + k2 * (0.5 * h), l, d), k4 = Derivative(q + k3 * h, l, d);
			q = q + (k1 + k2 * 2 + k3 * 2 + k4) * (h / 6);
		}
		return q;
	}
	static BodyCreationSettings Settings()
	{
		BodyCreationSettings s(new BoxShape(Vec4(0.35f, 0.65f, 1.05f, 1.6f), 0.0f), RVec4::sZero(),
			Rotor::sExp(Bivec(0.2f, -0.1f, 0.15f, 0.12f, -0.08f, 0.18f)), EMotionType::Dynamic, Layers::MOVING);
		s.mApplyGyroscopicForce = true;
		s.mLinearDamping = s.mAngularDamping = 0;
		s.mAllowSleeping = false;
		s.mGravityFactor = 0;
		s.mMaxLinearVelocity = s.mMaxAngularVelocity = 1.0e6f;
		s.mAngularVelocity = Bivec(1.1f, -0.7f, 0.9f, 0.6f, -1.2f, 0.8f);
		return s;
	}
}

TEST_SUITE("SpaceDynamicsTests")
{
	using namespace SpaceDynamicsReference;
	TEST_CASE("PrincipalSecondMomentsProduceSO4FrameThroughRotorRoundTrip")
	{
		std::mt19937 random(0x4d01);
		std::uniform_real_distribution<double> angle(-1.7, 1.7);
		for (int sample = 0; sample < 65; ++sample)
		{
			// Seeded general rotations, including a known diagonal regression.
			Cl4Ref::MV generator;
			for (int k = 0; k < 6; ++k)
				generator.mC[Cl4Ref::cBivecMask[k]] = sample == 0? 0 : angle(random);
			RefMatrix q = Rotation(Cl4Ref::ToRotor(Cl4Ref::Exp(generator)));
			RefMatrix diagonal; for (int i = 0; i < 4; ++i) diagonal.a[i][i] = double(1 << i);
			RefMatrix tensor = q * diagonal * q.Transposed();
			MassProperties mass; mass.mMass = 1; mass.mInertia = Mat44::sZero();
			for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) mass.mInertia(i, j) = float(tensor.a[i][j]);
			Mat44 frame; Vec4 moments;
			CHECK(mass.DecomposePrincipalMomentsOfInertia(frame, moments));
			// Independent Leibniz determinant, not the cross-product sign under test.
			double determinant = 0;
			int permutation[4] = { 0, 1, 2, 3 };
			do
			{
				int inversions = 0; double term = 1;
				for (int i = 0; i < 4; ++i)
				{
					term *= frame(i, permutation[i]);
					for (int j = i + 1; j < 4; ++j) inversions += permutation[i] > permutation[j];
				}
				determinant += (inversions & 1)? -term : term;
			} while (std::next_permutation(permutation, permutation + 4));
			CAPTURE(sample); CAPTURE(determinant);
			CHECK(abs(determinant - 1) < 2.0e-5);
			Rotor r = frame.GetRotor();
			CHECK(r.IsNormalized(2.0e-5f));
			RefMatrix rotor_frame = Rotation(r), principal;
			for (int i = 0; i < 4; ++i) principal.a[i][i] = moments[i];
			CHECK(Error(rotor_frame * principal * rotor_frame.Transposed(), tensor) < 5.0e-5);
		}
	}

	TEST_CASE("GyroscopicAdjointMatchesIndependentCliffordSandwich")
	{
		std::mt19937 random(0x4d02);
		std::uniform_real_distribution<double> angle(-1.7, 1.7);
		for (int sample = 0; sample < 64; ++sample)
		{
			Cl4Ref::MV generator;
			for (int k = 0; k < 6; ++k) generator.mC[Cl4Ref::cBivecMask[k]] = angle(random);
			Rotor r = Cl4Ref::ToRotor(Cl4Ref::Exp(generator));
			RefMatrix q = Rotation(r);
			// Normalize the vector sandwich by its measured basis-column length.
			double length_sq = 0;
			for (int i = 0; i < 4; ++i) length_sq += Square(q.a[i][0]);
			q = q * (1.0 / sqrt(length_sq));
			Six local, actual, recovered;
			for (int k = 0; k < 6; ++k) local[k] = angle(random);
			MotionPropertiesImpl::GyroscopicAdjoint adjoint(r);
			adjoint.ToWorld(local.data(), actual.data());
			adjoint.ToLocal(actual.data(), recovered.data());
			CHECK(Error(actual, Components(q * Skew(local) * q.Transposed())) < 1.0e-13);
			CHECK(Error(recovered, local) < 1.0e-13);
		}
	}

	TEST_CASE("ZeroTimeGyroscopicTransportDoesNotAmplifyRotorRoundoff")
	{
		PhysicsTestContext c; c.ZeroGravity();
		Body &b = c.CreateBody(Settings(), EActivation::Activate);
		Six w = Components(b.GetAngularVelocity());
		for (int i = 0; i < 1000; ++i) b.IntegrateRotationInternal(0);
		CHECK(Error(Components(b.GetAngularVelocity()), w) < 1.0e-5);
	}

	TEST_CASE("AsymmetricFreeSpinConservesWorldMomentumAndBoundsEnergy")
	{
		for (int hz : { 60, 120, 150 })
		{
			PhysicsTestContext c(1.0f / hz); c.ZeroGravity();
			Body &b = c.CreateBody(Settings(), EActivation::Activate);
			Six l = Momentum(b); double e = Energy(b), max_l = 0, max_e = 0;
			const int seconds = hz == 120? 600 : 60;
			for (int n = 0; n < seconds * hz; ++n)
			{
				CHECK(c.SimulateSingleStep() == EPhysicsUpdateError::None);
				max_l = std::max(max_l, Error(Momentum(b), l) / Norm(l));
				max_e = std::max(max_e, abs(Energy(b) / e - 1));
				CHECK(b.GetRotation().IsNormalized(1.0e-5f));
			}
			std::printf("space free hz=%d seconds=%d max_rel_L=%.9g max_rel_E=%.9g max_abs_L=%.9g max_abs_E=%.9g\n", hz, seconds, max_l, max_e, max_l * Norm(l), max_e * e);
			CAPTURE(hz); CAPTURE(max_l); CAPTURE(max_e);
			CHECK(max_l < 1.0e-4);
			CHECK(max_e < 2.0e-4);
			CHECK(b.IsActive());
			CHECK(b.GetAngularVelocity().Length() < 10);
		}
	}

	TEST_CASE("OffCenterForceAndImpulseHaveOneOwnerInAllSixPlanes")
	{
		for (int plane = 0; plane < 6; ++plane)
			for (bool force : { false, true })
			{
				PhysicsTestContext c(1.0f / 120); c.ZeroGravity();
				auto settings = Settings();
				// Keep the injected 1.4-unit angular impulse comparable to initial L;
				// a large default-density mass must not hide a missing torque lane.
				settings.mOverrideMassProperties = EOverrideMassProperties::CalculateInertia;
				settings.mMassPropertiesOverride.mMass = 1;
				Body &b = c.CreateBody(settings, EActivation::Activate);
				Six expected = Momentum(b); const double scale = Norm(expected);
				Vec4 impulse = Vec4::sZero(), offset = Vec4::sZero();
				impulse[cI[plane]] = 2; offset[cJ[plane]] = 0.7f;
				if (force) b.AddForce(impulse / c.GetDeltaTime(), b.GetCenterOfMassPosition() + RVec4(offset));
				else b.AddImpulse(impulse, b.GetCenterOfMassPosition() + RVec4(offset));
				expected[plane] += 1.4;
				for (int step = 0; step < 120; ++step)
				{
					CHECK(c.SimulateSingleStep() == EPhysicsUpdateError::None);
					CAPTURE(plane); CAPTURE(force); CAPTURE(step);
					CHECK(Error(Momentum(b), expected) / scale < 5.0e-5);
					CHECK((b.GetLinearVelocity() / b.GetMotionProperties()->GetInverseMass() - impulse).Length() < 1.0e-5f);
				}
			}
	}

	TEST_CASE("GeneralRotatedPhysicalInertiaMatchesIndependentFreeSpin")
	{
		std::mt19937 random(0x4d03);
		std::uniform_real_distribution<double> angle(-1.7, 1.7);
		for (int sample = 0; sample < 8; ++sample)
		{
			Cl4Ref::MV generator;
			for (int k = 0; k < 6; ++k) generator.mC[Cl4Ref::cBivecMask[k]] = angle(random);
			RefMatrix frame = Rotation(Cl4Ref::ToRotor(Cl4Ref::Exp(generator))), diagonal;
			for (int i = 0; i < 4; ++i) diagonal.a[i][i] = double(1 << i);
			RefMatrix tensor = frame * diagonal * frame.Transposed();
			auto settings = Settings(); settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
			settings.mMassPropertiesOverride.mMass = 1;
			for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) settings.mMassPropertiesOverride.mInertia(i, j) = float(tensor.a[i][j]);
			PhysicsTestContext c(1.0f / 120); c.ZeroGravity();
			Body &b = c.CreateBody(settings, EActivation::Activate);
			Six l = Momentum(b), d = Inertia(b); double e = Energy(b), max_l = 0, max_e = 0;
			RefMatrix reference = Reference(PrincipalRotation(b), l, d, 1.0 / 3840, 2 * 3840);
			for (int n = 0; n < 240; ++n)
			{
				CHECK(c.SimulateSingleStep() == EPhysicsUpdateError::None);
				max_l = std::max(max_l, Error(Momentum(b), l) / Norm(l));
				max_e = std::max(max_e, abs(Energy(b) / e - 1));
			}
			CAPTURE(sample); CAPTURE(max_l); CAPTURE(max_e);
			CHECK(max_l < 1.0e-5);
			CHECK(max_e < 2.0e-4);
			CHECK(Error(PrincipalRotation(b), reference) < 0.001);
		}
	}

	TEST_CASE("GyroscopicContactImpulsesConserveTotalMomentum")
	{
		PhysicsTestContext c(1.0f / 120); c.ZeroGravity();
		auto physics = c.GetSystem()->GetPhysicsSettings();
		physics.mNumPositionSteps = 0; // Separate velocity impulses from nonphysical position projection.
		c.GetSystem()->SetPhysicsSettings(physics);
		auto settings = Settings(); settings.mRotation = Rotor::sIdentity();
		settings.mFriction = 0; settings.mRestitution = 0.5f;
		settings.mAngularVelocity *= 0.15f;
		settings.mPosition = RVec4(-2, 0, 0.1f, 0.15f); settings.mLinearVelocity = Vec4(3, 0, 0, 0);
		Body &a = c.CreateBody(settings, EActivation::Activate);
		settings.mPosition = RVec4(2, 0, -0.1f, -0.15f); settings.mLinearVelocity = Vec4(-3, 0, 0, 0);
		Body &b = c.CreateBody(settings, EActivation::Activate);
		auto origin_momentum = [](const Body &body)
		{
			Six l = Momentum(body);
			Vec4 p = body.GetLinearVelocity() / body.GetMotionProperties()->GetInverseMass();
			Vec4 x(body.GetCenterOfMassPosition());
			for (int k = 0; k < 6; ++k) l[k] += double(p[cI[k]]) * x[cJ[k]] - double(p[cJ[k]]) * x[cI[k]];
			return l;
		};
		auto total_l = [&]() { Six l = origin_momentum(a), other = origin_momentum(b); for (int k = 0; k < 6; ++k) l[k] += other[k]; return l; };
		auto total_p = [&]() { return a.GetLinearVelocity() / a.GetMotionProperties()->GetInverseMass() + b.GetLinearVelocity() / b.GetMotionProperties()->GetInverseMass(); };
		auto total_e = [&]() { return Energy(a) + Energy(b) + 0.5 * (a.GetLinearVelocity().LengthSq() / a.GetMotionProperties()->GetInverseMass() + b.GetLinearVelocity().LengthSq() / b.GetMotionProperties()->GetInverseMass()); };
		Six initial_l = total_l(); Vec4 initial_p = total_p(); double initial_e = total_e();
		for (int n = 0; n < 240; ++n) CHECK(c.SimulateSingleStep() == EPhysicsUpdateError::None);
		CHECK(total_e() < 0.99 * initial_e); // Prove they actually collided.
		CHECK(Error(total_l(), initial_l) / Norm(initial_l) < 1.0e-4);
		CHECK((total_p() - initial_p).Length() * a.GetMotionProperties()->GetInverseMass() < 1.0e-4f);
	}

	TEST_CASE("CommutingDoubleRotationHasIndependentAnalyticalSolution")
	{
		PhysicsTestContext c(1.0f / 120); c.ZeroGravity();
		auto settings = Settings(); settings.mRotation = Rotor::sIdentity();
		Body &b = c.CreateBody(settings, EActivation::Activate);
		// In the actual principal frame two disjoint planes commute even for
		// unequal moments. Compare via Clifford power series, not PredictOrientation.
		Rotor frame = b.GetMotionProperties()->GetInertiaRotation();
		Bivec local_w(1.2f, 0, 0, 0, 0, -0.8f);
		b.SetAngularVelocity(frame.RotateBivec(local_w));
		Six w0 = Components(b.GetAngularVelocity());
		auto exp = Cl4Ref::Exp(Cl4Ref::FromBivec(b.GetAngularVelocity()) * 1.0);
		RefMatrix q0 = Rotation(b.GetRotation());
		for (int n = 0; n < 240; ++n) CHECK(c.SimulateSingleStep() == EPhysicsUpdateError::None);
		CHECK(Error(Rotation(b.GetRotation()), Rotation(Cl4Ref::ToRotor(exp)) * q0) < 1.0e-4);
		CHECK(Error(Components(b.GetAngularVelocity()), w0) < 1.0e-4);
	}

	TEST_CASE("AsymmetricFreeSpinRefinesAgainstIndependentRK4AtSameDuration")
	{
		double errors[3] = {};
		for (int level = 0; level < 3; ++level)
		{
			int hz = 30 << level;
			PhysicsTestContext c(1.0f / hz); c.ZeroGravity();
			Body &b = c.CreateBody(Settings(), EActivation::Activate);
			RefMatrix q0 = PrincipalRotation(b); Six l = Momentum(b), d = Inertia(b);
			RefMatrix reference = Reference(q0, l, d, 1.0 / 3840, 2 * 3840);
			RefMatrix finer = Reference(q0, l, d, 1.0 / 7680, 2 * 7680);
			CHECK(Error(reference, finer) < 1.0e-10);
			for (int n = 0; n < 2 * hz; ++n) CHECK(c.SimulateSingleStep() == EPhysicsUpdateError::None);
			errors[level] = Error(PrincipalRotation(b), reference);
			double w_error = Error(Components(b.GetAngularVelocity()), Velocity(reference, l, d));
			std::printf("space refine hz=%d seconds=2 orientation=%.9g omega=%.9g\n", hz, errors[level], w_error);
			CHECK(w_error < 0.005);
		}
		CHECK(errors[0] < 0.01);
		CHECK(errors[1] < 0.4 * errors[0]);
		CHECK(errors[2] < 0.4 * errors[1]);
		CHECK(errors[2] < 0.001);
	}
}
