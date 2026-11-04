// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2026 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include "Layers.h"
#include "Math/Cl4Reference.h"
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>
#include <Jolt/Physics/Collision/Shape/RotatedTranslatedShape.h>
#include <cmath>

TEST_SUITE("MassRoundTripTests")
{
	// Eight point masses, rotated with scalar Givens rotations (not Jolt's
	// matrix/rotor mass helpers). Their outer products are the reference tensor.
	static MassProperties sParticleProperties()
	{
		MassProperties result;
		result.mMass = 12.0f;
		result.mInertia = Mat44::sZero();
		const double radii[] = { 2.0, 3.0, 4.0, 5.0 };
		for (int axis = 0; axis < 4; ++axis)
			for (int sign : { -1, 1 })
			{
				double x[4] = {};
				x[axis] = sign * radii[axis];
				for (int i = 0; i < 4; ++i)
					for (int j = i + 1; j < 4; ++j)
					{
						double angle = 0.17 * (i + 1) - 0.23 * (j + 1);
						double c = std::cos(angle), s = std::sin(angle);
						double a = x[i], b = x[j];
						x[i] = c * a - s * b;
						x[j] = s * a + c * b;
					}
				for (int i = 0; i < 4; ++i)
					for (int j = 0; j < 4; ++j)
						result.mInertia(i, j) += float(1.5 * x[i] * x[j]);
			}
		return result;
	}

	static void sCheckTensor(Mat44Arg inActual, Mat44Arg inExpected)
	{
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
			{
				CAPTURE(i); CAPTURE(j);
				CHECK(abs(inActual(i, j) - inExpected(i, j)) < 2.0e-5f * max(1.0f, abs(inExpected(i, j))));
			}
	}

	// Reference I_ab = integral v_a(x).v_b(x) dm, where v_a is the
	// coordinate-plane unit angular velocity field. Invert its dense 6x6
	// matrix by Gaussian elimination, independently of principal moments.
	static void sCheckAngularResponse(const Body &inBody, Mat44Arg inSecond)
	{
		double rotation[4][4], world_second[4][4] = {};
		for (int col = 0; col < 4; ++col)
		{
			Cl4Ref::MV basis;
			basis.mC[Cl4Ref::cVecMask[col]] = 1.0;
			Cl4Ref::MV rotated = Cl4Ref::Sandwich(Cl4Ref::FromRotor(inBody.GetRotation()), basis);
			for (int row = 0; row < 4; ++row)
				rotation[row][col] = rotated.mC[Cl4Ref::cVecMask[row]];
		}
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				for (int k = 0; k < 4; ++k)
					for (int l = 0; l < 4; ++l)
						world_second[i][j] += rotation[i][k] * inSecond(k, l) * rotation[j][l];
		double generator[6][4][4] = {};
		int plane = 0;
		for (int i = 0; i < 4; ++i)
			for (int j = i + 1; j < 4; ++j, ++plane)
			{
				generator[plane][i][j] = 1.0;
				generator[plane][j][i] = -1.0;
			}
		double inertia[6][6] = {};
		for (int a = 0; a < 6; ++a)
			for (int b = 0; b < 6; ++b)
				for (int row = 0; row < 4; ++row)
					for (int k = 0; k < 4; ++k)
						for (int l = 0; l < 4; ++l)
							inertia[a][b] += generator[a][row][k] * generator[b][row][l] * world_second[k][l];
		const MotionProperties *motion = inBody.GetMotionProperties();
		uint mask = uint(motion->GetAllowedDOFs()) >> 4;
		for (int input = 0; input < 6; ++input)
		{
			double augmented[6][7] = {};
			for (int row = 0; row < 6; ++row)
			{
				for (int col = 0; col < 6; ++col)
					augmented[row][col] = inertia[row][col];
				augmented[row][6] = row == input && (mask & (1U << input)) != 0? 1.0 : 0.0;
			}
			for (int col = 0; col < 6; ++col)
			{
				int pivot = col;
				for (int row = col + 1; row < 6; ++row)
					if (std::abs(augmented[row][col]) > std::abs(augmented[pivot][col]))
						pivot = row;
				for (int k = 0; k < 7; ++k)
					std::swap(augmented[col][k], augmented[pivot][k]);
				CHECK(std::abs(augmented[col][col]) > 1.0e-10);
				if (std::abs(augmented[col][col]) <= 1.0e-10)
					return;
				double factor = augmented[col][col];
				for (int k = col; k < 7; ++k)
					augmented[col][k] /= factor;
				for (int row = 0; row < 6; ++row)
					if (row != col)
					{
						double f = augmented[row][col];
						for (int k = col; k < 7; ++k)
							augmented[row][k] -= f * augmented[col][k];
					}
			}
			float v[6] = {};
			v[input] = 1.0f;
			Bivec impulse(v[0], v[1], v[2], v[3], v[4], v[5]);
			Bivec actual = motion->MultiplyWorldSpaceInverseInertiaByVector(inBody.GetRotation(), impulse);
			for (int output = 0; output < 6; ++output)
			{
				double expected = (mask & (1U << output)) != 0? augmented[output][6] : 0.0;
				CAPTURE(input); CAPTURE(output);
				CHECK(std::abs(actual[output] - expected) < 1.0e-7);
			}
		}
	}

	TEST_CASE("BodyCreationSettingsPreserveFullParticleSecondMoment")
	{
		PhysicsTestContext context;
		// The shape origin deliberately differs from its COM in all four axes.
		BodyCreationSettings settings(new RotatedTranslatedShape(Vec4(0.5f, -1, 2, 3), Rotor::sIdentity(), new SphereShape(0.5f)), RVec4(1, -2, 3, 4),
			Rotor::sExp(Bivec(0.2f, -0.3f, 0.1f, 0.4f, -0.2f, 0.3f)), EMotionType::Dynamic, Layers::MOVING);
		settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
		for (float scale : { 1.0f, 7.3f })
			for (EAllowedDOFs dofs : { EAllowedDOFs::All, EAllowedDOFs::AllRotation, EAllowedDOFs::AllTranslation,
				EAllowedDOFs::TranslationW | EAllowedDOFs::RotationXW | EAllowedDOFs::RotationYW })
			{
				CAPTURE(scale); CAPTURE(uint(dofs));
				settings.mMassPropertiesOverride = sParticleProperties();
				settings.mMassPropertiesOverride.ScaleToMass(12.0f * scale);
				settings.mAllowedDOFs = dofs;
				Body &original = context.CreateBody(settings, EActivation::DontActivate);
				BodyCreationSettings exported = original.GetBodyCreationSettings();
				CHECK(exported.mAllowedDOFs == dofs);
				bool translates = (dofs & EAllowedDOFs::AllTranslation) != EAllowedDOFs::None;
				bool rotates = (dofs & EAllowedDOFs::AllRotation) != EAllowedDOFs::None;
				CHECK(exported.mMassPropertiesOverride.mMass == doctest::Approx(translates? 12.0f * scale : 1.0f));
				sCheckTensor(exported.mMassPropertiesOverride.mInertia, rotates? settings.mMassPropertiesOverride.mInertia : Mat44::sZero());
				Body &recreated = context.CreateBody(exported, EActivation::DontActivate);
				CHECK_APPROX_EQUAL(recreated.GetPosition(), original.GetPosition(), 1.0e-6f);
				CHECK_APPROX_EQUAL(recreated.GetCenterOfMassPosition(), original.GetCenterOfMassPosition(), 1.0e-6f);
				CHECK(recreated.GetMotionProperties()->GetInverseMass() == doctest::Approx(original.GetMotionProperties()->GetInverseMass()));
				sCheckAngularResponse(original, settings.mMassPropertiesOverride.mInertia);
				sCheckAngularResponse(recreated, settings.mMassPropertiesOverride.mInertia);
			}
	}

	TEST_CASE("BodyExportReconstructsKnownPhysicalPrincipalMoments")
	{
		PhysicsTestContext context;
		BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
		Body &body = context.CreateBody(settings, EActivation::DontActivate);
		Rotor principal = Rotor::sExp(Bivec(0.23f, -0.41f, 0.17f, 0.29f, 0.37f, -0.19f));
		const double d[4] = { 12.0, 27.0, 48.0, 75.0 };
		body.GetMotionProperties()->SetInverseMass(1.0f / 12.0f);
		body.GetMotionProperties()->SetInverseInertia(Bivec(1.0f / 39, 1.0f / 60, 1.0f / 87, 1.0f / 75, 1.0f / 102, 1.0f / 123), principal);
		// Build the physical tensor using independent Clifford sandwiches and
		// scalar outer products. This isolates export from eigensolver defects.
		Mat44 expected = Mat44::sZero();
		for (int axis = 0; axis < 4; ++axis)
		{
			Cl4Ref::MV basis;
			basis.mC[Cl4Ref::cVecMask[axis]] = 1.0;
			Cl4Ref::MV rotated = Cl4Ref::Sandwich(Cl4Ref::FromRotor(principal), basis);
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					expected(i, j) += float(d[axis] * rotated.mC[Cl4Ref::cVecMask[i]] * rotated.mC[Cl4Ref::cVecMask[j]]);
		}
		sCheckTensor(body.GetBodyCreationSettings().mMassPropertiesOverride.mInertia, expected);
		sCheckAngularResponse(body, expected);
	}

	TEST_CASE("BodyCreationSettingsPreserveStaticAndKinematicLatentMass")
	{
		PhysicsTestContext context;
		for (EMotionType type : { EMotionType::Static, EMotionType::Kinematic })
		{
			BodyCreationSettings settings(new SphereShape(0.5f), RVec4(4, 3, 2, 1),
				Rotor::sExp(Bivec(-0.1f, 0.2f, 0.3f, -0.4f, 0.1f, 0.2f)), type, Layers::MOVING);
			settings.mAllowDynamicOrKinematic = true;
			settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
			settings.mMassPropertiesOverride = sParticleProperties();
			Body &original = context.CreateBody(settings, EActivation::DontActivate);
			BodyCreationSettings exported = original.GetBodyCreationSettings();
			CHECK(exported.mMotionType == type);
			CHECK(exported.mAllowDynamicOrKinematic);
			CHECK(exported.mMassPropertiesOverride.mMass == doctest::Approx(12.0f));
			sCheckTensor(exported.mMassPropertiesOverride.mInertia, settings.mMassPropertiesOverride.mInertia);
			Body &recreated = context.CreateBody(exported, EActivation::DontActivate);
			context.GetSystem()->GetBodyInterface().SetMotionType(original.GetID(), EMotionType::Dynamic, EActivation::DontActivate);
			context.GetSystem()->GetBodyInterface().SetMotionType(recreated.GetID(), EMotionType::Dynamic, EActivation::DontActivate);
			sCheckAngularResponse(original, settings.mMassPropertiesOverride.mInertia);
			sCheckAngularResponse(recreated, settings.mMassPropertiesOverride.mInertia);
		}
		BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Static, Layers::NON_MOVING);
		Body &original = context.CreateBody(settings, EActivation::DontActivate);
		BodyCreationSettings exported = original.GetBodyCreationSettings();
		CHECK_FALSE(exported.mAllowDynamicOrKinematic);
		CHECK(exported.mMassPropertiesOverride.mMass == 1.0f);
		sCheckTensor(exported.mMassPropertiesOverride.mInertia, Mat44::sZero());
		Body &recreated = context.CreateBody(exported, EActivation::DontActivate);
		CHECK(recreated.IsStatic());
		CHECK(recreated.GetMotionPropertiesUnchecked() == nullptr);
	}

	TEST_CASE("BodyCreationSettingsPreserveZeroParticlePlanarMoments")
	{
		PhysicsTestContext context;
		BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
		Body &body = context.CreateBody(settings, EActivation::DontActivate);
		body.GetMotionProperties()->SetInverseInertia(Bivec(1.0f / 3, 1.0f / 3, 1.0f / 3, 0, 0, 0), Rotor::sIdentity());
		Mat44 expected = Mat44::sZero();
		expected(0, 0) = 3.0f; // Two equal particles on the X axis.
		BodyCreationSettings exported = body.GetBodyCreationSettings();
		sCheckTensor(exported.mMassPropertiesOverride.mInertia, expected);
		Body &recreated = context.CreateBody(exported, EActivation::DontActivate);
		sCheckTensor(recreated.GetBodyCreationSettings().mMassPropertiesOverride.mInertia, expected);
	}

	TEST_CASE("SphereSecondMomentMatchesRadialParticleQuadrature")
	{
		// Independent volume quadrature: shell measure 2*pi^2*r^3 dr and
		// eight equally weighted +/- coordinate-axis particles on each S^3.
		// This angular rule integrates every quadratic exactly; midpoint radial
		// integration converges to the solid ball, not its surface shell.
		for (float radius : { 0.7f, 2.3f })
			for (float density : { 1.0f, 17.0f })
			{
				SphereShapeSettings settings(radius);
				settings.mDensity = density;
				RefConst<Shape> shape = settings.Create().Get();
				MassProperties properties = shape->GetMassProperties();
				double mass = 0.0, second[4][4] = {}, planar[6] = {};
				constexpr int steps = 512;
				const double pi = std::acos(-1.0), dr = double(radius) / steps;
				for (int shell = 0; shell < steps; ++shell)
				{
					double r = (shell + 0.5) * dr;
					double particle_mass = density * 2.0 * pi * pi * r * r * r * dr / 8.0;
					for (int axis = 0; axis < 4; ++axis)
						for (int sign : { -1, 1 })
						{
							double x[4] = {};
							x[axis] = sign * r;
							mass += particle_mass;
							for (int i = 0; i < 4; ++i)
								for (int j = 0; j < 4; ++j)
									second[i][j] += particle_mass * x[i] * x[j];
							int plane = 0;
							for (int i = 0; i < 4; ++i)
								for (int j = i + 1; j < 4; ++j)
									planar[plane++] += particle_mass * (x[i] * x[i] + x[j] * x[j]);
						}
				}
				CAPTURE(radius); CAPTURE(density);
				CHECK(properties.mMass == doctest::Approx(mass).epsilon(1.0e-5));
				for (float scale : { 1.0f, 3.25f })
				{
					MassProperties scaled = properties;
					scaled.ScaleToMass(properties.mMass * scale);
					int plane = 0;
					for (int i = 0; i < 4; ++i)
					{
						CHECK(scaled.mInertia(i, i) == doctest::Approx(scaled.mMass * radius * radius / 6.0).epsilon(1.0e-6));
						for (int j = 0; j < 4; ++j)
							CHECK(scaled.mInertia(i, j) == doctest::Approx(scale * second[i][j]).epsilon(1.0e-5));
						for (int j = i + 1; j < 4; ++j)
							CHECK(scaled.mInertia(i, i) + scaled.mInertia(j, j) == doctest::Approx(scale * planar[plane++]).epsilon(1.0e-5));
					}
				}
			}
	}
}
