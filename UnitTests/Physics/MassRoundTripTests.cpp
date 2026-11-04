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
#include <cstdlib>
#include <vector>
#if defined(JPH_PLATFORM_LINUX)
#include <spawn.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <cerrno>
#include <csignal>
extern char **environ;
#endif

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

	TEST_CASE("RotatedTwoParticleMassDoesNotInventInverseInertia")
	{
		PhysicsTestContext context;
		BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(),
			Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
		settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
		settings.mMassPropertiesOverride.mMass = 2.0f;
		settings.mMassPropertiesOverride.mInertia = Mat44::sZero();
		// Two unit particles at +/-(0,0,1,1), not inverse-inertia setter overrides.
		settings.mMassPropertiesOverride.mInertia(2, 2) = settings.mMassPropertiesOverride.mInertia(2, 3)
			= settings.mMassPropertiesOverride.mInertia(3, 2) = settings.mMassPropertiesOverride.mInertia(3, 3) = 2.0f;
		Body &original = context.CreateBody(settings, EActivation::DontActivate);
		// A float matrix cannot generally encode an oblique null space.
		// The checked API must reject this state atomically, not invent a
		// huge inverse or freeze the body by silently changing its DOFs.
		BodyCreationSettings exported;
		exported.mUserData = 0x12345678;
		CHECK_FALSE(original.TryGetBodyCreationSettings(exported));
		CHECK(exported.mUserData == 0x12345678);
		CHECK(exported.GetShape() == nullptr);
		CHECK(original.GetMotionProperties()->GetAllowedDOFs() == EAllowedDOFs::All);
		for (const Body *body : { &original })
			for (int input = 0; input < 6; ++input)
			{
				float v[6] = {};
				v[input] = 1.0f;
				Bivec actual = body->GetMotionProperties()->MultiplyWorldSpaceInverseInertiaByVector(
					body->GetRotation(), Bivec(v[0], v[1], v[2], v[3], v[4], v[5]));
				// For n=(0,0,1,1)/sqrt(2), I=4 P on the image of
				// omega -> omega.Contract(n). Thus I^+=P/4: XZ/XW
				// and YZ/YW each have a 2x2 block of 1/8, ZW is 1/4.
				for (int output = 0; output < 6; ++output)
				{
					double expected = ((input == 1 || input == 2) && (output == 1 || output == 2))
						|| ((input == 3 || input == 4) && (output == 3 || output == 4))? 0.125
						: input == 5 && output == 5? 0.25 : 0.0;
					CAPTURE(input); CAPTURE(output);
					CHECK(std::abs(actual[output] - expected) < 1.0e-7);
				}
			}
	}

	TEST_CASE("ExactRankOneParticlesPreserveNullityUnderOrientationAndScale")
	{
		PhysicsTestContext context;
		const Vec4 directions[] = { Vec4(1, 0, 0, 0), Vec4(0, 1, 0, 0), Vec4(0, 0, 1, 1),
			Vec4(1, -1, 1, -1), Vec4(1, 2, -1, 0), Vec4(0, 0, 0, 2) };
		const Rotor orientations[] = { Rotor::sIdentity(), Rotor::sExp(Bivec(0.23f, -0.31f, 0.17f, 0.29f, -0.11f, 0.37f)) };
		for (Vec4 direction : directions)
			for (float scale : { 1.0f, 4.0f })
				for (Rotor orientation : orientations)
				{
					BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), orientation, EMotionType::Dynamic, Layers::MOVING);
					settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
					settings.mMassPropertiesOverride.mMass = 2.0f;
					settings.mMassPropertiesOverride.mInertia = Mat44::sZero();
					for (int i = 0; i < 4; ++i)
						for (int j = 0; j < 4; ++j)
							settings.mMassPropertiesOverride.mInertia(i, j) = 2.0f * Square(scale) * direction[i] * direction[j];
					Body &body = context.CreateBody(settings, EActivation::DontActivate);
					int zeros = 0;
					for (int plane = 0; plane < 6; ++plane)
						zeros += body.GetMotionProperties()->GetInverseInertiaDiagonal()[plane] == 0.0f;
					CHECK(zeros == 3);
					// Independent world particle position and velocity Gram matrix:
					// I = 2 C(x)^T C(x), I^+ = C(x)^T C(x) / (2 |x|^4).
					Cl4Ref::MV point;
					for (int i = 0; i < 4; ++i)
						point.mC[Cl4Ref::cVecMask[i]] = scale * direction[i];
					point = Cl4Ref::Sandwich(Cl4Ref::FromRotor(orientation), point);
					double x[4], r_sq = 0.0, field[6][4] = {};
					for (int i = 0; i < 4; ++i)
					{
						x[i] = point.mC[Cl4Ref::cVecMask[i]];
						r_sq += x[i] * x[i];
					}
					int plane = 0;
					for (int i = 0; i < 4; ++i)
						for (int j = i + 1; j < 4; ++j, ++plane)
						{
							field[plane][i] = x[j];
							field[plane][j] = -x[i];
						}
					BodyCreationSettings exported;
					bool accepted = body.TryGetBodyCreationSettings(exported);
					int occupied_axes = 0;
					for (int i = 0; i < 4; ++i)
						occupied_axes += direction[i] != 0.0f;
					if (occupied_axes == 1)
						CHECK(accepted);
					const Body *recreated = accepted? &context.CreateBody(exported, EActivation::DontActivate) : nullptr;
					for (const Body *candidate : { static_cast<const Body *>(&body), recreated })
						if (candidate != nullptr)
							for (int input = 0; input < 6; ++input)
							{
								float v[6] = {};
								v[input] = 1.0f;
								Bivec actual = candidate->GetMotionProperties()->MultiplyWorldSpaceInverseInertiaByVector(
									candidate->GetRotation(), Bivec(v[0], v[1], v[2], v[3], v[4], v[5]));
								for (int output = 0; output < 6; ++output)
								{
									double expected = 0.0;
									for (int i = 0; i < 4; ++i)
										expected += field[input][i] * field[output][i] / (2.0 * r_sq * r_sq);
									CAPTURE(direction); CAPTURE(scale); CAPTURE(input); CAPTURE(output); CAPTURE(accepted);
									CHECK(std::abs(actual[output] - expected) < 1.0e-7);
								}
							}
				}
	}

#if defined(JPH_PLATFORM_LINUX)
	// Only the parent test enables this case in a fresh executable. A mere
	// nonzero child exit is NOT success: only the legacy API's SIGABRT is.
	TEST_CASE("LegacyUnrepresentableMassExportDeathChild" * doctest::skip())
	{
		// Even --no-skip on the normal suite must not run an unrequested death.
		const char *marker = std::getenv("JPH_MASS_EXPORT_DEATH_CHILD");
		if (marker == nullptr || marker[0] != '1' || marker[1] != '\0')
			return;
		const rlimit no_core { 0, 0 };
		if (setrlimit(RLIMIT_CORE, &no_core) != 0)
			std::_Exit(92);
		JPH_IF_ENABLE_ASSERTS(AssertFailed = [](const char *, const char *, const char *, uint) { std::_Exit(91); return false; };)
		PhysicsTestContext context;
		BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
		settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
		settings.mMassPropertiesOverride.mMass = 2.0f;
		Mat44 &inertia = settings.mMassPropertiesOverride.mInertia;
		inertia = Mat44::sZero();
		inertia(2, 2) = inertia(2, 3) = inertia(3, 2) = inertia(3, 3) = 2.0f;
		Body &body = context.CreateBody(settings, EActivation::DontActivate);
		(void)body.GetBodyCreationSettings();
		CHECK(false); // Returning invented settings must not pass the death test.
	}

	TEST_CASE("LegacyUnrepresentableMassExportFailsClosedInSubprocess")
	{
		char executable[] = "/proc/self/exe";
		char filter[] = "--test-case=LegacyUnrepresentableMassExportDeathChild";
		char unskip[] = "--no-skip=true";
		char no_colors[] = "--no-colors=true";
		char quiet[] = "--quiet=true"; // Expected child SIGABRT is checked by waitpid, not its reporter.
		char *args[] = { executable, filter, unskip, no_colors, quiet, nullptr };
		char marker[] = "JPH_MASS_EXPORT_DEATH_CHILD=1";
		std::vector<char *> child_env;
		for (char **env = environ; *env != nullptr; ++env)
			child_env.push_back(*env);
		child_env.push_back(marker);
		child_env.push_back(nullptr);
		pid_t child;
		int spawned = posix_spawn(&child, executable, nullptr, nullptr, args, child_env.data());
		CHECK(spawned == 0);
		if (spawned != 0)
			return;
		int status = 0;
		pid_t waited;
		do { waited = waitpid(child, &status, 0); } while (waited == -1 && errno == EINTR);
		CHECK(waited == child);
		CHECK(WIFSIGNALED(status));
		if (WIFSIGNALED(status))
			CHECK(WTERMSIG(status) == SIGABRT);
	}
#endif

	TEST_CASE("ExactRankOneStructuralZerosSurviveLargeScalesAndEveryAxis")
	{
		PhysicsTestContext context;
		JPH_IF_ENABLE_ASSERTS(ExpectAssert no_asserts(0);)
		for (int zero_axis = 0; zero_axis < 4; ++zero_axis)
			for (int exponent : { -40, 0, 40 })
			{
				CAPTURE(zero_axis); CAPTURE(exponent);
				// This permutation family includes the reported x=(0,-2,10,-10),
				// scale=2^40 fixture. Every entry is exact binary32.
				float x[] = { 0, -2, 10, -10 };
				std::swap(x[0], x[zero_axis]);
				BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
				settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
				MassProperties &mass = settings.mMassPropertiesOverride;
				mass.mMass = 2.0f;
				for (int i = 0; i < 4; ++i)
					for (int j = 0; j < 4; ++j)
						mass.mInertia(i, j) = std::ldexp(1.0f, exponent) * x[i] * x[j];
				Mat44 rotation;
				Vec4 diagonal;
				bool decomposed = mass.DecomposePrincipalMomentsOfInertia(rotation, diagonal);
				CHECK(decomposed);
				if (!decomposed)
					return;
				CHECK(rotation(zero_axis, 0) == 0.0f);
				CHECK(diagonal[0] == std::ldexp(204.0f, exponent));
				CHECK(diagonal[1] == 0.0f);
				CHECK(diagonal[2] == 0.0f);
				CHECK(diagonal[3] == 0.0f);
				CHECK(rotation.GetDeterminant() == doctest::Approx(1.0f).epsilon(1.0e-6));
				Mat44 diagonal_matrix = Mat44::sZero();
				diagonal_matrix(0, 0) = diagonal[0];
				Mat44 reconstructed = rotation * diagonal_matrix * rotation.Transposed();
				for (int i = 0; i < 4; ++i)
					for (int j = 0; j < 4; ++j)
					{
						double dot = 0.0;
						for (int k = 0; k < 4; ++k)
							dot += double(rotation(k, i)) * rotation(k, j);
						CHECK(std::abs(dot - (i == j? 1.0 : 0.0)) < 2.0e-7);
						if (i == zero_axis || j == zero_axis)
							CHECK(reconstructed(i, j) == 0.0f);
					}
				Body &body = context.CreateBody(settings, EActivation::DontActivate);
				int zeros = 0;
				for (int plane = 0; plane < 6; ++plane)
					zeros += body.GetMotionProperties()->GetInverseInertiaDiagonal()[plane] == 0.0f;
				CHECK(zeros == 3);
			}
	}

	TEST_CASE("CheckedExportRejectsUnsupportedEigenpairAccuracyWithoutAssertions")
	{
		PhysicsTestContext context;
		BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
		settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
		settings.mMassPropertiesOverride.mMass = 2.0f;
		// Exactly represented SPD: x x^T + 2^-8 Identity. Initial creation
		// succeeds, but decoding its rounded export fails the eigenpair check.
		const float x[] = { -4, -3, 3, -3 };
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				settings.mMassPropertiesOverride.mInertia(i, j) = x[i] * x[j] + (i == j? std::ldexp(1.0f, -8) : 0.0f);
		JPH_IF_ENABLE_ASSERTS(ExpectAssert no_asserts(0);)
		Body &body = context.CreateBody(settings, EActivation::DontActivate);
		BodyCreationSettings exported(new SphereShape(0.75f), RVec4(9, 8, 7, 6), Rotor::sIdentity(), EMotionType::Kinematic, Layers::MOVING);
		exported.mUserData = 0xfeedbeef;
		exported.mMassPropertiesOverride.mMass = 53.0f;
		const Shape *shape = exported.GetShape();
		const Bivec inverse = body.GetMotionProperties()->GetInverseInertiaDiagonal();
		const Rotor principal = body.GetMotionProperties()->GetInertiaRotation();
		CHECK_FALSE(body.TryGetBodyCreationSettings(exported));
		CHECK(exported.GetShape() == shape);
		CHECK(exported.mPosition == RVec4(9, 8, 7, 6));
		CHECK(exported.mMotionType == EMotionType::Kinematic);
		CHECK(exported.mUserData == 0xfeedbeef);
		CHECK(exported.mMassPropertiesOverride.mMass == 53.0f);
		CHECK(body.GetMotionProperties()->GetInverseInertiaDiagonal() == inverse);
		CHECK(body.GetMotionProperties()->GetInertiaRotation() == principal);
	}

	TEST_CASE("DiagonalRankTwoMassPreservesEveryChoiceOfNullPlane")
	{
		PhysicsTestContext context;
		for (int a = 0; a < 4; ++a)
			for (int b = a + 1; b < 4; ++b)
			{
				BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
				settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
				settings.mMassPropertiesOverride.mMass = 2.0f;
				settings.mMassPropertiesOverride.mInertia = Mat44::sZero();
				settings.mMassPropertiesOverride.mInertia(a, a) = 2.0f;
				settings.mMassPropertiesOverride.mInertia(b, b) = 7.0f;
				Body &original = context.CreateBody(settings, EActivation::DontActivate);
				BodyCreationSettings exported;
				bool accepted = original.TryGetBodyCreationSettings(exported);
				CHECK(accepted);
				if (!accepted)
					continue;
				Body &recreated = context.CreateBody(exported, EActivation::DontActivate);
				int plane = 0;
				for (int i = 0; i < 4; ++i)
					for (int j = i + 1; j < 4; ++j, ++plane)
					{
						float moment = settings.mMassPropertiesOverride.mInertia(i, i) + settings.mMassPropertiesOverride.mInertia(j, j);
						float expected = moment > 0.0f? 1.0f / moment : 0.0f;
						CHECK(recreated.GetMotionProperties()->GetInverseInertiaDiagonal()[plane] == doctest::Approx(expected).epsilon(1.0e-6));
						if (moment == 0.0f)
							CHECK(recreated.GetMotionProperties()->GetInverseInertiaDiagonal()[plane] == 0.0f);
					}
			}
	}

	TEST_CASE("SkinnyPositiveSecondMomentsAreNeverClassifiedAsNull")
	{
		PhysicsTestContext context;
		for (int exponent : { 4, 12, 20, 40, 80 })
		{
			float small = std::ldexp(1.0f, -exponent);
			BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
			settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
			settings.mMassPropertiesOverride.mMass = 2.0f;
			settings.mMassPropertiesOverride.mInertia = Mat44::sZero();
			settings.mMassPropertiesOverride.mInertia(0, 0) = 4.0f;
			for (int i = 1; i < 4; ++i)
				settings.mMassPropertiesOverride.mInertia(i, i) = small;
			Body &original = context.CreateBody(settings, EActivation::DontActivate);
			Bivec before = original.GetMotionProperties()->GetInverseInertiaDiagonal();
			for (int plane = 0; plane < 6; ++plane)
			{
				CAPTURE(exponent); CAPTURE(plane);
				CHECK(before[plane] > 0.0f);
				CHECK(before[plane] == doctest::Approx(1.0f / (plane < 3? 4.0f + small : 2.0f * small)));
			}
			BodyCreationSettings exported;
			bool accepted = original.TryGetBodyCreationSettings(exported);
			CAPTURE(exponent); CAPTURE(accepted);
			if (exponent <= 20)
				CHECK(accepted);
			// Unrecoverable reciprocal cancellation may reject an extreme
			// aspect ratio; it must NEVER silently turn a positive mode off.
			if (accepted)
			{
				Body &recreated = context.CreateBody(exported, EActivation::DontActivate);
				Bivec after = recreated.GetMotionProperties()->GetInverseInertiaDiagonal();
				for (int plane = 0; plane < 6; ++plane)
				{
					CHECK(after[plane] > 0.0f);
					CHECK(after[plane] == doctest::Approx(before[plane]).epsilon(1.0e-6));
				}
			}
		}
	}

	TEST_CASE("ZeroAndSmallPhysicalMomentsDoNotAcquireSphereFallbackInertia")
	{
		PhysicsTestContext context;
		for (float scale : { 0.0f, std::ldexp(1.0f, -40), std::ldexp(1.0f, -60) })
		{
			BodyCreationSettings settings(new SphereShape(0.5f), RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
			settings.mOverrideMassProperties = EOverrideMassProperties::MassAndInertiaProvided;
			settings.mMassPropertiesOverride.mMass = 2.0f;
			settings.mMassPropertiesOverride.mInertia = Mat44::sZero();
			for (int i = 0; i < 4; ++i)
				settings.mMassPropertiesOverride.mInertia(i, i) = (i + 1) * scale;
			Body &original = context.CreateBody(settings, EActivation::DontActivate);
			BodyCreationSettings exported;
			bool accepted = original.TryGetBodyCreationSettings(exported);
			CHECK(accepted);
			if (!accepted)
				continue;
			Body &recreated = context.CreateBody(exported, EActivation::DontActivate);
			for (const Body *body : { &original, &recreated })
			{
				int plane = 0;
				for (int i = 0; i < 4; ++i)
					for (int j = i + 1; j < 4; ++j, ++plane)
					{
						float expected = scale == 0.0f? 0.0f : 1.0f / ((i + j + 2) * scale);
						CAPTURE(scale); CAPTURE(plane);
						CHECK(body->GetMotionProperties()->GetInverseInertiaDiagonal()[plane] == doctest::Approx(expected).epsilon(1.0e-6));
					}
			}
		}
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
