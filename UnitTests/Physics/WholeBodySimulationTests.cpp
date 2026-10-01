// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include "PhysicsTestContext.h"
#include "Layers.h"
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>

// TODO(4D): whole-scene tests. Everything else in this suite checks one stage in isolation - a
// product table, a closest point, a single constraint solve. These drive the real PhysicsSystem
// update loop, which is the only place integration order, islands, sleeping, the contact cache and
// the manifold pipeline are exercised together, and the only place a conservation law can be
// stated at all. They are also where an asymmetry between W and the other three axes shows up:
// almost every formula in the port is written to look dimension-agnostic, so a dropped fourth
// component usually compiles, passes its own unit test, and only appears as a body drifting in a
// direction nothing pushes it.
TEST_SUITE("WholeBodySimulationTests")
{
	using namespace JPH::literals;

	// Forward inertia operator I = Ad(R) D Ad(R^-1), the inverse of
	// MotionProperties::MultiplyWorldSpaceInverseInertiaByVector. Needed to state angular momentum;
	// the body only ever needs I^-1 so it does not store this.
	static Bivec sMultiplyWorldSpaceInertiaByVector(const Body &inBody, BivecArg inV)
	{
		const MotionProperties *mp = inBody.GetMotionProperties();
		Rotor combined = inBody.GetRotation() * mp->GetInertiaRotation();
		Bivec inv_d = mp->GetInverseInertiaDiagonal();
		Bivec local = combined.InverseRotateBivec(inV);
		Bivec scaled(local[0] / inv_d[0], local[1] / inv_d[1], local[2] / inv_d[2],
					 local[3] / inv_d[3], local[4] / inv_d[4], local[5] / inv_d[5]);
		return combined.RotateBivec(scaled);
	}

	// Angular momentum about the origin: L = I * Omega + m * (v ^ x). The wedge operand order is the
	// one fixed by DynamicsConventionsTest - see the torque convention note in Bivec.h.
	static Bivec sAngularMomentum(const Body &inBody)
	{
		const MotionProperties *mp = inBody.GetMotionProperties();
		Bivec spin = sMultiplyWorldSpaceInertiaByVector(inBody, inBody.GetAngularVelocity());
		Bivec orbit = Bivec::sWedge(inBody.GetLinearVelocity(), Vec4(inBody.GetCenterOfMassPosition())) * (1.0f / mp->GetInverseMass());
		return spin + orbit;
	}

	static float sKineticEnergy(const Body &inBody)
	{
		const MotionProperties *mp = inBody.GetMotionProperties();
		Bivec w = inBody.GetAngularVelocity();
		return 0.5f * (inBody.GetLinearVelocity().LengthSq() / mp->GetInverseMass() + w.Dot(sMultiplyWorldSpaceInertiaByVector(inBody, w)));
	}

	// A body in free fall must follow the integrator's own prediction exactly, in every one of the
	// four directions. Running the same scene once per axis is what catches a W that is handled
	// differently from X / Y / Z - the failure mode this whole port is most exposed to.
	TEST_CASE("TestFreeFallAlongEachAxis")
	{
		for (int axis = 0; axis < 4; ++axis)
		{
			PhysicsTestContext c;

			Vec4 gravity = Vec4::sZero();
			gravity[axis] = -9.81f;
			c.GetSystem()->SetGravity(gravity);

			// Gravity really is 4D now; a Vec3 would have mirrored Z into W
			CHECK(c.GetSystem()->GetGravity() == gravity);

			const RVec4 cInitialPos(1.0_r, 2.0_r, 3.0_r, 4.0_r);
			Body &body = c.CreateBox(cInitialPos, Rotor::sIdentity(), EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4::sReplicate(0.5f));

			const float cTime = 1.0f;
			c.Simulate(cTime);

			RVec4 expected = c.PredictPosition(cInitialPos, Vec4::sZero(), gravity, cTime);
			CHECK_APPROX_EQUAL(body.GetPosition(), expected, 1.0e-4f);

			// It must not pick up a rotation out of nowhere either
			CHECK(body.GetAngularVelocity().Length() < 1.0e-6f);
		}
	}

	// Free rotation. A 4D angular velocity need not lie in a single plane, so this also exercises
	// the general rotor exponential rather than the simple-rotation special case.
	TEST_CASE("TestFreeRotationMatchesIntegrator")
	{
		PhysicsTestContext c;
		c.ZeroGravity();

		// Deliberately not a single rotation plane, and deliberately not inertially symmetric, so
		// the gyroscopic term is active
		const Bivec cAngularVelocity(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.7f);

		Body &body = c.CreateBox(RVec4::sZero(), Rotor::sIdentity(), EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4(0.5f, 0.5f, 0.5f, 0.5f));
		body.SetAngularVelocity(cAngularVelocity);

		const float cTime = 0.5f;
		c.Simulate(cTime);

		// An inertially symmetric body (a tesseract with equal extents is) has no gyroscopic torque,
		// so the angular velocity is unchanged and the orientation is the free integration of it
		CHECK_APPROX_EQUAL(body.GetAngularVelocity(), cAngularVelocity, 1.0e-4f);
		Rotor expected = c.PredictOrientation(Rotor::sIdentity(), cAngularVelocity, Bivec::sZero(), cTime);
		CHECK_APPROX_EQUAL(body.GetRotation(), expected, 1.0e-4f);
	}

	// A tesseract dropped flat lands on a cubic cell. This is the case the 4D contact manifold
	// exists for: in 3D the contact patch is a polygon, here it is a 3D polytope, and a single
	// contact point would let the body tip over.
	TEST_CASE("TestTesseractRestsOnFloor")
	{
		PhysicsTestContext c;
		c.CreateFloor();

		const float cHalfExtent = 0.5f;
		const RVec4 cInitialPos(0.0_r, 2.0_r, 0.0_r, 0.0_r);
		Body &body = c.CreateBox(cInitialPos, Rotor::sIdentity(), EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4::sReplicate(cHalfExtent));

		c.Simulate(2.0f);

		// Resting on the floor, within the solver's allowed penetration
		RVec4 pos = body.GetPosition();
		CHECK(abs(float(pos.GetY()) - cHalfExtent) < 1.1f * c.GetSystem()->GetPhysicsSettings().mPenetrationSlop);

		// Nothing acts along X, Z or W. Any motion there is the fourth dimension leaking.
		CHECK(abs(float(pos.GetX())) < 1.0e-3f);
		CHECK(abs(float(pos.GetZ())) < 1.0e-3f);
		CHECK(abs(float(pos.GetW())) < 1.0e-3f);

		// A flat landing on a flat floor must not induce a spin, in any of the six planes
		CHECK(body.GetAngularVelocity().Length() < 1.0e-3f);

		// And it must settle rather than jitter
		CHECK(!body.IsActive());
	}

	// The manifold must actually span the contact cell, not collapse to a point. A single contact
	// point under a body still holds it up, so the resting-height check above passes either way;
	// what it cannot do is resist a torque. Press down on one corner of the top cell and check the
	// tesseract stays flat. (Verified to discriminate: forcing the single-point fallback in
	// ManifoldBetweenTwoFaces tips it through 0.40 rad.)
	TEST_CASE("TestRestingTesseractResistsTipping")
	{
		PhysicsTestContext c;
		c.CreateFloor();

		const float cHalfExtent = 0.5f;
		Body &body = c.CreateBox(RVec4(0.0_r, cHalfExtent, 0.0_r, 0.0_r), Rotor::sIdentity(), EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4::sReplicate(cHalfExtent));
		c.Simulate(0.5f);

		// Go through the body interface so the force also wakes the body up - it has settled and
		// fallen asleep by now, and Body::AddForce on a sleeping body is silently ignored
		BodyInterface &bi = c.GetBodyInterface();
		const BodyID id = body.GetID();
		const float cWeight = 9.81f / body.GetMotionProperties()->GetInverseMass();
		c.Simulate(0.5f, [&bi, id, cHalfExtent, cWeight]() {
			bi.AddForce(id, Vec4(0.0f, -20.0f * cWeight, 0.0f, 0.0f), RVec4(Real(cHalfExtent), Real(cHalfExtent), Real(cHalfExtent), Real(cHalfExtent)));
		});

		// Still flat on the floor. The rotor's scalar part is cos of half the rotation angle (for
		// each of the two SU(2) factors), so this bounds the tilt to about 0.09 rad.
		CHECK(abs(abs(body.GetRotation().GetScalar()) - 1.0f) < 1.0e-3f);
		CHECK(abs(float(body.GetPosition().GetY()) - cHalfExtent) < 0.05f);
	}

	// Stacking is the standard stability benchmark, and the one the manifold reduction is for: the
	// patch between two tesseracts is 3D, so the points kept must span it in all three directions
	// or the stack slowly shears apart.
	TEST_CASE("TestTesseractStackIsStable")
	{
		PhysicsTestContext c;
		c.CreateFloor();

		const float cHalfExtent = 0.5f;
		const int cStackHeight = 5;

		Body *bodies[cStackHeight];
		for (int i = 0; i < cStackHeight; ++i)
		{
			// A small gap so they settle onto each other rather than starting interpenetrating
			Real y = Real(cHalfExtent + i * (2.0f * cHalfExtent + 0.01f));
			bodies[i] = &c.CreateBox(RVec4(0.0_r, y, 0.0_r, 0.0_r), Rotor::sIdentity(), EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4::sReplicate(cHalfExtent));
		}

		c.Simulate(3.0f);

		for (int i = 0; i < cStackHeight; ++i)
		{
			RVec4 pos = bodies[i]->GetPosition();

			// Each body is still one box above the one below it (allowing the solver's penetration
			// at every interface below it)
			float expected_y = cHalfExtent + i * 2.0f * cHalfExtent;
			CHECK(abs(float(pos.GetY()) - expected_y) < (i + 1) * 1.1f * c.GetSystem()->GetPhysicsSettings().mPenetrationSlop);

			// And has not slid out of the stack in any of the three horizontal directions
			CHECK(abs(float(pos.GetX())) < 0.05f);
			CHECK(abs(float(pos.GetZ())) < 0.05f);
			CHECK(abs(float(pos.GetW())) < 0.05f);
		}
	}

	// A 4D pyramid. The 5-high column above is one contact per interface; a pyramid is a single
	// large island where every body is held by several neighbours at once, which is what shakes out
	// a manifold that spans the patch unevenly, and it is the only test here that drives EPA hard
	// enough to reach its fallback path (see sCollideConvexVsConvex).
	TEST_CASE("TestTesseractPyramidIsStable")
	{
		PhysicsTestContext c;
		c.CreateFloor();

		// A 4D pyramid has three horizontal axes, so its body count is a sum of cubes: height 4
		// gives (4*5/2)^2 = 100 bodies, which is enough to be a real island and still quick
		const float cHalfExtent = 0.5f;
		const float cSize = 2.0f * cHalfExtent;
		const int cHeight = 4;

		Array<Body *> bodies;
		Array<int> level;
		for (int i = 0; i < cHeight; ++i)
			for (int j = i / 2; j < cHeight - (i + 1) / 2; ++j)
				for (int k = i / 2; k < cHeight - (i + 1) / 2; ++k)
					for (int l = i / 2; l < cHeight - (i + 1) / 2; ++l)
					{
						float offset = (i & 1)? cHalfExtent : 0.0f;
						// A small vertical gap so they settle onto each other rather than starting interpenetrating
						RVec4 pos(Real(-cHeight + cSize * j + offset), Real(cHalfExtent + (cSize + 0.01f) * i), Real(-cHeight + cSize * k + offset), Real(-cHeight + cSize * l + offset));
						bodies.push_back(&c.CreateBox(pos, Rotor::sIdentity(), EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4::sReplicate(cHalfExtent)));
						level.push_back(i);
					}
		CHECK(bodies.size() == 100);

		c.Simulate(3.0f);

		const float cSlop = c.GetSystem()->GetPhysicsSettings().mPenetrationSlop;
		for (size_t n = 0; n < bodies.size(); ++n)
		{
			RVec4 pos = bodies[n]->GetPosition();

			// Nothing sank through the floor or through the body below it. Allow the solver's
			// penetration at each of the interfaces underneath, plus a little settling.
			float expected_y = cHalfExtent + level[n] * cSize;
			CHECK(float(pos.GetY()) > expected_y - (level[n] + 1) * 2.0f * cSlop);
			CHECK(float(pos.GetY()) < expected_y + cSlop);
		}
	}

	// Conservation laws across a real collision. These cannot be checked at the level of a single
	// constraint part: they are statements about the whole solve.
	TEST_CASE("TestCollisionConservesMomentum")
	{
		PhysicsTestContext c;
		c.ZeroGravity();

		// Offset the impact in W as well as Z so the collision is genuinely four dimensional and
		// imparts angular velocity in planes containing W
		Body &body1 = c.CreateBox(RVec4(-2.0_r, 0.0_r, 0.1_r, 0.15_r), Rotor::sIdentity(), EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4::sReplicate(0.5f));
		Body &body2 = c.CreateBox(RVec4(2.0_r, 0.0_r, -0.1_r, -0.15_r), Rotor::sIdentity(), EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4::sReplicate(0.5f));
		body1.SetLinearVelocity(Vec4(2.0f, 0.0f, 0.0f, 0.0f));
		body2.SetLinearVelocity(Vec4(-2.0f, 0.0f, 0.0f, 0.0f));

		float mass = 1.0f / body1.GetMotionProperties()->GetInverseMass();
		Vec4 p_before = mass * (body1.GetLinearVelocity() + body2.GetLinearVelocity());
		Bivec l_before = sAngularMomentum(body1) + sAngularMomentum(body2);
		float e_before = sKineticEnergy(body1) + sKineticEnergy(body2);

		// Long enough for them to meet and separate again
		c.Simulate(2.0f);

		Vec4 p_after = mass * (body1.GetLinearVelocity() + body2.GetLinearVelocity());
		Bivec l_after = sAngularMomentum(body1) + sAngularMomentum(body2);
		float e_after = sKineticEnergy(body1) + sKineticEnergy(body2);

		// They did collide, otherwise this proves nothing
		CHECK(e_after < 0.99f * e_before);

		// Contacts are internal forces, so they can only move momentum between the two bodies.
		// Angular momentum is checked relative to its own magnitude: |L| is 721 here (the bodies are
		// 1000 kg tesseracts two metres off the axis), and the measured drift is 1.2e-3, i.e. 1.7e-6
		// relative - float round-off over 120 steps, not a leak. Linear momentum comes out exact.
		CHECK_APPROX_EQUAL(p_after, p_before, 1.0e-3f);
		CHECK(Bivec(l_after - l_before).Length() < 1.0e-5f * l_before.Length());

		// And with the default restitution of 0 the collision can only take energy out
		CHECK(e_after <= e_before + 1.0e-3f);
	}

	// Determinism. Jolt's central guarantee, and in 4D it is also a check that nothing reads
	// uninitialized memory in the new code paths - the EPA cells and the clipped manifold are both
	// built in scratch buffers that are not cleared between steps.
	TEST_CASE("TestSceneIsDeterministic")
	{
		auto run = [](Array<RVec4> &outPositions, Array<Rotor> &outRotations) {
			PhysicsTestContext c(1.0f / 60.0f, 1, 4);
			c.CreateFloor();

			Array<Body *> bodies;
			for (int i = 0; i < 6; ++i)
			{
				// Spread out in W as well, and start them rotated so the solve is not trivially symmetric
				RVec4 pos(Real(0.3f * i), Real(1.0f + 1.5f * i), Real(-0.2f * i), Real(0.25f * i));
				Rotor rot = Rotor::sExp(Bivec(0.1f * i, 0.0f, 0.05f * i, 0.0f, 0.0f, 0.0f));
				bodies.push_back(&c.CreateBox(pos, rot, EMotionType::Dynamic, EMotionQuality::Discrete, Layers::MOVING, Vec4::sReplicate(0.5f)));
			}

			c.Simulate(2.0f);

			for (Body *b : bodies)
			{
				outPositions.push_back(b->GetPosition());
				outRotations.push_back(b->GetRotation());
			}
		};

		Array<RVec4> pos1, pos2;
		Array<Rotor> rot1, rot2;
		run(pos1, rot1);
		run(pos2, rot2);

		CHECK(pos1.size() == pos2.size());
		for (size_t i = 0; i < pos1.size(); ++i)
		{
			// Bit exact, not approximately equal - that is the whole point
			CHECK(pos1[i] == pos2[i]);
			CHECK(rot1[i] == rot2[i]);
		}
	}
}
