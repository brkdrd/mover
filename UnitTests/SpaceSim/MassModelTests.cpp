// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#include <Jolt/SpaceSim/MassModel.h>
#include "Math/Cl4Reference.h"
#include <limits>

#define MASS_REQUIRE(expression) do { bool success = (expression); CHECK(success); if (!success) return; } while (false)

using namespace JPH::SpaceSim;

namespace
{
ComponentMass Component(uint64 id, float mass, MassVector center, MassVector diagonal)
{
    ComponentMass c;
    c.mID = id;
    c.mLocalCOM = center;
    c.mProperties.mMass = mass;
    for (int i = 0; i < 4; ++i)
        c.mProperties.mInertia(i, i) = float(diagonal[size_t(i)]);
    return c;
}

MassVector ReferenceRotate(const Cl4Ref::MV &rotation, const MassVector &v)
{
    Cl4Ref::MV point;
    for (size_t i = 0; i < 4; ++i) point.mC[Cl4Ref::cVecMask[i]] = v[i];
    auto result = Cl4Ref::Sandwich(rotation, point);
    MassVector out;
    for (size_t i = 0; i < 4; ++i) out[i] = result.mC[Cl4Ref::cVecMask[i]];
    return out;
}

MassVector ReferencePointVelocity(const MassBody &body, const MassVector &point, const MassVector &com)
{
    // Deliberately use a full antisymmetric 4x4 matrix, independent of production Contract.
    double w[4][4] {};
    size_t k = 0;
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = i + 1; j < 4; ++j)
        {
            w[i][j] = body.mAngularVelocity[k]; w[j][i] = -w[i][j]; ++k;
        }
    MassVector v = body.mLinearVelocity;
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 4; ++j) v[i] += w[i][j] * (point[j] - com[j]);
    return v;
}

void CheckGlobalConservation(const std::vector<MassBody> &before, const MassTransaction &tx, double tolerance = 1e-7)
{
    MassMomentum a, b;
    auto collect = [](const std::vector<MassBody> &bodies, MassMomentum &sum)
    {
        for (const auto &body : bodies)
        {
            BodyMassData data;
            CHECK(MassModel::Derive(body, data));
            sum.mMass += data.mMomentum.mMass;
            sum.mKineticEnergy += data.mMomentum.mKineticEnergy;
            for (size_t i = 0; i < 4; ++i) sum.mLinear[i] += data.mMomentum.mLinear[i];
            for (size_t i = 0; i < 6; ++i) sum.mAngular[i] += data.mMomentum.mAngular[i];
        }
    };
    collect(before, a); collect(tx.mBodies, b);
    for (const auto &gas : tx.mExhaust)
    {
        b.mMass += gas.mMass;
        for (size_t i = 0; i < 4; ++i)
        {
            b.mLinear[i] += gas.mMass * gas.mWorldVelocity[i];
            b.mKineticEnergy += 0.5 * gas.mMass * gas.mWorldVelocity[i] * gas.mWorldVelocity[i];
        }
        size_t k = 0;
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = i + 1; j < 4; ++j)
                b.mAngular[k++] += gas.mMass * (gas.mWorldVelocity[i] * gas.mWorldPosition[j] - gas.mWorldVelocity[j] * gas.mWorldPosition[i]);
    }
    CHECK_APPROX_EQUAL(a.mMass, b.mMass, tolerance);
    for (size_t i = 0; i < 4; ++i) CHECK_APPROX_EQUAL(a.mLinear[i], b.mLinear[i], tolerance);
    for (size_t i = 0; i < 6; ++i) CHECK_APPROX_EQUAL(a.mAngular[i], b.mAngular[i], tolerance);
    CHECK_APPROX_EQUAL(b.mKineticEnergy - a.mKineticEnergy, tx.mLedger.mAddedEnergy - tx.mLedger.mDissipatedEnergy, tolerance);
}
}

TEST_SUITE("SpaceSimMass")
{
    TEST_CASE("AggregateRejectsInvalidStateAtomically")
    {
        auto good = Component(1, 2, { 0, 0, 0, 0 }, { 1, 2, 3, 4 });
        AggregateMass out;
        out.mMass = 123;
        auto rejected = [&](const std::vector<ComponentMass> &input)
        {
            CHECK_FALSE(MassModel::Aggregate(input, out));
            CHECK(out.mMass == 123);
        };
        rejected({});
        rejected({ good, good });
        auto bad = good; bad.mID = 0; rejected({ bad });
        bad = good; bad.mProperties.mMass = -1; rejected({ bad });
        bad = good; bad.mProperties.mMass = std::numeric_limits<float>::infinity(); rejected({ bad });
        bad = good; bad.mLocalCOM[3] = std::numeric_limits<double>::quiet_NaN(); rejected({ bad });
        bad = good; bad.mLocalCOM[0] = 1e100; rejected({ bad });
        bad = good; bad.mLocalRotation = Rotor(0.70710677f, 0, 0, 0, 0, 0, 0, 0.70710677f); rejected({ bad });
        bad = good; bad.mProperties.mInertia(0, 1) = 2; rejected({ bad });
        bad = good; bad.mProperties.mInertia(2, 2) = -0.1f; rejected({ bad });
        bad = good; bad.mProperties.mInertia(0, 0) = std::numeric_limits<float>::quiet_NaN(); rejected({ bad });
        bad = good; bad.mProperties.mInertia = Mat44::sZero(); rejected({ bad });
        bad = good; bad.mProperties.mInertia(0, 0) = 1e-20f; bad.mProperties.mInertia(1, 1) = 1e-20f;
        bad.mProperties.mInertia(2, 2) = 1e-20f; rejected({ bad });
    }

    TEST_CASE("BodyMomentumUsesPointVelocityAndGlobalPWedgeX")
    {
        MassBody body;
        body.mID = 44;
        body.mOrigin = { 10, 20, 30, 40 };
        body.mComponents = std::vector<ComponentMass> { Component(1, 2, { 1, 0, 0, 0 }, { 1, 2, 3, 4 }) };
        body.mLinearVelocity = { 3, 4, 5, 6 };
        body.mAngularVelocity = { 2, 3, 4, 5, 6, 7 };
        BodyMassData data;
        CHECK(MassModel::Derive(body, data));
        CHECK(data.mWorldCOM == MassVector { 11, 20, 30, 40 });
        const double spin[] = { 6, 12, 20, 25, 36, 49 };
        const size_t pi[] = { 0, 0, 0, 1, 1, 2 }, pj[] = { 1, 2, 3, 2, 3, 3 };
        double energy = 0;
        for (size_t k = 0; k < 6; ++k)
        {
            CHECK_APPROX_EQUAL(data.mSpinMomentum[k], spin[k], 1e-12);
            double orbital = 2 * (body.mLinearVelocity[pi[k]] * data.mWorldCOM[pj[k]] - body.mLinearVelocity[pj[k]] * data.mWorldCOM[pi[k]]);
            CHECK_APPROX_EQUAL(data.mMomentum.mAngular[k], spin[k] + orbital, 1e-12);
            energy += 0.5 * spin[k] * body.mAngularVelocity[k];
        }
        CHECK_APPROX_EQUAL(data.mMomentum.mKineticEnergy, 86.0 + energy, 1e-12);
        MassVector pointVelocity;
        CHECK(MassModel::PointVelocity(body, { 13, 19, 30, 40 }, pointVelocity));
        CHECK(pointVelocity == MassVector { 1, 0, 4, 4 });
    }

    TEST_CASE("RecenterPreservesGeometryMomentumAndPointVelocity")
    {
        MassBody body;
        body.mID = 44;
        body.mOrigin = { 1e8, -2e8, 30, 40 };
        body.mRotation = Rotor::sExp(Bivec(0.2f, 0.3f, -0.4f, 0.5f, 0.1f, 0.7f));
        body.mComponents = std::vector<ComponentMass> { Component(1, 2, { 5, -2, 3, 1 }, { 1, 2, 3, 4 }) };
        body.mLinearVelocity = { 3, 4, 5, 6 };
        body.mAngularVelocity = { 2, 3, 4, 5, 6, 7 };
        MassBody centered;
        CHECK(MassModel::Recenter(body, centered));
        BodyMassData before, after;
        CHECK(MassModel::Derive(body, before));
        CHECK(MassModel::Derive(centered, after));
        CHECK(centered.mComponents[0].mLocalCOM == MassVector {});
        for (size_t i = 0; i < 4; ++i) CHECK_APPROX_EQUAL(before.mWorldCOM[i], after.mWorldCOM[i], 2e-8);
        CHECK(centered.mLinearVelocity == body.mLinearVelocity);
        CHECK(centered.mAngularVelocity == body.mAngularVelocity);
        MassVector v1, v2;
        CHECK(MassModel::PointVelocity(body, { 1e8 + 7, -2e8 + 1, 32, 40 }, v1));
        CHECK(MassModel::PointVelocity(centered, { 1e8 + 7, -2e8 + 1, 32, 40 }, v2));
        for (size_t i = 0; i < 4; ++i) CHECK_APPROX_EQUAL(v1[i], v2[i], 3e-7);
    }

    TEST_CASE("InternalRedistributionConservesMomentumAndAccountsWork")
    {
        MassBody body;
        body.mID = 1;
        body.mOrigin = { 10, 20, 30, 40 };
        body.mComponents = std::vector<ComponentMass> { Component(7, 2, { 0, 0, 0, 0 }, { 1, 2, 3, 4 }) };
        body.mLinearVelocity = { 3, 4, 5, 6 };
        body.mAngularVelocity = { 2, 3, 4, 5, 6, 7 };
        auto replacement = body.mComponents;
        replacement[0].mLocalCOM = { 1, 2, 3, 4 };
        replacement[0].mProperties.mInertia *= 2.0f;
        MassTransaction tx;
        CHECK(MassModel::Redistribute(body, replacement, tx));
        CHECK(tx.mBodies.size() == 1);
        if (tx.mBodies.size() != 1) return;
        BodyMassData a, b;
        CHECK(MassModel::Derive(body, a));
        CHECK(MassModel::Derive(tx.mBodies[0], b));
        CHECK(a.mWorldCOM == b.mWorldCOM);
        CHECK(body.mLinearVelocity == tx.mBodies[0].mLinearVelocity);
        for (size_t i = 0; i < 6; ++i)
        {
            CHECK_APPROX_EQUAL(tx.mBodies[0].mAngularVelocity[i], 0.5 * body.mAngularVelocity[i], 1e-12);
            CHECK_APPROX_EQUAL(a.mMomentum.mAngular[i], b.mMomentum.mAngular[i], 1e-10);
        }
        CHECK_APPROX_EQUAL(tx.mLedger.mDissipatedEnergy, 0.5 * (a.mMomentum.mKineticEnergy - 86), 1e-10);
        CHECK(tx.mLedger.mAddedEnergy == 0);
        CHECK(MassModel::Redistribute(tx.mBodies[0], body.mComponents, tx));
        CHECK(tx.mLedger.mAddedEnergy > 0);
        CHECK(tx.mLedger.mDissipatedEnergy == 0);
    }

    TEST_CASE("FreeFragmentationInheritsRigidMotionWithoutDoubleMass")
    {
        MassBody body;
        body.mID = 1;
        body.mOrigin = { 10, 20, -30, 40 };
        body.mRotation = Rotor::sExp(Bivec(0.2f, 0.3f, 0.1f, -0.1f, 0.4f, 0.5f));
        body.mComponents = std::vector<ComponentMass> { Component(7, 2, { -2, 1, 0, 3 }, { 1, 2, 3, 4 }), Component(8, 6, { 2, -1, 2, 0 }, { 2, 3, 4, 5 }) };
        body.mLinearVelocity = { 3, 4, -5, 6 };
        body.mAngularVelocity = { 2, 3, -4, 5, 6, -7 };
        MassTransaction tx;
        CHECK(MassModel::Fragment(body, { { 11, { 7 } }, { 12, { 8 } } }, tx));
        CHECK(tx.mBodies.size() == 2);
        BodyMassData initial;
        CHECK(MassModel::Derive(body, initial));
        MassVector p {};
        MassBivector l {};
        double energy = 0;
        for (const auto &piece : tx.mBodies)
        {
            BodyMassData d;
            CHECK(MassModel::Derive(piece, d));
            MassVector v;
            CHECK(MassModel::PointVelocity(body, d.mWorldCOM, v));
            for (size_t i = 0; i < 4; ++i)
            {
                CHECK_APPROX_EQUAL(piece.mLinearVelocity[i], v[i], 1e-12);
                p[i] += d.mMomentum.mLinear[i];
            }
            CHECK(piece.mAngularVelocity == body.mAngularVelocity);
            CHECK(piece.mOrigin == body.mOrigin);
            CHECK(piece.mRotation == body.mRotation);
            for (size_t i = 0; i < 6; ++i) l[i] += d.mMomentum.mAngular[i];
            energy += d.mMomentum.mKineticEnergy;
        }
        for (size_t i = 0; i < 4; ++i) CHECK_APPROX_EQUAL(p[i], initial.mMomentum.mLinear[i], 1e-10);
        for (size_t i = 0; i < 6; ++i) CHECK_APPROX_EQUAL(l[i], initial.mMomentum.mAngular[i], 1e-9);
        CHECK_APPROX_EQUAL(energy, initial.mMomentum.mKineticEnergy, 1e-9);
        CHECK_FALSE(MassModel::Fragment(body, { { 11, { 7 } }, { 12, { 7, 8 } } }, tx));
        CHECK(tx.mBodies.size() == 2);
        CHECK_FALSE(MassModel::Fragment(body, { { 11, { 7 } }, { 11, { 8 } } }, tx));
        CHECK_FALSE(MassModel::Fragment(body, { { 11, { 7 } }, { 12, { 999 } } }, tx));
    }

    TEST_CASE("RigidDockingConservesOrbitalMomentumAndDissipatesEnergy")
    {
        MassBody a, b;
        a.mID = 1; b.mID = 2;
        a.mComponents = std::vector<ComponentMass> { Component(11, 2, {}, { 1, 1, 1, 1 }) };
        b.mComponents = std::vector<ComponentMass> { Component(12, 2, {}, { 1, 1, 1, 1 }) };
        a.mOrigin = { -1, 0, 0, 0 }; b.mOrigin = { 1, 0, 0, 0 };
        a.mLinearVelocity = { 0, 1, 0, 0 }; b.mLinearVelocity = { 0, -1, 0, 0 };
        MassTransaction tx;
        CHECK(MassModel::Dock({ a, b }, 90, tx));
        CHECK(tx.mBodies.size() == 1);
        if (tx.mBodies.size() != 1) return;
        const auto &merged = tx.mBodies[0];
        CHECK(merged.mID == 90);
        CHECK(merged.mComponents.size() == 2);
        CHECK(merged.mLinearVelocity == MassVector {});
        CHECK_APPROX_EQUAL(merged.mAngularVelocity[0], 0.5, 1e-12);
        for (size_t i = 1; i < 6; ++i) CHECK_APPROX_EQUAL(merged.mAngularVelocity[i], 0.0, 1e-12);
        CHECK_APPROX_EQUAL(tx.mLedger.mBefore.mKineticEnergy, 2.0, 1e-12);
        CHECK_APPROX_EQUAL(tx.mLedger.mAfter.mKineticEnergy, 1.0, 1e-12);
        CHECK_APPROX_EQUAL(tx.mLedger.mDissipatedEnergy, 1.0, 1e-12);
        CHECK(tx.mLedger.mAddedEnergy == 0);
        MassTransaction reverse;
        CHECK(MassModel::Dock({ b, a }, 90, reverse));
        CHECK(reverse.mBodies[0].mAngularVelocity == merged.mAngularVelocity);
        CHECK(reverse.mBodies[0].mOrigin == merged.mOrigin);
        b.mComponents[0].mID = 11;
        CHECK_FALSE(MassModel::Dock({ a, b }, 90, tx));
        CHECK(tx.mBodies[0].mID == 90);
    }

    TEST_CASE("LumpedExhaustTransfersLinearAndOffCenterAngularMomentum")
    {
        MassBody body;
        body.mID = 1;
        body.mComponents = std::vector<ComponentMass> { Component(7, 4, {}, { 1, 1, 1, 1 }), Component(8, 2, { 0, 1, 0, 0 }, {}) };
        LumpedExhaust gas;
        gas.mID = 100;
        gas.mMass = 2;
        gas.mWorldPosition = { 0, 1, 0, 0 };
        gas.mWorldVelocity = { -3, 0, 0, 0 };
        MassTransaction tx;
        CHECK(MassModel::Expel(body, { body.mComponents[0] }, { gas }, tx));
        CHECK(tx.mBodies.size() == 1);
        if (tx.mBodies.size() != 1) return;
        const auto &remaining = tx.mBodies[0];
        CHECK_APPROX_EQUAL(remaining.mLinearVelocity[0], 1.5, 1e-12);
        CHECK_APPROX_EQUAL(remaining.mAngularVelocity[0], 3.0, 1e-12);
        BodyMassData data;
        CHECK(MassModel::Derive(remaining, data));
        CHECK_APPROX_EQUAL(data.mMomentum.mAngular[0], 6.0, 1e-12);
        CHECK_APPROX_EQUAL(tx.mLedger.mExported.mMass, 2.0, 1e-12);
        CHECK_APPROX_EQUAL(tx.mLedger.mAddedEnergy, 22.5, 1e-12);
        CHECK(tx.mLedger.mDissipatedEnergy == 0);
        CHECK(tx.mExhaust.size() == 1);
        CHECK(tx.mExhaust[0].mID == gas.mID);
        gas.mMass = -1;
        CHECK_FALSE(MassModel::Expel(body, { body.mComponents[0] }, { gas }, tx));
        CHECK(tx.mBodies[0].mComponents[0].mID == 7);
        gas.mMass = 2; gas.mWorldPosition[1] = 100;
        CHECK_FALSE(MassModel::Expel(body, { body.mComponents[0] }, { gas }, tx));
    }

    TEST_CASE("RandomPointCloudReferenceValidatesCoupledInertiaAndTopology")
    {
        UnitTestRandom random(0x4DCA55);
        std::uniform_real_distribution<float> u(-1.0f, 1.0f);
        for (int trial = 0; trial < 100; ++trial)
        {
            CAPTURE(trial);
            auto rotor = [&]() { return Cl4Ref::Exp(Cl4Ref::FromBivec(Bivec(u(random), u(random), u(random), u(random), u(random), u(random)))); };
            auto bodyRotation = rotor();
            MassBody body;
            body.mID = 1;
            body.mRotation = Cl4Ref::ToRotor(bodyRotation);
            body.mOrigin = { 10*u(random), 10*u(random), 10*u(random), 10*u(random) };
            body.mLinearVelocity = { 3*u(random), 3*u(random), 3*u(random), 3*u(random) };
            for (double &w : body.mAngularVelocity) w = 2*u(random);
            std::vector<std::pair<double, MassVector>> particles;
            for (uint64 id = 1; id <= 4; ++id)
            {
                MassVector center { 2*u(random), 2*u(random), 2*u(random), 2*u(random) };
                auto c = Component(id, float(id + 1), center, { 1.0 + 0.5*u(random), 1.0 + 0.5*u(random), 1.0 + 0.5*u(random), 1.0 + 0.5*u(random) });
                auto rotation = rotor();
                c.mLocalRotation = Cl4Ref::ToRotor(rotation);
                body.mComponents.push_back(c);
                // Eight equally weighted axial point masses reproduce each diagonal tensor.
                for (size_t axis = 0; axis < 4; ++axis)
                    for (double sign : { -1.0, 1.0 })
                    {
                        MassVector x {};
                        x[axis] = sign * std::sqrt(4 * double(c.mProperties.mInertia(uint(axis), uint(axis))) / c.mProperties.mMass);
                        x = ReferenceRotate(rotation, x);
                        for (size_t i = 0; i < 4; ++i) x[i] += center[i];
                        x = ReferenceRotate(bodyRotation, x);
                        for (size_t i = 0; i < 4; ++i) x[i] += body.mOrigin[i];
                        particles.push_back({ double(c.mProperties.mMass) / 8, x });
                    }
            }
            BodyMassData data;
            MASS_REQUIRE(MassModel::Derive(body, data));
            double mass = 0;
            MassVector com {};
            for (const auto &p : particles)
            {
                mass += p.first;
                for (size_t i = 0; i < 4; ++i) com[i] += p.first*p.second[i];
            }
            for (double &x : com) x /= mass;
            MassMatrix second {};
            MassBivector angular {};
            double energy = 0;
            for (const auto &p : particles)
            {
                auto v = ReferencePointVelocity(body, p.second, com);
                for (size_t i = 0; i < 4; ++i)
                {
                    energy += 0.5*p.first*v[i]*v[i];
                    for (size_t j = 0; j < 4; ++j) second[i][j] += p.first*(p.second[i] - com[i])*(p.second[j] - com[j]);
                }
                size_t k = 0;
                for (size_t i = 0; i < 4; ++i)
                    for (size_t j = i + 1; j < 4; ++j) angular[k++] += p.first*(v[i]*p.second[j] - v[j]*p.second[i]);
            }
            for (size_t i = 0; i < 4; ++i)
            {
                CHECK_APPROX_EQUAL(data.mWorldCOM[i], com[i], 2e-6);
                for (size_t j = 0; j < 4; ++j) CHECK_APPROX_EQUAL(data.mWorldSecondMoment[i][j], second[i][j], 2e-5);
            }
            for (size_t i = 0; i < 6; ++i) CHECK_APPROX_EQUAL(data.mMomentum.mAngular[i], angular[i], 5e-5);
            CHECK_APPROX_EQUAL(data.mMomentum.mKineticEnergy, energy, 1e-4);

            MassTransaction split, docked, redistributed;
            MASS_REQUIRE(MassModel::Fragment(body, { { 11, { 1, 3 } }, { 12, { 2, 4 } } }, split));
            CheckGlobalConservation({ body }, split);
            MASS_REQUIRE(MassModel::Dock(split.mBodies, 30, docked));
            CheckGlobalConservation(split.mBodies, docked);
            for (size_t i = 0; i < 6; ++i) CHECK_APPROX_EQUAL(docked.mBodies[0].mAngularVelocity[i], body.mAngularVelocity[i], 3e-6);
            auto replacement = body.mComponents;
            for (auto &c : replacement)
            {
                for (double &x : c.mLocalCOM) x += u(random);
                c.mProperties.mInertia *= 1.5f;
            }
            MASS_REQUIRE(MassModel::Redistribute(body, replacement, redistributed));
            CheckGlobalConservation({ body }, redistributed);
        }
    }

    TEST_CASE("LargeOffsetAggregationIsCenteredAndStableByID")
    {
        auto a = Component(0x100000001ULL, 3, { 1e8 - 0.25, -1e8, 0, 0 }, { 1, 2, 3, 4 });
        auto b = Component(0x200000001ULL, 1, { 1e8 + 0.75, -1e8, 0, 0 }, { 2, 3, 4, 5 });
        AggregateMass x, y;
        MASS_REQUIRE(MassModel::Aggregate({ a, b }, x));
        MASS_REQUIRE(MassModel::Aggregate({ b, a }, y));
        CHECK(x.mLocalCOM == y.mLocalCOM);
        CHECK(x.mSecondMoment == y.mSecondMoment);
        CHECK(x.mLocalCOM[0] == 1e8);
        CHECK_APPROX_EQUAL(x.mSecondMoment[0][0], 3.75, 1e-12);
    }

    TEST_CASE("RandomLumpedExhaustConservesGlobalMomentum")
    {
        UnitTestRandom random(0xE4A057);
        std::uniform_real_distribution<float> u(-1.0f, 1.0f);
        for (int trial = 0; trial < 100; ++trial)
        {
            MassBody body;
            body.mID = 1;
            body.mOrigin = { 10*u(random), 10*u(random), 10*u(random), 10*u(random) };
            body.mComponents = std::vector<ComponentMass> { Component(1, 10, {}, { 2, 3, 4, 5 }) };
            body.mLinearVelocity = { 3*u(random), 3*u(random), 3*u(random), 3*u(random) };
            for (double &x : body.mAngularVelocity) x = u(random);
            std::vector<LumpedExhaust> exhaust;
            for (uint64 id = 2; id < 6; ++id)
            {
                auto c = Component(id, 1, { u(random), u(random), u(random), u(random) }, {});
                body.mComponents.push_back(c);
                LumpedExhaust gas;
                gas.mID = id; gas.mMass = 1;
                for (size_t i = 0; i < 4; ++i)
                {
                    gas.mWorldPosition[i] = body.mOrigin[i] + c.mLocalCOM[i];
                    gas.mWorldVelocity[i] = 5*u(random);
                }
                exhaust.push_back(gas);
            }
            MassTransaction tx;
            MASS_REQUIRE(MassModel::Expel(body, { body.mComponents[0] }, exhaust, tx));
            CheckGlobalConservation({ body }, tx);
        }
    }

    TEST_CASE("CheckedNativeMassPropertiesConversionDoesNotPublishDegenerateRounding")
    {
        AggregateMass aggregate;
        MASS_REQUIRE(MassModel::Aggregate({ Component(1, 2, {}, { 1, 2, 3, 4 }) }, aggregate));
        MassProperties native;
        MASS_REQUIRE(MassModel::ToMassProperties(aggregate, native));
        CHECK(native.mMass == 2);
        for (uint i = 0; i < 4; ++i) CHECK(native.mInertia(i, i) == float(i + 1));
        aggregate.mSecondMoment = {};
        for (size_t i = 0; i < 4; ++i) aggregate.mSecondMoment[i][i] = 1e-100;
        CHECK_FALSE(MassModel::ToMassProperties(aggregate, native));
        CHECK(native.mMass == 2);
        CHECK(native.mInertia(3, 3) == 4);
        aggregate.mMass = std::numeric_limits<double>::infinity();
        CHECK_FALSE(MassModel::ToMassProperties(aggregate, native));
    }

    TEST_CASE("AllPrepareFailuresPreserveOutputAndRejectBrokenOwnership")
    {
        MassBody body;
        body.mID = 1;
        body.mComponents.push_back(Component(1, 4, {}, { 1, 1, 1, 1 }));
        body.mComponents.push_back(Component(2, 2, {}, { 1, 1, 1, 1 }));
        MassTransaction out;
        out.mBodies.push_back(body);
        out.mLedger.mAddedEnergy = 123;
        auto unchanged = [&]() { CHECK(out.mBodies.size() == 1); CHECK(out.mBodies[0].mID == 1); CHECK(out.mLedger.mAddedEnergy == 123); };
        auto replacement = body.mComponents;
        replacement[0].mID = 999;
        CHECK_FALSE(MassModel::Redistribute(body, replacement, out)); unchanged();
        replacement = body.mComponents; replacement[0].mProperties.mMass = 5;
        CHECK_FALSE(MassModel::Redistribute(body, replacement, out)); unchanged();
        CHECK_FALSE(MassModel::Fragment(body, { { 1, { 1 } }, { 2, {} } }, out)); unchanged();
        CHECK_FALSE(MassModel::Fragment(body, { { 1, { 1 } }, { 0, { 2 } } }, out)); unchanged();
        CHECK_FALSE(MassModel::Dock({ body, body }, 5, out)); unchanged();
        CHECK_FALSE(MassModel::Dock({ body }, 5, out)); unchanged();
        CHECK_FALSE(MassModel::Dock({ body, body }, 0, out)); unchanged();
        LumpedExhaust gas;
        gas.mID = 3; gas.mMass = 2;
        CHECK_FALSE(MassModel::Expel(body, { body.mComponents[0] }, { gas, gas }, out)); unchanged();
        gas.mMass = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(MassModel::Expel(body, { body.mComponents[0] }, { gas }, out)); unchanged();
        MassBody sentinel = body;
        auto bad = body; bad.mAngularVelocity[5] = std::numeric_limits<double>::infinity();
        CHECK_FALSE(MassModel::Recenter(bad, sentinel));
        CHECK(sentinel.mID == 1);
        CHECK(sentinel.mComponents.size() == 2);
        BodyMassData data; data.mMomentum.mMass = 123;
        CHECK_FALSE(MassModel::Derive(bad, data));
        CHECK(data.mMomentum.mMass == 123);
        MassVector v { 1, 2, 3, 4 };
        CHECK_FALSE(MassModel::PointVelocity(bad, {}, v));
        CHECK(v == MassVector { 1, 2, 3, 4 });
    }

    TEST_CASE("DockingEnergyLedgerResolvesSmallLossAtLargeCommonVelocity")
    {
        MassBody a, b;
        a.mID = 1; b.mID = 2;
        a.mComponents.push_back(Component(1, 2, {}, { 1, 1, 1, 1 }));
        b.mComponents.push_back(Component(2, 2, {}, { 1, 1, 1, 1 }));
        a.mOrigin = { -1, 0, 0, 0 }; b.mOrigin = { 1, 0, 0, 0 };
        a.mLinearVelocity = { 1e8 + 1, 0, 0, 0 }; b.mLinearVelocity = { 1e8 - 1, 0, 0, 0 };
        MassTransaction tx;
        MASS_REQUIRE(MassModel::Dock({ a, b }, 3, tx));
        CHECK_APPROX_EQUAL(tx.mLedger.mDissipatedEnergy, 2.0, 1e-12);
        CHECK(tx.mLedger.mAddedEnergy == 0);
    }

    TEST_CASE("ComovingUnequalMassDockingDoesNotInventSpinUnderBoost")
    {
        MassBody a, b;
        a.mID = 1; b.mID = 2;
        a.mComponents.push_back(Component(1, 2, {}, { 1, 1, 1, 1 }));
        b.mComponents.push_back(Component(2, 3, {}, { 1, 1, 1, 1 }));
        b.mOrigin = { 0, 0.1, 0, 0 };
        MassTransaction rest, boosted;
        MASS_REQUIRE(MassModel::Dock({ a, b }, 3, rest));
        a.mLinearVelocity = b.mLinearVelocity = { 1e8, 0, 0, 0 };
        MASS_REQUIRE(MassModel::Dock({ a, b }, 3, boosted));
        CHECK(boosted.mBodies[0].mLinearVelocity == a.mLinearVelocity);
        CHECK(boosted.mBodies[0].mAngularVelocity == MassBivector {});
        CHECK(boosted.mLedger.mAddedEnergy == 0);
        CHECK(boosted.mLedger.mDissipatedEnergy == rest.mLedger.mDissipatedEnergy);
    }

    TEST_CASE("NonbinaryOffsetFragmentationSurvivesAbsoluteVelocityRounding")
    {
        MassBody body;
        body.mID = 1;
        body.mComponents = { Component(1, 2, {}, { 1, 1, 1, 1 }), Component(2, 3, { 0, 0.3, 0, 0 }, { 1, 1, 1, 1 }) };
        body.mAngularVelocity[0] = 1;
        MassTransaction rest, boosted;
        MASS_REQUIRE(MassModel::Fragment(body, { { 2, { 1 } }, { 3, { 2 } } }, rest));
        body.mLinearVelocity = { 1e8, 0, 0, 0 };
        MASS_REQUIRE(MassModel::Fragment(body, { { 2, { 1 } }, { 3, { 2 } } }, boosted));
        CHECK(boosted.mLedger.mEnergyFrameVelocity == body.mLinearVelocity);
        // Output is an absolute double velocity: compare to the correctly rounded sum,
        // not to an unrealistically exact subtraction of the common boost.
        CHECK(boosted.mBodies[0].mLinearVelocity[0] == 1e8 - 0.18);
        CHECK(boosted.mBodies[1].mLinearVelocity[0] == 1e8 + 0.12);
        for (const auto &piece : boosted.mBodies) CHECK(piece.mAngularVelocity == body.mAngularVelocity);
        CHECK_APPROX_EQUAL(boosted.mLedger.mBefore.mKineticEnergy, 2.054, 1e-14);
        CHECK(rest.mLedger.mVelocityRoundingEnergyBound == 0);
        CHECK(boosted.mLedger.mVelocityRoundingEnergyBound > 0);
        CHECK(boosted.mLedger.mVelocityRoundingEnergyBound < 1e-8);
        double error = boosted.mLedger.mAddedEnergy + boosted.mLedger.mDissipatedEnergy;
        CHECK(error > 1e-9 * boosted.mLedger.mBefore.mKineticEnergy); // Old blanket relative test rejects.
        CHECK(error <= boosted.mLedger.mVelocityRoundingEnergyBound + 1e-14);
    }

    TEST_CASE("ComovingPointExpulsionDoesNotInventRetainedSpinUnderBoost")
    {
        MassBody body;
        body.mID = 1;
        body.mComponents = { Component(1, 2, {}, { 1, 1, 1, 1 }), Component(2, 3, { 0, 0.7, 0, 0 }, {}) };
        LumpedExhaust gas;
        gas.mID = 2; gas.mMass = 3; gas.mWorldPosition = { 0, 0.7, 0, 0 };
        MassTransaction rest, boosted;
        MASS_REQUIRE(MassModel::Expel(body, { body.mComponents[0] }, { gas }, rest));
        body.mLinearVelocity = gas.mWorldVelocity = { 1e8, 0, 0, 0 };
        MASS_REQUIRE(MassModel::Expel(body, { body.mComponents[0] }, { gas }, boosted));
        CHECK(boosted.mBodies[0].mLinearVelocity == body.mLinearVelocity);
        CHECK(boosted.mBodies[0].mAngularVelocity == rest.mBodies[0].mAngularVelocity);
        CHECK(boosted.mLedger.mAddedEnergy == 0);
        CHECK(boosted.mLedger.mDissipatedEnergy == 0);
    }

    TEST_CASE("SeededGeneralCoordinatesAndBoostsPreserveCentralTopologyMotion")
    {
        UnitTestRandom random(0xB00574D);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        for (int trial = 0; trial < 100; ++trial)
        {
            CAPTURE(trial);
            MassVector boost;
            for (double &x : boost) x = 1e8 * u(random);
            MassBody body;
            body.mID = 1;
            body.mRotation = Rotor::sExp(Bivec(float(u(random)), float(u(random)), float(u(random)), float(u(random)), float(u(random)), float(u(random))));
            for (double &x : body.mOrigin) x = 10 * u(random);
            for (double &x : body.mAngularVelocity) x = 2 * u(random);
            for (uint64 id = 1; id <= 4; ++id)
                body.mComponents.push_back(Component(id, float(id + 1), { 2*u(random), 2*u(random), 2*u(random), 2*u(random) }, { 1, 2, 3, 4 }));
            MassTransaction rest, boosted;
            MASS_REQUIRE(MassModel::Fragment(body, { { 11, { 1, 3 } }, { 12, { 2, 4 } } }, rest));
            body.mLinearVelocity = boost;
            MASS_REQUIRE(MassModel::Fragment(body, { { 11, { 1, 3 } }, { 12, { 2, 4 } } }, boosted));
            // Independent long-double energy perturbation bound, using the zero-speed
            // outputs as the ideal local velocities, not the production TwoSum residue.
            long double bound = 0;
            for (size_t n = 0; n < rest.mBodies.size(); ++n)
            {
                BodyMassData data;
                MASS_REQUIRE(MassModel::Derive(rest.mBodies[n], data));
                CHECK(boosted.mBodies[n].mAngularVelocity == rest.mBodies[n].mAngularVelocity);
                for (size_t i = 0; i < 4; ++i)
                {
                    double relative = rest.mBodies[n].mLinearVelocity[i];
                    double stored = boosted.mBodies[n].mLinearVelocity[i];
                    CHECK(stored == boost[i] + relative);
                    long double e = (static_cast<long double>(stored) - boost[i]) - relative;
                    bound += data.mLocal.mMass * (std::abs(static_cast<long double>(relative) * e) + 0.5L * e * e);
                }
            }
            CHECK_APPROX_EQUAL(boosted.mLedger.mVelocityRoundingEnergyBound, double(bound), 1e-20);
            CHECK(boosted.mLedger.mAddedEnergy + boosted.mLedger.mDissipatedEnergy
                <= double(bound) + 1e-9 * boosted.mLedger.mBefore.mKineticEnergy);
            CHECK(rest.mLedger.mVelocityRoundingEnergyBound == 0);

            // A boost quantizes input velocities too. Express the SAME represented inputs
            // in a nearby frame; central spin should be bit-identical, not merely close.
            auto nearby = boosted.mBodies;
            const MassVector anchor = nearby.front().mLinearVelocity;
            for (auto &piece : nearby)
                for (size_t i = 0; i < 4; ++i) piece.mLinearVelocity[i] -= anchor[i];
            MassTransaction docked, nearbyDocked;
            MASS_REQUIRE(MassModel::Dock(boosted.mBodies, 30, docked));
            MASS_REQUIRE(MassModel::Dock(nearby, 30, nearbyDocked));
            CHECK(docked.mBodies[0].mAngularVelocity == nearbyDocked.mBodies[0].mAngularVelocity);
            for (size_t i = 0; i < 4; ++i)
                CHECK(docked.mBodies[0].mLinearVelocity[i] == anchor[i] + nearbyDocked.mBodies[0].mLinearVelocity[i]);

            // Co-moving unequal masses and arbitrary four-lane offsets must never
            // generate spin or kinetic work, regardless of first-moment roundoff.
            for (auto &piece : nearby)
            {
                piece.mLinearVelocity = boost;
                piece.mAngularVelocity = {};
            }
            MASS_REQUIRE(MassModel::Dock(nearby, 30, docked));
            CHECK(docked.mBodies[0].mAngularVelocity == MassBivector {});
            CHECK(docked.mLedger.mAddedEnergy == 0);
            CHECK(docked.mLedger.mDissipatedEnergy == 0);
        }
    }

    TEST_CASE("SeededPointExpulsionUsesTheSameStableVelocityReferenceUnderBoost")
    {
        UnitTestRandom random(0xE4B0057);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        for (int trial = 0; trial < 100; ++trial)
        {
            CAPTURE(trial);
            MassBody body;
            body.mID = 1;
            for (double &x : body.mOrigin) x = 10 * u(random);
            for (double &x : body.mLinearVelocity) x = 1e8 * u(random);
            for (double &x : body.mAngularVelocity) x = u(random);
            body.mComponents.push_back(Component(1, 2, {}, { 1, 2, 3, 4 }));
            std::vector<LumpedExhaust> exhaust;
            for (uint64 id = 2; id <= 4; ++id)
            {
                auto c = Component(id, float(id + 1), { u(random), u(random), u(random), u(random) }, {});
                body.mComponents.push_back(c);
                LumpedExhaust gas;
                gas.mID = id; gas.mMass = c.mProperties.mMass;
                for (size_t i = 0; i < 4; ++i)
                {
                    gas.mWorldPosition[i] = body.mOrigin[i] + c.mLocalCOM[i];
                    gas.mWorldVelocity[i] = body.mLinearVelocity[i] + 3*u(random);
                }
                exhaust.push_back(gas);
            }
            auto nearbyBody = body;
            nearbyBody.mLinearVelocity = {};
            auto nearbyExhaust = exhaust;
            for (auto &gas : nearbyExhaust)
                for (size_t i = 0; i < 4; ++i) gas.mWorldVelocity[i] -= body.mLinearVelocity[i];
            MassTransaction boosted, nearby;
            MASS_REQUIRE(MassModel::Expel(body, { body.mComponents[0] }, exhaust, boosted));
            MASS_REQUIRE(MassModel::Expel(nearbyBody, { body.mComponents[0] }, nearbyExhaust, nearby));
            CHECK(boosted.mBodies[0].mAngularVelocity == nearby.mBodies[0].mAngularVelocity);
            for (size_t i = 0; i < 4; ++i)
                CHECK(boosted.mBodies[0].mLinearVelocity[i] == body.mLinearVelocity[i] + nearby.mBodies[0].mLinearVelocity[i]);
            body.mAngularVelocity = {};
            for (auto &gas : exhaust) gas.mWorldVelocity = body.mLinearVelocity;
            MASS_REQUIRE(MassModel::Expel(body, { body.mComponents[0] }, exhaust, boosted));
            CHECK(boosted.mBodies[0].mAngularVelocity == MassBivector {});
            CHECK(boosted.mBodies[0].mLinearVelocity == body.mLinearVelocity);
            CHECK(boosted.mLedger.mAddedEnergy == 0);
            CHECK(boosted.mLedger.mDissipatedEnergy == 0);
        }
    }

    TEST_CASE("SymmetricIndefiniteTensorsAndSingularFragmentsReject")
    {
        auto c = Component(1, 2, {}, { 2, 2, 2, 2 });
        c.mProperties.mInertia(0, 1) = c.mProperties.mInertia(1, 0) = 1;
        AggregateMass out;
        MASS_REQUIRE(MassModel::Aggregate({ c }, out));
        CHECK_APPROX_EQUAL(out.mSecondMoment[0][1], 1.0, 1e-12);
        c.mProperties.mInertia(0, 1) = c.mProperties.mInertia(1, 0) = 3;
        CHECK_FALSE(MassModel::Aggregate({ c }, out));
        CHECK(out.mSecondMoment[0][1] == 1);
        auto point = Component(2, 1, { 1, 0, 0, 0 }, {});
        MassBody body;
        body.mID = 1;
        body.mComponents.push_back(Component(1, 2, {}, { 2, 2, 2, 2 }));
        body.mComponents.push_back(point);
        MassTransaction tx;
        tx.mLedger.mAddedEnergy = 123;
        CHECK_FALSE(MassModel::Fragment(body, { { 1, { 1 } }, { 2, { 2 } } }, tx));
        CHECK(tx.mBodies.empty());
        CHECK(tx.mLedger.mAddedEnergy == 123);
        std::vector<ComponentMass> tooMany(4097, body.mComponents[0]);
        CHECK_FALSE(MassModel::Aggregate(tooMany, out));
    }

    TEST_CASE("HighMassRatioDockingRetainsSmallBodyMomentum")
    {
        MassBody a, b;
        a.mID = 1; b.mID = 2;
        a.mComponents.push_back(Component(1, 1e6f, {}, { 1e6, 1e6, 1e6, 1e6 }));
        b.mComponents.push_back(Component(2, 1, {}, { 1, 1, 1, 1 }));
        a.mOrigin = { -1, 0, 0, 0 }; b.mOrigin = { 1, 0, 0, 0 };
        b.mLinearVelocity = { 3, 0, 0, 0 };
        MassTransaction tx;
        MASS_REQUIRE(MassModel::Dock({ a, b }, 3, tx));
        CHECK_APPROX_EQUAL(tx.mBodies[0].mLinearVelocity[0], 3.0 / 1000001.0, 1e-18);
        CHECK_APPROX_EQUAL(tx.mLedger.mDissipatedEnergy, 4.5 * (1000000.0 / 1000001.0), 1e-12);
        CheckGlobalConservation({ a, b }, tx);
    }

    TEST_CASE("DockingPreservesWorldComponentGeometryUnderGeneralRotations")
    {
        MassBody a, b;
        a.mID = 1; b.mID = 2;
        a.mComponents.push_back(Component(1, 1, { 1, 2, 3, 4 }, { 1, 2, 3, 4 }));
        b.mComponents.push_back(Component(2, 2, { -2, 1, 4, 3 }, { 2, 1, 4, 3 }));
        a.mRotation = Rotor::sExp(Bivec(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f));
        b.mRotation = Rotor::sExp(Bivec(-0.3f, 0.7f, 0.4f, -0.1f, 0.2f, -0.8f));
        a.mComponents[0].mLocalRotation = b.mRotation;
        b.mComponents[0].mLocalRotation = a.mRotation;
        a.mOrigin = { -3, 5, 10, 2 }; b.mOrigin = { 7, 2, 6, -4 };
        MassTransaction tx;
        MASS_REQUIRE(MassModel::Dock({ a, b }, 3, tx));
        auto &merged = tx.mBodies[0];
        for (const MassBody *original : { &a, &b })
        {
            const auto &c = original->mComponents[0];
            auto found = std::find_if(merged.mComponents.begin(), merged.mComponents.end(), [&](const ComponentMass &x) { return c.mID == x.mID; });
            MASS_REQUIRE(found != merged.mComponents.end());
            auto before = ReferenceRotate(Cl4Ref::FromRotor(c.mLocalRotation), { 0.25, -0.5, 0.75, 1.0 });
            for (size_t i = 0; i < 4; ++i) before[i] += c.mLocalCOM[i];
            before = ReferenceRotate(Cl4Ref::FromRotor(original->mRotation), before);
            auto after = ReferenceRotate(Cl4Ref::FromRotor(found->mLocalRotation), { 0.25, -0.5, 0.75, 1.0 });
            for (size_t i = 0; i < 4; ++i)
            {
                before[i] += original->mOrigin[i];
                after[i] += found->mLocalCOM[i] + merged.mOrigin[i];
                CHECK_APPROX_EQUAL(before[i], after[i], 3e-6);
            }
        }
    }

    TEST_CASE("AggregateAnalyticalRotatedComponents")
    {
        auto a = Component(7, 2, { -2, 0, 0, 0 }, { 1, 2, 3, 4 });
        auto b = Component(19, 6, { 2, 0, 0, 0 }, { 3, 4, 5, 6 });
        b.mLocalRotation = Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), 0.5f * JPH_PI);
        AggregateMass out;
        CHECK(MassModel::Aggregate({ a, b }, out));
        CHECK_APPROX_EQUAL(out.mMass, 8.0, 1e-12);
        CHECK_APPROX_EQUAL(out.mLocalCOM[0], 1.0, 1e-12);
        const double expected[] = { 29, 5, 8, 10 };
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                CHECK_APPROX_EQUAL(out.mSecondMoment[size_t(i)][size_t(j)], i == j? expected[i] : 0.0, 2e-6);
    }
}
