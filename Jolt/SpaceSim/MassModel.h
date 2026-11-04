// SPDX-License-Identifier: MIT
#pragma once

#include <Jolt/Math/Rotor.h>
#include <Jolt/Physics/Body/MassProperties.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

JPH_NAMESPACE_BEGIN
namespace SpaceSim
{
using MassVector = std::array<double, 4>;
using MassBivector = std::array<double, 6>;
using MassMatrix = std::array<MassVector, 4>;

/// Component coordinates are relative to a body geometry origin, NOT its aggregate COM.
/// mProperties contains central second moments in the component's own rotated frame.
struct ComponentMass
{
    uint64 mID = 0;
    MassVector mLocalCOM {};
    Rotor mLocalRotation = Rotor::sIdentity();
    MassProperties mProperties;
};

struct AggregateMass
{
    double mMass = 0;
    MassVector mLocalCOM {};
    MassMatrix mSecondMoment {};
};

/// Input/output velocity snapshot, not an additional authoritative momentum store.
/// The physics adapter must atomically install velocities, transformed collision geometry and mass.
struct MassBody
{
    uint64 mID = 0;
    std::vector<ComponentMass> mComponents;
    MassVector mOrigin {}; ///< World position of the fixed geometry origin (not COM).
    Rotor mRotation = Rotor::sIdentity();
    MassVector mLinearVelocity {}; ///< World COM velocity.
    MassBivector mAngularVelocity {}; ///< World bivector [XY,XZ,XW,YZ,YW,ZW].
};

/// Derived diagnostics only; never authoritative alongside MassBody velocities.
struct MassMomentum
{
    double mMass = 0;
    MassVector mLinear {};
    MassBivector mAngular {}; ///< Spin + p wedge (COM - reference origin).
    double mKineticEnergy = 0;
};

struct BodyMassData
{
    AggregateMass mLocal;
    MassVector mWorldCOM {};
    MassMatrix mWorldSecondMoment {};
    MassBivector mSpinMomentum {};
    MassMomentum mMomentum;
};

struct FragmentRequest
{
    uint64 mBodyID = 0;
    std::vector<uint64> mComponentIDs;
};

struct MassLedger
{
    MassVector mReferenceOrigin {}; ///< Shared world reference for all angular entries.
    MassVector mEnergyFrameVelocity {}; ///< Kinetic energies ONLY are measured in this inertial frame.
    MassMomentum mBefore, mAfter, mExported, mResidual;
    double mAddedEnergy = 0; ///< Positive net kinetic work required by this model.
    double mDissipatedEnergy = 0; ///< Positive kinetic energy removed (heat/deformation etc.).
    double mVelocityRoundingEnergyBound = 0; ///< Fragment-only bound from storing absolute double velocities; not work.
};

/// Point-mass flux sample at the transaction epoch. World velocity is ABSOLUTE, not nozzle
/// relative speed. Intrinsic exhaust spin / unresolved thermal energy are outside this model.
struct LumpedExhaust
{
    uint64 mID = 0;
    double mMass = 0;
    MassVector mWorldPosition {}, mWorldVelocity {};
};

/// Fully prepared value result. No body, physics world, or caller output is mutated on failure.
struct MassTransaction
{
    std::vector<MassBody> mBodies;
    std::vector<LumpedExhaust> mExhaust;
    MassLedger mLedger;
};

namespace MassDetail
{
struct Sum
{
    double mValue = 0, mCorrection = 0;
    void Add(double x)
    {
        double next = mValue + x;
        mCorrection += std::abs(mValue) >= std::abs(x)? (mValue - next) + x : (x - next) + mValue;
        mValue = next;
    }
    double Get() const { return mValue + mCorrection; }
};

inline bool Bounded(const MassVector &v, double limit = 1e9)
{
    for (double x : v) if (!std::isfinite(x) || std::abs(x) > limit) return false;
    return true;
}

/// Input float tensors may have roundoff asymmetry / tiny negative eigenvalues only.
inline bool Tensor(MassMatrix &m)
{
    double scale = 0;
    for (const auto &row : m)
        for (double x : row)
        {
            if (!std::isfinite(x) || std::abs(x) > 1e36) return false;
            scale = std::max(scale, std::abs(x));
        }
    if (scale == 0) return true; // Point components are allowed; aggregate must have full inertia.
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = i + 1; j < 4; ++j)
        {
            if (std::abs(m[i][j] - m[j][i]) > 1e-6 * scale) return false;
            m[i][j] = m[j][i] = 0.5 * (m[i][j] + m[j][i]);
        }
    MassMatrix a = m;
    for (int iteration = 0; iteration < 64; ++iteration)
    {
        size_t p = 0, q = 1;
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = i + 1; j < 4; ++j)
                if (std::abs(a[i][j]) > std::abs(a[p][q])) { p = i; q = j; }
        if (std::abs(a[p][q]) <= 1e-14 * scale) break;
        double angle = 0.5 * std::atan2(2 * a[p][q], a[q][q] - a[p][p]);
        double c = std::cos(angle), s = std::sin(angle);
        double pp = c*c*a[p][p] - 2*s*c*a[p][q] + s*s*a[q][q];
        double qq = s*s*a[p][p] + 2*s*c*a[p][q] + c*c*a[q][q];
        for (size_t k = 0; k < 4; ++k)
            if (k != p && k != q)
            {
                double kp = c*a[k][p] - s*a[k][q], kq = s*a[k][p] + c*a[k][q];
                a[k][p] = a[p][k] = kp; a[k][q] = a[q][k] = kq;
            }
        a[p][p] = pp; a[q][q] = qq; a[p][q] = a[q][p] = 0;
    }
    for (size_t i = 0; i < 4; ++i) if (a[i][i] < -1e-12 * scale) return false;
    return true;
}

inline constexpr size_t cPlaneI[6] = { 0, 0, 0, 1, 1, 2 };
inline constexpr size_t cPlaneJ[6] = { 1, 2, 3, 2, 3, 3 };
using InertiaMatrix = std::array<MassBivector, 6>;

inline InertiaMatrix Inertia(const MassMatrix &m)
{
    InertiaMatrix result {};
    for (size_t k = 0; k < 6; ++k)
    {
        MassMatrix omega {};
        omega[cPlaneI[k]][cPlaneJ[k]] = 1;
        omega[cPlaneJ[k]][cPlaneI[k]] = -1;
        for (size_t n = 0; n < 6; ++n)
            for (size_t j = 0; j < 4; ++j)
                result[n][k] += omega[cPlaneI[n]][j] * m[j][cPlaneJ[n]] + m[cPlaneI[n]][j] * omega[j][cPlaneJ[n]];
    }
    return result;
}

inline bool Factor(const MassMatrix &m, InertiaMatrix &l)
{
    InertiaMatrix a = Inertia(m);
    double scale = 0;
    for (size_t i = 0; i < 6; ++i) scale = std::max(scale, a[i][i]);
    if (!(scale > 0) || !std::isfinite(scale)) return false;
    l = {};
    for (size_t i = 0; i < 6; ++i)
        for (size_t j = 0; j <= i; ++j)
        {
            double x = a[i][j];
            for (size_t k = 0; k < j; ++k) x -= l[i][k] * l[j][k];
            if (i == j)
            {
                if (!std::isfinite(x) || x <= 1e-10 * scale) return false;
                l[i][j] = std::sqrt(x);
            }
            else l[i][j] = x / l[j][j];
        }
    return true;
}

inline MassBivector Wedge(const MassVector &a, const MassVector &b)
{
    MassBivector r;
    for (size_t k = 0; k < 6; ++k) r[k] = a[cPlaneI[k]] * b[cPlaneJ[k]] - a[cPlaneJ[k]] * b[cPlaneI[k]];
    return r;
}

inline MassVector Contract(const MassBivector &omega, const MassVector &r)
{
    MassVector out {};
    for (size_t k = 0; k < 6; ++k)
    {
        out[cPlaneI[k]] += omega[k] * r[cPlaneJ[k]];
        out[cPlaneJ[k]] -= omega[k] * r[cPlaneI[k]];
    }
    return out;
}

inline bool AngularBounded(const MassBivector &v)
{
    for (double x : v) if (!std::isfinite(x) || std::abs(x) > 1e9) return false;
    return true;
}

inline MassVector Sub(const MassVector &a, const MassVector &b)
{
    MassVector r;
    for (size_t i = 0; i < 4; ++i) r[i] = a[i] - b[i];
    return r;
}

inline MassVector Transform(const MassMatrix &a, const MassVector &v)
{
    MassVector r {};
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 4; ++j) r[i] += a[i][j] * v[j];
    return r;
}

inline MassMatrix Rotate(const MassMatrix &r, const MassMatrix &m)
{
    MassMatrix out {};
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 4; ++j)
            for (size_t k = 0; k < 4; ++k)
                for (size_t l = 0; l < 4; ++l) out[i][j] += r[i][k] * m[k][l] * r[j][l];
    return out;
}

/// Normalize the two SU(2) factors in double, after validating the supplied Spin(4) rotor.
/// Avoids a float rotation matrix's spurious scale in mass and topology calculations.
inline bool Rotation(const Rotor &rotor, MassMatrix &out)
{
    double v[8];
    for (size_t i = 0; i < 8; ++i)
    {
        v[i] = double(rotor.GetVec8()[uint(i)]);
        if (!std::isfinite(v[i])) return false;
    }
    MassVector l { v[0] + v[7], v[1] - v[6], v[2] + v[5], v[3] - v[4] };
    MassVector r { v[0] - v[7], v[1] + v[6], v[2] - v[5], v[3] + v[4] };
    double ll = 0, rr = 0;
    for (size_t i = 0; i < 4; ++i) { ll += l[i] * l[i]; rr += r[i] * r[i]; }
    if (std::abs(ll - 1) > 2e-5 || std::abs(rr - 1) > 2e-5) return false;
    for (size_t i = 0; i < 4; ++i) { l[i] /= std::sqrt(ll); r[i] /= std::sqrt(rr); }
    double s = 0.5 * (l[0] + r[0]), a = 0.5 * (l[1] + r[1]), b = 0.5 * (l[2] + r[2]);
    double c = 0.5 * (l[3] + r[3]), d = 0.5 * (r[3] - l[3]), e = 0.5 * (l[2] - r[2]);
    double f = 0.5 * (r[1] - l[1]), p = 0.5 * (l[0] - r[0]);
    // The four columns are the sandwich R e_i reverse(R).
    for (size_t i = 0; i < 4; ++i)
    {
        double x = i == 0? 1 : 0, y = i == 1? 1 : 0, z = i == 2? 1 : 0, w = i == 3? 1 : 0;
        double t1 = s*x + a*y + b*z + c*w, t2 = s*y - a*x + d*z + e*w;
        double t3 = s*z - b*x - d*y + f*w, t4 = s*w - c*x - e*y - f*z;
        double t123 = d*x - b*y + a*z + p*w, t124 = e*x - c*y + a*w - p*z;
        double t134 = f*x + p*y - c*z + b*w, t234 = -p*x + f*y - e*z + d*w;
        out[0][i] = t1*s + t2*a + t3*b + t4*c + t123*d + t124*e + t134*f + t234*p;
        out[1][i] = t2*s - t1*a + t3*d + t4*e - t123*b - t124*c + t234*f - t134*p;
        out[2][i] = t3*s - t1*b - t2*d + t4*f + t123*a - t134*c - t234*e + t124*p;
        out[3][i] = t4*s - t1*c - t2*e - t3*f + t124*a + t134*b + t234*d - t123*p;
    }
    return true;
}
}

/// Pure, instantaneous mass/topology preparation in one nonrotating frame and epoch.
/// Units are caller-consistent (typically kg, m, s). IDs are nonzero stable uint64 values;
/// component IDs must be unique across a docking set. Exhaust IDs use a separate namespace.
///
/// Numerical envelope / rejection policy:
/// - 1..4096 components per body/result, <=4096 docking bodies/exhaust lumps; no unbounded
///   topology expansion. Positive input/total mass in [1e-12, 1e24], |tensor entry| <=1e36.
/// - Finite coordinate, linear-velocity and angular-velocity components bounded by 1e9;
///   world COM and queried point velocities must also be inside that bound.
/// - Spin(4) factor norm-squared error <=2e-5, projected in double for calculations.
///   Tensor asymmetry <=1e-6 of max entry is symmetrized. Eigenvalues below -1e-12 of
///   that scale reject; smaller PSD roundoff is tolerated. Aggregate bivector inertia
///   must have every Cholesky pivot >1e-10 of its largest diagonal. No pseudo-inverse.
/// - Mass closure 1e-12 relative; momentum and exhaust first moment closure 1e-9 of
///   sums of absolute participant scales (floor 1e-12). Residuals are published, not hidden.
///   Free fragmentation accepts |delta K| <=1e-9*max(K_before, 1e-12) + B_v. B_v is
///   sum_{fragment,axis} m*(|u|*|e| + e^2/2), where u = omega*(COM_fragment-COM_before)
///   and e is the actual rounding residue of storing V_before+u as an absolute double.
///   TwoSum measures that residue; mVelocityRoundingEnergyBound publishes the bound.
///   Actual energy differences remain in the work/dissipation entries, never clamped away.
///   Docking may add <=1e-6*max(K_before, 1e-12) as a float-Rotor composition budget.
/// - Ledger kinetic energies use the incoming system's COM-velocity frame to resolve
///   small work at large common speed; momenta remain world components. Derive() alone
///   reports world-frame kinetic energy. Central spin solves use relative velocities and
///   nearby COM offsets, never a difference of world orbital angular momenta. Absolute
///   global angular momentum is inherently ill-conditioned far from its origin; ledger
///   origins are nearby world COMs. Already-quantized input velocities cannot be recovered.
///   Rounding analysis assumes binary64 round-to-nearest without fast-math reassociation;
///   subnormal-scale accuracy and arbitrary ill-conditioned cancellations are not promised.
///
/// Evidence envelope (MassModelTests): independent 32-particle, general Spin(4) references
/// and 100 seeded redistribution/fragment/remerge trials at component masses 2..5,
/// offsets <=2, speeds <=3, spin components <=2; 100 multi-lump exhaust trials; exact
/// analytical examples; 1e6:1 mass-ratio docking; 1e8 coordinate recentering/quarter-unit
/// aggregation; unequal-mass 0.1-offset docking / 0.3-offset fragmentation at speed 1e8;
/// 100 additional general-Spin(4) fragment/remerge trials with double-distributed offsets
/// <=2, origins <=10 and independent common boosts <=1e8 in all four lanes; 100 boosted
/// multi-point exhaust trials (masses 2..5, offsets <=1, relative speeds <=3). These check
/// absolute-velocity quantization against an independent long-double bound and central
/// spin covariance for identical represented relative velocities, including exact zero
/// spin/work for co-moving participants. Limits above prevent overflow, NOT a uniform
/// accuracy guarantee over every mass ratio/conditioning. Rebase tiny distant features.
///
/// Not implemented here: live physics commit/locking/cache/shape replacement, collision
/// separation, joints, chemical/thermal resource authorization, separation impulses,
/// finite-duration fuel flow or gyroscopic integration. Input central tensors must match
/// actual component geometry; a PSD check cannot certify that. All allocations complete
/// before output replacement; process-level allocation failure follows the host's policy
/// (Jolt normally disables exceptions). No cross-platform bitwise determinism claim.
class MassModel
{
public:
    /// referenceOrigin is common to every participant in a ledger. Default is global zero.
    static bool Derive(const MassBody &body, BodyMassData &out, const MassVector &referenceOrigin = {})
    {
        BodyMassData result;
        MassMatrix rotation;
        if (body.mID == 0 || !MassDetail::Bounded(body.mOrigin) || !MassDetail::Bounded(referenceOrigin)
            || !MassDetail::Bounded(body.mLinearVelocity) || !MassDetail::AngularBounded(body.mAngularVelocity)
            || !MassDetail::Rotation(body.mRotation, rotation) || !Aggregate(body.mComponents, result.mLocal)) return false;
        result.mWorldCOM = MassDetail::Transform(rotation, result.mLocal.mLocalCOM);
        for (size_t i = 0; i < 4; ++i) result.mWorldCOM[i] += body.mOrigin[i];
        if (!MassDetail::Bounded(result.mWorldCOM)) return false;
        result.mWorldSecondMoment = MassDetail::Rotate(rotation, result.mLocal.mSecondMoment);
        auto inertia = MassDetail::Inertia(result.mWorldSecondMoment);
        for (size_t i = 0; i < 6; ++i)
            for (size_t j = 0; j < 6; ++j) result.mSpinMomentum[i] += inertia[i][j] * body.mAngularVelocity[j];
        result.mMomentum.mMass = result.mLocal.mMass;
        for (size_t i = 0; i < 4; ++i)
        {
            result.mMomentum.mLinear[i] = result.mLocal.mMass * body.mLinearVelocity[i];
            result.mMomentum.mKineticEnergy += 0.5 * result.mMomentum.mLinear[i] * body.mLinearVelocity[i];
        }
        MassBivector orbital = MassDetail::Wedge(result.mMomentum.mLinear, MassDetail::Sub(result.mWorldCOM, referenceOrigin));
        for (size_t i = 0; i < 6; ++i)
        {
            result.mMomentum.mAngular[i] = result.mSpinMomentum[i] + orbital[i];
            result.mMomentum.mKineticEnergy += 0.5 * result.mSpinMomentum[i] * body.mAngularVelocity[i];
        }
        if (!std::isfinite(result.mMomentum.mKineticEnergy) || result.mMomentum.mKineticEnergy < 0) return false;
        out = result;
        return true;
    }

    static bool PointVelocity(const MassBody &body, const MassVector &worldPoint, MassVector &out)
    {
        BodyMassData data;
        if (!MassDetail::Bounded(worldPoint) || !Derive(body, data)) return false;
        MassVector velocity = MassDetail::Contract(body.mAngularVelocity, MassDetail::Sub(worldPoint, data.mWorldCOM));
        for (size_t i = 0; i < 4; ++i) velocity[i] += body.mLinearVelocity[i];
        if (!MassDetail::Bounded(velocity)) return false;
        out = velocity;
        return true;
    }

    /// Pure reference-point change: shifts component local positions opposite to the origin.
    /// COM velocity is NOT geometry-origin velocity and must not receive an angular correction.
    static bool Recenter(const MassBody &body, MassBody &out)
    {
        BodyMassData data;
        if (!Derive(body, data)) return false;
        MassBody result = body;
        result.mOrigin = data.mWorldCOM;
        for (auto &c : result.mComponents) c.mLocalCOM = MassDetail::Sub(c.mLocalCOM, data.mLocal.mLocalCOM);
        if (!Derive(result, data)) return false;
        out = std::move(result);
        return true;
    }
    /// Isolated internal rearrangement: keep world COM fixed (the geometry origin recoils).
    /// Component identity set and total mass must not change; per-component fuel transfers may.
    /// Energy difference is explicit required work / dissipation, not a conservation claim.
    static bool Redistribute(const MassBody &body, const std::vector<ComponentMass> &components, MassTransaction &out)
    {
        BodyMassData before;
        AggregateMass target;
        if (!Derive(body, before) || !Aggregate(components, target)
            || !SameIDs(body.mComponents, components) || !NearMass(before.mLocal.mMass, target.mMass)) return false;
        MassTransaction result;
        MassBody changed = body;
        changed.mComponents = components;
        MassMatrix rotation;
        if (!MassDetail::Rotation(body.mRotation, rotation)) return false;
        MassVector shift = MassDetail::Transform(rotation, MassDetail::Sub(before.mLocal.mLocalCOM, target.mLocalCOM));
        for (size_t i = 0; i < 4; ++i) changed.mOrigin[i] += shift[i];
        if (!SetCentralMomentum(changed, body.mLinearVelocity, before.mSpinMomentum)) return false;
        result.mBodies.push_back(std::move(changed));
        if (!Ledger({ body }, before.mWorldCOM, {}, result)) return false;
        out = std::move(result);
        return true;
    }

    /// Force-free, impulse-free partition only. Geometry frame/orientation are inherited exactly.
    /// Every component must appear exactly once; each fragment must have invertible inertia.
    static bool Fragment(const MassBody &body, const std::vector<FragmentRequest> &requests, MassTransaction &out)
    {
        BodyMassData initial;
        if (!Derive(body, initial) || requests.size() < 2 || requests.size() > body.mComponents.size()) return false;
        size_t count = 0;
        for (const auto &request : requests)
        {
            if (request.mComponentIDs.size() > body.mComponents.size() - count) return false;
            count += request.mComponentIDs.size();
        }
        if (count != body.mComponents.size()) return false;
        MassTransaction result;
        std::vector<uint64> used, bodyIDs;
        MassDetail::Sum velocityRoundingEnergyBound;
        used.reserve(count); bodyIDs.reserve(requests.size());
        for (const auto &request : requests)
        {
            if (request.mBodyID == 0 || request.mComponentIDs.empty()) return false;
            bodyIDs.push_back(request.mBodyID);
            MassBody piece = body;
            piece.mID = request.mBodyID;
            piece.mComponents.clear();
            for (uint64 id : request.mComponentIDs)
            {
                auto found = std::find_if(body.mComponents.begin(), body.mComponents.end(), [id](const ComponentMass &c) { return c.mID == id; });
                if (found == body.mComponents.end()) return false;
                used.push_back(id);
                piece.mComponents.push_back(*found);
            }
            BodyMassData data;
            if (!Derive(piece, data)) return false;
            MassVector velocity = MassDetail::Contract(body.mAngularVelocity, MassDetail::Sub(data.mWorldCOM, initial.mWorldCOM));
            for (size_t i = 0; i < 4; ++i)
            {
                double a = body.mLinearVelocity[i], b = velocity[i], sum = a + b;
                piece.mLinearVelocity[i] = sum;
                // TwoSum recovers the actual rounding residue of the absolute-velocity
                // addition (no overflow in the validated envelope). For e = stored - ideal,
                // |delta K| <= m (|v_relative| |e| + e^2 / 2), in the incoming COM frame.
                double bv = sum - a;
                double e = std::abs((a - (sum - bv)) + (b - bv));
                velocityRoundingEnergyBound.Add(data.mLocal.mMass * (std::abs(b) * e + 0.5 * e * e));
            }
            result.mBodies.push_back(std::move(piece));
        }
        std::sort(used.begin(), used.end()); std::sort(bodyIDs.begin(), bodyIDs.end());
        if (used.size() != body.mComponents.size() || std::adjacent_find(used.begin(), used.end()) != used.end()
            || std::adjacent_find(bodyIDs.begin(), bodyIDs.end()) != bodyIDs.end()) return false;
        std::sort(result.mBodies.begin(), result.mBodies.end(), [](const MassBody &a, const MassBody &b) { return a.mID < b.mID; });
        if (!Ledger({ body }, initial.mWorldCOM, {}, result)) return false;
        result.mLedger.mVelocityRoundingEnergyBound = velocityRoundingEnergyBound.Get();
        double energyError = result.mLedger.mAddedEnergy + result.mLedger.mDissipatedEnergy;
        if (energyError > 1e-9 * std::max(result.mLedger.mBefore.mKineticEnergy, 1e-12)
            + result.mLedger.mVelocityRoundingEnergyBound) return false;
        out = std::move(result);
        return true;
    }

    /// Inelastic rigid merge, NOT a docking constraint. New frame is world-aligned at the
    /// lowest-ID body's COM. Component rotations are composed once in float Rotor storage.
    static bool Dock(const std::vector<MassBody> &bodies, uint64 newBodyID, MassTransaction &out)
    {
        if (bodies.size() < 2 || bodies.size() > 4096 || newBodyID == 0) return false;
        size_t count = 0;
        for (const auto &body : bodies)
        {
            if (body.mComponents.empty() || body.mComponents.size() > 4096 - count) return false;
            count += body.mComponents.size();
        }
        std::vector<MassBody> ordered = bodies;
        std::sort(ordered.begin(), ordered.end(), [](const MassBody &a, const MassBody &b) { return a.mID < b.mID; });
        BodyMassData first;
        if (!Derive(ordered.front(), first)) return false;
        MassBody merged;
        merged.mID = newBodyID;
        merged.mOrigin = first.mWorldCOM;
        MomentumSum sum;
        std::vector<BodyMassData> incoming;
        uint64 previous = 0;
        for (const auto &body : ordered)
        {
            BodyMassData data;
            if (body.mID == previous || !Derive(body, data, first.mWorldCOM)) return false;
            previous = body.mID;
            incoming.push_back(data);
            MassMomentum relative = data.mMomentum;
            for (size_t i = 0; i < 4; ++i)
                relative.mLinear[i] = data.mLocal.mMass * (body.mLinearVelocity[i] - ordered.front().mLinearVelocity[i]);
            sum.Add(relative);
            MassMatrix rotation;
            if (!MassDetail::Rotation(body.mRotation, rotation)) return false;
            for (const auto &component : body.mComponents)
            {
                if (merged.mComponents.size() >= 4096) return false;
                ComponentMass c = component;
                MassVector offset = MassDetail::Transform(rotation, c.mLocalCOM);
                for (size_t i = 0; i < 4; ++i) c.mLocalCOM[i] = (body.mOrigin[i] - merged.mOrigin[i]) + offset[i];
                c.mLocalRotation = (body.mRotation.Normalized() * component.mLocalRotation.Normalized()).Normalized();
                merged.mComponents.push_back(c);
            }
        }
        std::sort(merged.mComponents.begin(), merged.mComponents.end(), [](const ComponentMass &a, const ComponentMass &b) { return a.mID < b.mID; });
        MassMomentum target = sum.Get();
        MassVector mean {}, velocity {};
        for (size_t i = 0; i < 4; ++i)
        {
            mean[i] = target.mLinear[i] / target.mMass;
            velocity[i] = ordered.front().mLinearVelocity[i] + mean[i];
        }
        BodyMassData combined;
        if (!Derive(merged, combined)) return false;
        MassDetail::Sum central[6];
        for (size_t n = 0; n < ordered.size(); ++n)
        {
            MassVector p;
            for (size_t i = 0; i < 4; ++i)
                p[i] = incoming[n].mLocal.mMass * ((ordered[n].mLinearVelocity[i] - ordered.front().mLinearVelocity[i]) - mean[i]);
            auto orbital = MassDetail::Wedge(p, MassDetail::Sub(incoming[n].mWorldCOM, combined.mWorldCOM));
            for (size_t i = 0; i < 6; ++i) central[i].Add(incoming[n].mSpinMomentum[i] + orbital[i]);
        }
        MassBivector spin;
        for (size_t i = 0; i < 6; ++i) spin[i] = central[i].Get();
        if (!SetCentralMomentum(merged, velocity, spin)) return false;
        MassTransaction result;
        result.mBodies.push_back(std::move(merged));
        if (!Ledger(ordered, first.mWorldCOM, {}, result)) return false;
        // Float Rotor recomposition can perturb inertia at its storage precision.
        if (result.mLedger.mAddedEnergy > 1e-6 * std::max(result.mLedger.mBefore.mKineticEnergy, 1e-12)) return false;
        out = std::move(result);
        return true;
    }

    /// Instantaneous lumped flux: subtract exported relative p and p wedge x in the incoming
    /// COM-velocity frame, then solve central spin; world momenta are a separate ledger.
    /// Do NOT also apply thrust/recoil for the same exhaust. Retained component masses
    /// cannot increase; removed mass and first spatial moment must close exactly to
    /// the declared tolerance. Central shape changes are caller's explicit lump approximation.
    static bool Expel(const MassBody &body, const std::vector<ComponentMass> &retained,
        const std::vector<LumpedExhaust> &exhaust, MassTransaction &out)
    {
        BodyMassData initial, rest;
        if (!Derive(body, initial) || exhaust.empty() || exhaust.size() > 4096 || retained.empty() || retained.size() > 4096) return false;
        MassBody remaining = body;
        remaining.mComponents = retained;
        if (!Derive(remaining, rest) || !(rest.mLocal.mMass < initial.mLocal.mMass)) return false;
        for (const auto &c : retained)
        {
            auto found = std::find_if(body.mComponents.begin(), body.mComponents.end(), [&](const ComponentMass &old) { return old.mID == c.mID; });
            if (found == body.mComponents.end() || c.mProperties.mMass > found->mProperties.mMass) return false;
        }
        std::vector<LumpedExhaust> ordered = exhaust;
        std::sort(ordered.begin(), ordered.end(), [](const LumpedExhaust &a, const LumpedExhaust &b) { return a.mID < b.mID; });
        MomentumSum sum;
        MassDetail::Sum first[4];
        MassDetail::Sum relativeLinear[4], relativeAngular[6];
        double firstScale = 0;
        MassVector restOffset = MassDetail::Sub(rest.mWorldCOM, initial.mWorldCOM);
        for (size_t i = 0; i < 4; ++i)
        {
            first[i].Add(rest.mLocal.mMass * restOffset[i]);
            firstScale += std::abs(rest.mLocal.mMass * restOffset[i]);
        }
        uint64 previous = 0;
        for (const auto &gas : ordered)
        {
            if (gas.mID == 0 || gas.mID == previous || !std::isfinite(gas.mMass) || gas.mMass < 1e-12 || gas.mMass > 1e24
                || !MassDetail::Bounded(gas.mWorldPosition) || !MassDetail::Bounded(gas.mWorldVelocity)) return false;
            previous = gas.mID;
            MassMomentum m;
            m.mMass = gas.mMass;
            MassVector offset = MassDetail::Sub(gas.mWorldPosition, initial.mWorldCOM);
            for (size_t i = 0; i < 4; ++i)
            {
                m.mLinear[i] = gas.mMass * gas.mWorldVelocity[i];
                m.mKineticEnergy += 0.5 * m.mLinear[i] * gas.mWorldVelocity[i];
                first[i].Add(gas.mMass * offset[i]);
                firstScale += std::abs(gas.mMass * offset[i]);
            }
            m.mAngular = MassDetail::Wedge(m.mLinear, offset);
            sum.Add(m);
            MassVector p;
            for (size_t i = 0; i < 4; ++i)
            {
                p[i] = gas.mMass * (gas.mWorldVelocity[i] - body.mLinearVelocity[i]);
                relativeLinear[i].Add(p[i]);
            }
            auto orbital = MassDetail::Wedge(p, offset);
            for (size_t i = 0; i < 6; ++i) relativeAngular[i].Add(orbital[i]);
        }
        MassMomentum exported = sum.Get();
        if (!NearMass(rest.mLocal.mMass + exported.mMass, initial.mLocal.mMass)) return false;
        for (const auto &moment : first)
            if (std::abs(moment.Get()) > 1e-9 * std::max(firstScale, 1e-12)) return false;
        MassVector linear, velocity;
        MassBivector angular;
        for (size_t i = 0; i < 4; ++i)
        {
            linear[i] = -relativeLinear[i].Get();
            velocity[i] = body.mLinearVelocity[i] + linear[i] / rest.mLocal.mMass;
        }
        auto orbital = MassDetail::Wedge(linear, restOffset);
        for (size_t i = 0; i < 6; ++i)
            angular[i] = initial.mSpinMomentum[i] - relativeAngular[i].Get() - orbital[i];
        if (!SetCentralMomentum(remaining, velocity, angular)) return false;
        MassTransaction result;
        result.mBodies.push_back(std::move(remaining));
        result.mExhaust = std::move(ordered);
        if (!Ledger({ body }, initial.mWorldCOM, exported, result)) return false;
        out = std::move(result);
        return true;
    }

    /// Checked narrowing for the physics adapter. Does not diagonalize or encode an inertia
    /// tensor as a Rotor. Adapter must budget float quantization when installing velocities.
    static bool ToMassProperties(const AggregateMass &aggregate, MassProperties &out)
    {
        if (!std::isfinite(aggregate.mMass) || aggregate.mMass < 1e-12 || aggregate.mMass > 1e24
            || !MassDetail::Bounded(aggregate.mLocalCOM)) return false;
        MassMatrix m = aggregate.mSecondMoment;
        MassDetail::InertiaMatrix factor;
        if (!MassDetail::Tensor(m) || !MassDetail::Factor(m, factor)) return false;
        MassProperties result;
        result.mMass = float(aggregate.mMass);
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = 0; j < 4; ++j)
            {
                result.mInertia(uint(i), uint(j)) = float(m[i][j]);
                m[i][j] = double(result.mInertia(uint(i), uint(j)));
            }
        if (!MassDetail::Tensor(m) || !MassDetail::Factor(m, factor)) return false;
        out = result;
        return true;
    }

    /// Two-pass centered accumulation avoids subtracting two enormous raw second moments.
    /// All prepare APIs leave their output untouched on rejection.
    static bool Aggregate(const std::vector<ComponentMass> &components, AggregateMass &out)
    {
        if (components.empty() || components.size() > 4096) return false;
        std::vector<const ComponentMass *> ordered;
        ordered.reserve(components.size());
        for (const auto &c : components) ordered.push_back(&c);
        std::sort(ordered.begin(), ordered.end(), [](const ComponentMass *a, const ComponentMass *b) { return a->mID < b->mID; });
        uint64 previous = 0;
        for (const ComponentMass *cp : ordered)
        {
            const ComponentMass &c = *cp;
            if (c.mID == 0 || c.mID == previous || !MassDetail::Bounded(c.mLocalCOM)) return false;
            previous = c.mID;
        }
        AggregateMass result;
        MassDetail::Sum mass, first[4];
        const MassVector anchor = ordered.front()->mLocalCOM;
        for (const ComponentMass *cp : ordered)
        {
            const ComponentMass &c = *cp;
            if (!std::isfinite(c.mProperties.mMass) || c.mProperties.mMass < 1e-12 || c.mProperties.mMass > 1e24) return false;
            mass.Add(double(c.mProperties.mMass));
            for (size_t i = 0; i < 4; ++i)
                first[i].Add(double(c.mProperties.mMass) * (c.mLocalCOM[i] - anchor[i]));
        }
        result.mMass = mass.Get();
        if (!std::isfinite(result.mMass) || result.mMass > 1e24) return false;
        for (size_t i = 0; i < 4; ++i) result.mLocalCOM[i] = anchor[i] + first[i].Get() / result.mMass;
        MassDetail::Sum second[4][4];
        for (const ComponentMass *cp : ordered)
        {
            const ComponentMass &c = *cp;
            MassMatrix r, m;
            if (!MassDetail::Rotation(c.mLocalRotation, r)) return false;
            for (size_t i = 0; i < 4; ++i)
                for (size_t j = 0; j < 4; ++j) m[i][j] = double(c.mProperties.mInertia(uint(i), uint(j)));
            if (!MassDetail::Tensor(m)) return false;
            m = MassDetail::Rotate(r, m);
            MassVector delta = MassDetail::Sub(c.mLocalCOM, result.mLocalCOM);
            for (size_t i = 0; i < 4; ++i)
                for (size_t j = 0; j < 4; ++j)
                    second[i][j].Add(m[i][j] + double(c.mProperties.mMass) * delta[i] * delta[j]);
        }
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = 0; j < 4; ++j) result.mSecondMoment[i][j] = second[i][j].Get();
        MassDetail::InertiaMatrix factor;
        if (!MassDetail::Tensor(result.mSecondMoment) || !MassDetail::Factor(result.mSecondMoment, factor)) return false;
        out = result;
        return true;
    }

private:
    static bool NearMass(double a, double b) { return std::abs(a - b) <= 1e-12 * std::max(a, b); }

    static bool SameIDs(const std::vector<ComponentMass> &a, const std::vector<ComponentMass> &b)
    {
        if (a.size() != b.size()) return false;
        std::vector<uint64> x, y;
        for (const auto &c : a) x.push_back(c.mID);
        for (const auto &c : b) y.push_back(c.mID);
        std::sort(x.begin(), x.end()); std::sort(y.begin(), y.end());
        return x == y;
    }

    /// Solve central spin directly; never recover it by subtracting world orbital momenta.
    static bool SetCentralMomentum(MassBody &body, const MassVector &velocity, const MassBivector &spin)
    {
        BodyMassData data;
        if (!Derive(body, data)) return false;
        MassDetail::InertiaMatrix l;
        if (!MassDetail::Factor(data.mWorldSecondMoment, l)) return false;
        MassBivector y {}, omega {};
        for (size_t i = 0; i < 6; ++i)
        {
            double v = spin[i];
            for (size_t j = 0; j < i; ++j) v -= l[i][j] * y[j];
            y[i] = v / l[i][i];
        }
        for (size_t n = 6; n-- > 0;)
        {
            double v = y[n];
            for (size_t j = n + 1; j < 6; ++j) v -= l[j][n] * omega[j];
            omega[n] = v / l[n][n];
        }
        body.mAngularVelocity = omega;
        body.mLinearVelocity = velocity;
        return Derive(body, data);
    }

    struct MomentumSum
    {
        MassDetail::Sum mMass, mLinear[4], mAngular[6], mEnergy;
        void Add(const MassMomentum &m)
        {
            mMass.Add(m.mMass); mEnergy.Add(m.mKineticEnergy);
            for (size_t i = 0; i < 4; ++i) mLinear[i].Add(m.mLinear[i]);
            for (size_t i = 0; i < 6; ++i) mAngular[i].Add(m.mAngular[i]);
        }
        MassMomentum Get() const
        {
            MassMomentum m;
            m.mMass = mMass.Get(); m.mKineticEnergy = mEnergy.Get();
            for (size_t i = 0; i < 4; ++i) m.mLinear[i] = mLinear[i].Get();
            for (size_t i = 0; i < 6; ++i) m.mAngular[i] = mAngular[i].Get();
            return m;
        }
    };

    static bool Ledger(const std::vector<MassBody> &before, const MassVector &reference, const MassMomentum &exported, MassTransaction &result)
    {
        MassLedger ledger;
        ledger.mReferenceOrigin = reference;
        ledger.mExported = exported;
        MomentumSum initial;
        for (const auto &body : before)
        {
            BodyMassData data;
            if (!Derive(body, data, reference)) return false;
            for (size_t i = 0; i < 4; ++i)
                data.mMomentum.mLinear[i] = data.mLocal.mMass * (body.mLinearVelocity[i] - before.front().mLinearVelocity[i]);
            initial.Add(data.mMomentum);
        }
        MassMomentum total = initial.Get();
        if (!(total.mMass > 0)) return false;
        for (size_t i = 0; i < 4; ++i)
            ledger.mEnergyFrameVelocity[i] = before.front().mLinearVelocity[i] + total.mLinear[i] / total.mMass;
        // Momentum closure makes kinetic-work differences invariant under this common boost.
        // Evaluate squared relative velocities directly, not by subtracting huge world energies.
        ledger.mExported.mKineticEnergy = 0;
        MassDetail::Sum exhaustEnergy;
        for (const auto &gas : result.mExhaust)
            for (size_t i = 0; i < 4; ++i)
            {
                double v = gas.mWorldVelocity[i] - ledger.mEnergyFrameVelocity[i];
                exhaustEnergy.Add(0.5 * gas.mMass * v * v);
            }
        ledger.mExported.mKineticEnergy = exhaustEnergy.Get();
        MomentumSum a, b;
        double linearScale = 0, angularScale = 0;
        auto collect = [&](const std::vector<MassBody> &bodies, MomentumSum &sum)
        {
            for (const auto &body : bodies)
            {
                BodyMassData data;
                if (!Derive(body, data, reference)) return false;
                MassDetail::Sum energy;
                for (size_t i = 0; i < 4; ++i)
                {
                    double v = body.mLinearVelocity[i] - ledger.mEnergyFrameVelocity[i];
                    energy.Add(0.5 * data.mLocal.mMass * v * v);
                }
                for (size_t i = 0; i < 6; ++i) energy.Add(0.5 * data.mSpinMomentum[i] * body.mAngularVelocity[i]);
                data.mMomentum.mKineticEnergy = energy.Get();
                sum.Add(data.mMomentum);
                for (double x : data.mMomentum.mLinear) linearScale += std::abs(x);
                for (double x : data.mMomentum.mAngular) angularScale += std::abs(x);
                // A zero net momentum must not imply a zero numerical error budget.
                for (double x : data.mSpinMomentum) angularScale += std::abs(x);
            }
            return true;
        };
        if (!collect(before, a) || !collect(result.mBodies, b)) return false;
        ledger.mBefore = a.Get(); ledger.mAfter = b.Get();
        for (double x : exported.mLinear) linearScale += std::abs(x);
        for (double x : exported.mAngular) angularScale += std::abs(x);
        ledger.mResidual.mMass = ledger.mAfter.mMass + exported.mMass - ledger.mBefore.mMass;
        if (!NearMass(ledger.mAfter.mMass + exported.mMass, ledger.mBefore.mMass)) return false;
        for (size_t i = 0; i < 4; ++i)
        {
            double residual = ledger.mAfter.mLinear[i] + exported.mLinear[i] - ledger.mBefore.mLinear[i];
            ledger.mResidual.mLinear[i] = residual;
            if (!std::isfinite(residual) || std::abs(residual) > 1e-9 * std::max(linearScale, 1e-12)) return false;
        }
        for (size_t i = 0; i < 6; ++i)
        {
            double residual = ledger.mAfter.mAngular[i] + exported.mAngular[i] - ledger.mBefore.mAngular[i];
            ledger.mResidual.mAngular[i] = residual;
            if (!std::isfinite(residual) || std::abs(residual) > 1e-9 * std::max(angularScale, 1e-12)) return false;
        }
        double delta = ledger.mAfter.mKineticEnergy + ledger.mExported.mKineticEnergy - ledger.mBefore.mKineticEnergy;
        if (!std::isfinite(delta)) return false;
        ledger.mAddedEnergy = std::max(delta, 0.0);
        ledger.mDissipatedEnergy = std::max(-delta, 0.0);
        ledger.mResidual.mKineticEnergy = delta - ledger.mAddedEnergy + ledger.mDissipatedEnergy;
        result.mLedger = ledger;
        return true;
    }
};
}
JPH_NAMESPACE_END
