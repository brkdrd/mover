// SPDX-License-Identifier: MIT
// Diagnostic acceptance gate, deliberately separate from the green regression suite.
#include <Jolt/Jolt.h>
#include <Jolt/Core/Memory.h>
#include <Jolt/Geometry/ConvexSupport.h>
#include <Jolt/Geometry/GJKClosestPoint.h>
#include <Jolt/Geometry/AABox.h>
#include <Jolt/Geometry/Sphere.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

using namespace JPH;

int main()
{
    RegisterDefaultAllocator();
    std::mt19937 random(0x4d6a6bu);
    const Vec4 axes[] = { Vec4::sZero(), Vec4::sAxisX(), Vec4::sAxisW(), Vec4(-1, 2, -3, 4) };
    unsigned queries = 0, misses = 0, false_hits = 0, distance_failures = 0, oracle_failures = 0;
    unsigned selected = 0, selected_misses = 0;
    std::cout << "sample,scale,margin,tolerance,overlap,axis,reversed,hit,distance_sq\n" << std::setprecision(17);
    for (int sample = 0; sample < 64; ++sample)
    {
        float angles[6];
        for (float &angle : angles) angle = float(int(random() & 1023u) - 512) / 512.0f;
        const Rotor rotor = Rotor::sExp(Bivec(angles[0], angles[1], angles[2], angles[3], angles[4], angles[5]));
        const Mat44 rotation = Mat44::sRotation(rotor);
        Vec4 extent, center;
        double corner_distance_sq = 0;
        for (int axis = 0; axis < 4; ++axis)
        {
            extent[axis] = 0.5f + float(random() & 31u) / 32.0f;
            const float delta = 0.25f + float(random() & 15u) / 64.0f;
            center[axis] = (random() & 1u ? 1.0f : -1.0f) * (extent[axis] + delta);
            corner_distance_sq += double(delta) * delta;
        }
        for (float scale : { 0.001f, 0.003f, 0.01f, 0.03f, 0.1f, 1.0f })
        {
            const AABox local_box(-extent * scale, extent * scale);
            const TransformedConvexObject<AABox> box(rotation, Vec4::sZero(), local_box);
            for (float margin : { -0.125f, -0.001f, -0.00005f, 0.00005f, 0.001f, 0.125f })
            {
                const Sphere sphere(rotor * (center * scale), float(std::sqrt(corner_distance_sq) + double(margin)) * scale);
                double residual[4], distance_sq = 0;
                for (int j = 0; j < 4; ++j)
                {
                    double corner = 0;
                    for (int k = 0; k < 4; ++k)
                        corner += double(rotation.GetColumn4(k)[j]) * double(center[k] > 0 ? local_box.mMax[k] : local_box.mMin[k]);
                    residual[j] = double(sphere.GetCenter()[j]) - corner;
                    distance_sq += residual[j] * residual[j];
                }
                // Verify that the prescribed corner really is the closest point:
                // the residual points outward along every transformed edge.
                for (int k = 0; k < 4; ++k)
                {
                    double outward = 0;
                    for (int j = 0; j < 4; ++j) outward += residual[j] * double(rotation.GetColumn4(k)[j]);
                    if (!(outward * (center[k] > 0 ? 1 : -1) > 0)) ++oracle_failures;
                }
                const double overlap = double(sphere.GetRadius()) - std::sqrt(distance_sq);
                const float tolerance = 1.0e-4f * scale;
                for (unsigned axis = 0; axis < 4; ++axis)
                    for (bool reverse : { false, true })
                    {
                        GJKClosestPoint gjk;
                        Vec4 v = axes[axis], pa, pb;
                        const bool hit = reverse ? gjk.Intersects(box, sphere, tolerance, v) : gjk.Intersects(sphere, box, tolerance, v);
                        v = axes[axis];
                        const float closest_sq = reverse
                            ? gjk.GetClosestPoints(box, sphere, tolerance, cLargeFloat, v, pa, pb)
                            : gjk.GetClosestPoints(sphere, box, tolerance, cLargeFloat, v, pa, pb);
                        ++queries;
                        if (overlap > 2.0 * tolerance && !hit) ++misses;
                        if (overlap < -2.0 * tolerance && hit) ++false_hits;
                        const double expected_distance = std::max(0.0, -overlap);
                        if (!std::isfinite(closest_sq) || closest_sq < 0
                            || std::abs(std::sqrt(double(closest_sq)) - expected_distance) > 2.0 * tolerance) ++distance_failures;
                        if (scale == 0.01f && margin == 0.001f)
                        {
                            ++selected;
                            if (overlap > 2.0 * tolerance && !hit) ++selected_misses;
                        }
                        std::cout << sample << ',' << scale << ',' << margin << ',' << tolerance << ',' << overlap << ','
                            << axis << ',' << int(reverse) << ',' << int(hit) << ',' << closest_sq << '\n';
                    }
            }
        }
    }
    std::cout.flush();
    if (!std::cout) return 2;
    std::cerr << "queries=" << queries << " material_misses=" << misses << " material_false_hits=" << false_hits
        << " distance_failures=" << distance_failures << " oracle_failures=" << oracle_failures
        << " selected_scale_0.01_margin_0.001=" << selected_misses << '/' << selected << '\n';
    if (queries != 18432 || selected != 512 || oracle_failures != 0) return 2;
    return misses || false_hits || distance_failures ? 1 : 0;
}
