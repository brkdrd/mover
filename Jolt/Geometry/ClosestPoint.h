// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Math/Bivec.h>

JPH_NAMESPACE_BEGIN

// Turn off fused multiply add instruction because it makes the equations of the form a * b - c * d inaccurate below
JPH_PRECISE_MATH_ON

/// Helper utils to find the closest point to a line segment, triangle or tetrahedron
namespace ClosestPoint
{
	/// The 2x2 Gram determinant |u|^2 |v|^2 - (u . v)^2 of two vectors, i.e. the squared area of the
	/// parallelogram they span.
	///
	/// TODO(4D): compute it from the wedge product rather than directly. 3D Jolt gets this quantity as
	/// |u x v|^2, and the cross product forms each component as a difference of products of individual
	/// coordinates, which stays small. Evaluating |u|^2 |v|^2 - (u . v)^2 literally instead subtracts
	/// two huge nearly equal numbers: for the sliver triangle in TestLongTriangle it is
	/// 40000 * 40001 - 40000^2, where both products round to a multiple of 128 in float and the
	/// difference lands ~0.16% off. The wedge product is the dimension-agnostic generalization of the
	/// cross product and has exactly the same cancellation-free structure.
	JPH_INLINE float GramDeterminant(Vec4Arg inU, Vec4Arg inV)
	{
		return Bivec::sWedge(inU, inV).LengthSq();
	}

	/// Compute barycentric coordinates of closest point to origin for infinite line defined by (inA, inB)
	/// Point can then be computed as inA * outU + inB * outV
	/// Returns false if the points inA, inB do not form a line (are at the same point)
	inline bool GetBaryCentricCoordinates(Vec4Arg inA, Vec4Arg inB, float &outU, float &outV)
	{
		Vec4 ab = inB - inA;
		float denominator = ab.LengthSq();
		if (denominator < Square(FLT_EPSILON))
		{
			// Degenerate line segment, fallback to points
			if (inA.LengthSq() < inB.LengthSq())
			{
				// A closest
				outU = 1.0f;
				outV = 0.0f;
			}
			else
			{
				// B closest
				outU = 0.0f;
				outV = 1.0f;
			}
			return false;
		}
		else
		{
			outV = -inA.Dot(ab) / denominator;
			outU = 1.0f - outV;
		}
		return true;
	}

	/// Compute barycentric coordinates of closest point to origin for plane defined by (inA, inB, inC)
	/// Point can then be computed as inA * outU + inB * outV + inC * outW
	/// Returns false if the points inA, inB, inC do not form a plane (are on the same line or at the same point)
	inline bool GetBaryCentricCoordinates(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, float &outU, float &outV, float &outW)
	{
		// Taken from: Real-Time Collision Detection - Christer Ericson (Section: Barycentric Coordinates)
		// With p = 0
		// Adjusted to always include the shortest edge of the triangle in the calculation to improve numerical accuracy

		// First calculate the three edges
		Vec4 v0 = inB - inA;
		Vec4 v1 = inC - inA;
		Vec4 v2 = inC - inB;

		// Make sure that the shortest edge is included in the calculation to keep the products a * b - c * d as small as possible to preserve accuracy
		float d00 = v0.LengthSq();
		float d11 = v1.LengthSq();
		float d22 = v2.LengthSq();
		if (d00 <= d22)
		{
			// Use v0 and v1 to calculate barycentric coordinates
			float d01 = v0.Dot(v1);

			// Denominator must be positive:
			// |v0|^2 * |v1|^2 - (v0 . v1)^2 = |v0|^2 * |v1|^2 * (1 - cos(angle)^2) >= 0
			float denominator = GramDeterminant(v0, v1);
			if (denominator < 1.0e-12f)
			{
				// Degenerate triangle, return coordinates along longest edge
				if (d00 > d11)
				{
					GetBaryCentricCoordinates(inA, inB, outU, outV);
					outW = 0.0f;
				}
				else
				{
					GetBaryCentricCoordinates(inA, inC, outU, outW);
					outV = 0.0f;
				}
				return false;
			}
			else
			{
				float a0 = inA.Dot(v0);
				float a1 = inA.Dot(v1);
				outV = (d01 * a1 - d11 * a0) / denominator;
				outW = (d01 * a0 - d00 * a1) / denominator;
				outU = 1.0f - outV - outW;
			}
		}
		else
		{
			// Use v1 and v2 to calculate barycentric coordinates
			float d12 = v1.Dot(v2);

			float denominator = GramDeterminant(v1, v2);
			if (denominator < 1.0e-12f)
			{
				// Degenerate triangle, return coordinates along longest edge
				if (d11 > d22)
				{
					GetBaryCentricCoordinates(inA, inC, outU, outW);
					outV = 0.0f;
				}
				else
				{
					GetBaryCentricCoordinates(inB, inC, outV, outW);
					outU = 0.0f;
				}
				return false;
			}
			else
			{
				float c1 = inC.Dot(v1);
				float c2 = inC.Dot(v2);
				outU = (d22 * c1 - d12 * c2) / denominator;
				outV = (d11 * c2 - d12 * c1) / denominator;
				outW = 1.0f - outU - outV;
			}
		}
		return true;
	}

	/// Compute barycentric coordinates of closest point to origin for hyperplane through (inA, inB, inC, inD)
	/// Point can then be computed as inA * outU + inB * outV + inC * outW + inD * outX
	/// Returns false if the points are degenerate (co-planar in 4D)
	inline bool GetBaryCentricCoordinates(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, Vec4Arg inD, float &outU, float &outV, float &outW, float &outX)
	{
		// Edge vectors from A; solve the 3x3 Gram system to project the origin onto the affine hull
		Vec4 e1 = inB - inA, e2 = inC - inA, e3 = inD - inA;

		// 3x3 Gram matrix G: G_ij = e_i . e_j
		float d00 = e1.Dot(e1), d01 = e1.Dot(e2), d02 = e1.Dot(e3);
		float d11 = e2.Dot(e2), d12 = e2.Dot(e3), d22 = e3.Dot(e3);

		// Cofactors of the 3x3 Gram matrix
		float C00 = d11 * d22 - d12 * d12;
		float C01 = d02 * d12 - d01 * d22;
		float C02 = d01 * d12 - d02 * d11;
		float C11 = d00 * d22 - d02 * d02;
		float C12 = d01 * d02 - d00 * d12;
		float C22 = d00 * d11 - d01 * d01;

		// Gram determinant (always >= 0)
		float det = d00 * C00 + d01 * C01 + d02 * C02;

		if (det < 1.0e-12f)
		{
			// Degenerate tetrahedron, fall back to best triangle using largest Gram sub-determinant
			if (C22 >= C00 && C22 >= C11)
			{
				GetBaryCentricCoordinates(inA, inB, inC, outU, outV, outW);
				outX = 0.0f;
			}
			else if (C11 >= C00)
			{
				GetBaryCentricCoordinates(inA, inB, inD, outU, outV, outX);
				outW = 0.0f;
			}
			else
			{
				GetBaryCentricCoordinates(inA, inC, inD, outU, outW, outX);
				outV = 0.0f;
			}
			return false;
		}

		// RHS: b_i = -A . e_i
		float b0 = -inA.Dot(e1), b1 = -inA.Dot(e2), b2 = -inA.Dot(e3);

		// Cramer's rule: lambda = adj(G) * b / det(G)
		float inv_det = 1.0f / det;
		outV = (b0 * C00 + b1 * C01 + b2 * C02) * inv_det;
		outW = (b0 * C01 + b1 * C11 + b2 * C12) * inv_det;
		outX = (b0 * C02 + b1 * C12 + b2 * C22) * inv_det;
		outU = 1.0f - outV - outW - outX;
		return true;
	}

	/// Get the closest point to the origin of line (inA, inB)
	/// outSet describes which features are closest: 1 = a, 2 = b, 3 = line segment ab
	inline Vec4	GetClosestPointOnLine(Vec4Arg inA, Vec4Arg inB, uint32 &outSet)
	{
		float u, v;
		GetBaryCentricCoordinates(inA, inB, u, v);
		if (v <= 0.0f)
		{
			// inA is closest point
			outSet = 0b0001;
			return inA;
		}
		else if (u <= 0.0f)
		{
			// inB is closest point
			outSet = 0b0010;
			return inB;
		}
		else
		{
			// Closest point lies on line inA inB
			outSet = 0b0011;
			return u * inA + v * inB;
		}
	}

	/// Get the closest point to the origin of triangle (inA, inB, inC)
	/// outSet describes which features are closest: 1 = a, 2 = b, 4 = c, 5 = line segment ac, 7 = triangle interior etc.
	/// If MustIncludeC is true, the function assumes that C is part of the closest feature (vertex, edge, face) and does less work, if the assumption is not true then a closest point to the other features is returned.
	template <bool MustIncludeC = false>
	inline Vec4	GetClosestPointOnTriangle(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, uint32 &outSet)
	{
		// Taken from: Real-Time Collision Detection - Christer Ericson (Section: Closest Point on Triangle to Point)
		// With p = 0

		// The most accurate normal is calculated by using the two shortest edges
		// See: https://box2d.org/posts/2014/01/troublesome-triangle/
		// The difference in normals is most pronounced when one edge is much smaller than the others (in which case the other 2 must have roughly the same length).
		// Therefore we can suffice by just picking the shortest from 2 edges and use that with the 3rd edge to calculate the normal.
		// We first check which of the edges is shorter and if bc is shorter than ac then we swap a with c to a is always on the shortest edge
		UVec4 swap_ac;
		{
			Vec4 ac = inC - inA;
			Vec4 bc = inC - inB;
			swap_ac = Vec4::sLess(bc.DotV(bc), ac.DotV(ac));
		}
		Vec4 a = Vec4::sSelect(inA, inC, swap_ac);
		Vec4 c = Vec4::sSelect(inC, inA, swap_ac);

		// Calculate edges and Gram determinant (dimension-agnostic degeneracy check)
		Vec4 ab = inB - a;
		Vec4 ac = c - a;
		float n_len_sq = GramDeterminant(ab, ac); // = |ab x ac|^2 in 3D

		// Check degenerate
		if (n_len_sq < 1.0e-10f) // Square(FLT_EPSILON) was too small and caused numerical problems, see test case TestCollideParallelTriangleVsCapsule
		{
			// Degenerate, fallback to vertices and edges

			// Start with vertex C being the closest
			uint32 closest_set = 0b0100;
			Vec4 closest_point = inC;
			float best_dist_sq = inC.LengthSq();

			// If the closest point must include C then A or B cannot be closest
			// Note that we test vertices first because we want to prefer a closest vertex over a closest edge (this results in an outSet with fewer bits set)
			if constexpr (!MustIncludeC)
			{
				// Try vertex A
				float a_len_sq = inA.LengthSq();
				if (a_len_sq < best_dist_sq)
				{
					closest_set = 0b0001;
					closest_point = inA;
					best_dist_sq = a_len_sq;
				}

				// Try vertex B
				float b_len_sq = inB.LengthSq();
				if (b_len_sq < best_dist_sq)
				{
					closest_set = 0b0010;
					closest_point = inB;
					best_dist_sq = b_len_sq;
				}
			}

			// Edge AC
			float ac_len_sq = ac.LengthSq();
			if (ac_len_sq > Square(FLT_EPSILON))
			{
				float v = Clamp(-a.Dot(ac) / ac_len_sq, 0.0f, 1.0f);
				Vec4 q = a + v * ac;
				float dist_sq = q.LengthSq();
				if (dist_sq < best_dist_sq)
				{
					closest_set = 0b0101;
					closest_point = q;
					best_dist_sq = dist_sq;
				}
			}

			// Edge BC
			Vec4 bc = inC - inB;
			float bc_len_sq = bc.LengthSq();
			if (bc_len_sq > Square(FLT_EPSILON))
			{
				float v = Clamp(-inB.Dot(bc) / bc_len_sq, 0.0f, 1.0f);
				Vec4 q = inB + v * bc;
				float dist_sq = q.LengthSq();
				if (dist_sq < best_dist_sq)
				{
					closest_set = 0b0110;
					closest_point = q;
					best_dist_sq = dist_sq;
				}
			}

			// If the closest point must include C then AB cannot be closest
			if constexpr (!MustIncludeC)
			{
				// Edge AB
				ab = inB - inA;
				float ab_len_sq = ab.LengthSq();
				if (ab_len_sq > Square(FLT_EPSILON))
				{
					float v = Clamp(-inA.Dot(ab) / ab_len_sq, 0.0f, 1.0f);
					Vec4 q = inA + v * ab;
					float dist_sq = q.LengthSq();
					if (dist_sq < best_dist_sq)
					{
						closest_set = 0b0011;
						closest_point = q;
						best_dist_sq = dist_sq;
					}
				}
			}

			outSet = closest_set;
			return closest_point;
		}

		// Check if P in vertex region outside A
		Vec4 ap = -a;
		float d1 = ab.Dot(ap);
		float d2 = ac.Dot(ap);
		if (d1 <= 0.0f && d2 <= 0.0f)
		{
			outSet = swap_ac.GetX()? 0b0100 : 0b0001;
			return a; // barycentric coordinates (1,0,0)
		}

		// Check if P in vertex region outside B
		Vec4 bp = -inB;
		float d3 = ab.Dot(bp);
		float d4 = ac.Dot(bp);
		if (d3 >= 0.0f && d4 <= d3)
		{
			outSet = 0b0010;
			return inB; // barycentric coordinates (0,1,0)
		}

		// Check if P in edge region of AB, if so return projection of P onto AB
		if (d1 * d4 <= d3 * d2 && d1 >= 0.0f && d3 <= 0.0f)
		{
			float v = d1 / (d1 - d3);
			outSet = swap_ac.GetX()? 0b0110 : 0b0011;
			return a + v * ab; // barycentric coordinates (1-v,v,0)
		}

		// Check if P in vertex region outside C
		Vec4 cp = -c;
		float d5 = ab.Dot(cp);
		float d6 = ac.Dot(cp);
		if (d6 >= 0.0f && d5 <= d6)
		{
			outSet = swap_ac.GetX()? 0b0001 : 0b0100;
			return c; // barycentric coordinates (0,0,1)
		}

		// Check if P in edge region of AC, if so return projection of P onto AC
		if (d5 * d2 <= d1 * d6 && d2 >= 0.0f && d6 <= 0.0f)
		{
			float w = d2 / (d2 - d6);
			outSet = 0b0101;
			return a + w * ac; // barycentric coordinates (1-w,0,w)
		}

		// Check if P in edge region of BC, if so return projection of P onto BC
		float d4_d3 = d4 - d3;
		float d5_d6 = d5 - d6;
		if (d3 * d6 <= d5 * d4 && d4_d3 >= 0.0f && d5_d6 >= 0.0f)
		{
			float w = d4_d3 / (d4_d3 + d5_d6);
			outSet = swap_ac.GetX()? 0b0011 : 0b0110;
			return inB + w * (c - inB); // barycentric coordinates (0,1-w,w)
		}

		// P inside face region.
		// Like 3D Jolt we deviate from Christer Ericson's article here to improve accuracy: solving for
		// barycentric coordinates and then rebuilding the point from them is much less accurate,
		// because the 2x2 Gram solve subtracts two large nearly equal products.
		//
		// TODO(4D): 3D Jolt projects onto the single triangle normal n. A 4D triangle spans a 2-flat
		// whose orthogonal complement is a 2-plane, so there is no single normal - but the wedge
		// bivector B = ab ^ ac plays the same role. For a simple bivector, contracting twice is
		// (minus) the projection onto its plane, |B|^2 * P_plane(v) = -B _| (B _| v), so the closest
		// point on the affine hull is the part of a perpendicular to the plane:
		//     closest = a - P_plane(a) = a + B _| (B _| a) / |B|^2
		// This keeps every product small (it is the same structure as the 3D cross product form) and
		// is exact for the TestLongTriangle sliver, where the barycentric form lost ~4 digits.
		outSet = 0b0111;
		Bivec plane = Bivec::sWedge(ab, ac);
		return a + plane.Contract(plane.Contract(a)) * (1.0f / n_len_sq);
	}

	/// Check if the origin is outside the face (inA, inB, inC) of tetrahedron (inA, inB, inC, inD).
	/// Dimension-agnostic: uses Gram determinant / barycentric coordinates.
	inline bool OriginOutsideOfPlane(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, Vec4Arg inD)
	{
		// Express -A in the basis {e1=B-A, e2=C-A, e3=D-A} by solving the Gram system.
		// Origin is outside face ABC if the barycentric coordinate for D (lambda3) is negative.
		Vec4 e1 = inB - inA, e2 = inC - inA, e3 = inD - inA;

		float d00 = e1.Dot(e1), d01 = e1.Dot(e2), d02 = e1.Dot(e3);
		float d11 = e2.Dot(e2), d12 = e2.Dot(e3), d22 = e3.Dot(e3);

		float b0 = -inA.Dot(e1), b1 = -inA.Dot(e2), b2 = -inA.Dot(e3);

		// Cofactor for the (3,3) minor (used in both det_G and det_G3)
		float c01_12 = d01 * d12 - d11 * d02;

		float det_G = d00 * (d11 * d22 - d12 * d12) - d01 * (d01 * d22 - d12 * d02) + d02 * c01_12;
		float det_G3 = d00 * (d11 * b2 - d12 * b1) - d01 * (d01 * b2 - d02 * b1) + b0 * c01_12;

		// lambda3 = det_G3 / det_G, and det_G >= 0 (it is a Gram determinant), so the origin is
		// outside face ABC exactly when det_G3 < 0.
		//
		// TODO(4D): compare det_G3 against a *relative* slop, not an absolute FLT_EPSILON. Both
		// determinants scale as the 6th power of the simplex size, so testing det_G3 * det_G against
		// an absolute epsilon silently reports "outside" for every face once the simplex shrinks -
		// (1e-2)^6 is already below FLT_EPSILON. That matters far more in 4D than in 3D: a 4D simplex
		// needs 5 points to enclose the origin, so it goes through more (and smaller) intermediate
		// configurations, and once every face reads "outside" the tetrahedron case can never return
		// the full set and GJK can never grow to a pentachoron.
		if (det_G <= 1.0e-6f * d00 * d11 * d22)
			return true; // Degenerate: no reliable inside
		return det_G3 < FLT_EPSILON * det_G;
	}

	/// Returns for each of the planes of the tetrahedron if the origin is outside it.
	/// Dimension-agnostic: uses barycentric coordinates via Gram matrix / Cramer's rule.
	/// Roughly equivalent to:
	///	[OriginOutsideOfPlane(inA, inB, inC, inD),
	///	 OriginOutsideOfPlane(inA, inC, inD, inB),
	///	 OriginOutsideOfPlane(inA, inD, inB, inC),
	///	 OriginOutsideOfPlane(inB, inD, inC, inA)]
	inline UVec4 OriginOutsideOfTetrahedronPlanes(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, Vec4Arg inD)
	{
		// Edge vectors from A
		Vec4 ab = inB - inA;
		Vec4 ac = inC - inA;
		Vec4 ad = inD - inA;

		// Gram matrix G of {ab, ac, ad}: G_ij = e_i . e_j (symmetric, positive semi-definite)
		float d00 = ab.Dot(ab), d01 = ab.Dot(ac), d02 = ab.Dot(ad);
		float d11 = ac.Dot(ac), d12 = ac.Dot(ad), d22 = ad.Dot(ad);

		// RHS: b_i = -A . e_i (projection of -A onto each edge)
		float b0 = -inA.Dot(ab), b1 = -inA.Dot(ac), b2 = -inA.Dot(ad);

		// Cofactors of the symmetric 3x3 Gram matrix (cofactor matrix is also symmetric)
		float C00 = d11 * d22 - d12 * d12;
		float C01 = d02 * d12 - d01 * d22;
		float C02 = d01 * d12 - d02 * d11;
		float C11 = d00 * d22 - d02 * d02;
		float C12 = d01 * d02 - d00 * d12;
		float C22 = d00 * d11 - d01 * d01;

		// Gram determinant (always >= 0; zero iff edges are linearly dependent = degenerate)
		float det_G = d00 * C00 + d01 * C01 + d02 * C02;

		// Cramer's rule: barycentric coordinate numerators = adj(G) * b
		float det_G1 = b0 * C00 + b1 * C01 + b2 * C02;  // lambda_B numerator
		float det_G2 = b0 * C01 + b1 * C11 + b2 * C12;  // lambda_C numerator
		float det_G3 = b0 * C02 + b1 * C12 + b2 * C22;  // lambda_D numerator

		// lambda_A = 1 - lambda_B - lambda_C - lambda_D
		float det_GA = det_G - det_G1 - det_G2 - det_G3;

		// Origin is outside face i iff the opposite vertex's barycentric coordinate is negative.
		// Since det_G >= 0 (Gram), lambda_i < 0 iff det_Gi < 0.
		//
		// TODO(4D): the slop has to be relative to det_G - see the note in OriginOutsideOfPlane. The
		// previous form compared det_Gi * det_G against an absolute FLT_EPSILON, which reports every
		// face as "outside" as soon as the simplex is small, so the tetrahedron never returns the full
		// set and GJK never grows to 5 points. A degenerate simplex (det_G == 0) still reports all
		// faces outside, which is the safe answer.
		// A numerically degenerate simplex has no reliable inside: report every face outside so the
		// caller falls back to the sub-simplices. det_G <= d00*d11*d22 (Hadamard), so the ratio is a
		// scale-invariant conditioning measure.
		if (det_G <= 1.0e-6f * d00 * d11 * d22)
			return UVec4::sReplicate(0xffffffff);

		// Order: [face ABC/D opposite, face ACD/B opposite, face ADB/C opposite, face BDC/A opposite]
		Vec4 numerators(det_G3, det_G1, det_G2, det_GA);
		return Vec4::sLess(numerators, Vec4::sReplicate(FLT_EPSILON * det_G));
	}

	/// Get the closest point between tetrahedron (inA, inB, inC, inD) to the origin
	/// outSet specifies which feature was closest, 1 = a, 2 = b, 4 = c, 8 = d. Edges have 2 bits set, triangles 3 and if the point is in the interior 4 bits are set.
	/// If MustIncludeD is true, the function assumes that D is part of the closest feature (vertex, edge, face, tetrahedron) and does less work, if the assumption is not true then a closest point to the other features is returned.
	template <bool MustIncludeD = false>
	inline Vec4	GetClosestPointOnTetrahedron(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, Vec4Arg inD, uint32 &outSet)
	{
		// Taken from: Real-Time Collision Detection - Christer Ericson (Section: Closest Point on Tetrahedron to Point)
		// With p = 0

		// Start out assuming the origin projects into the tetrahedron.
		// TODO(4D): in 3D a tetrahedron is a solid and "inside all 4 halfspaces" means the closest
		// point IS the origin. In 4D four points span only a 3-flat, so the tetrahedron has no
		// 4-volume and the origin is generally off that flat: the closest point is the orthogonal
		// projection of the origin onto the affine hull, not the origin itself. (It is the
		// pentachoron - 5 points - that can actually contain the origin in 4D, see
		// GetClosestPointOnPentachoron.) GetBaryCentricCoordinates solves exactly that projection
		// via the 3x3 Gram system, and fills the coordinates even in the degenerate case.
		uint32 closest_set = 0b1111;
		float bary_u, bary_v, bary_w, bary_x;
		GetBaryCentricCoordinates(inA, inB, inC, inD, bary_u, bary_v, bary_w, bary_x);
		Vec4 closest_point = bary_u * inA + bary_v * inB + bary_w * inC + bary_x * inD;
		float best_dist_sq = FLT_MAX;

		// Determine for each of the faces of the tetrahedron if the origin is in front of the plane
		UVec4 origin_out_of_planes = OriginOutsideOfTetrahedronPlanes(inA, inB, inC, inD);

		// If point outside face abc then compute closest point on abc
		if (origin_out_of_planes.GetX()) // OriginOutsideOfPlane(inA, inB, inC, inD)
		{
			if constexpr (MustIncludeD)
			{
				// If the closest point must include D then ABC cannot be closest but the closest point
				// cannot be an interior point either so we return A as closest point
				closest_set = 0b0001;
				closest_point = inA;
			}
			else
			{
				// Test the face normally
				closest_point = GetClosestPointOnTriangle<false>(inA, inB, inC, closest_set);
			}
			best_dist_sq = closest_point.LengthSq();
		}

		// Repeat test for face acd
		if (origin_out_of_planes.GetY()) // OriginOutsideOfPlane(inA, inC, inD, inB)
		{
			uint32 set;
			Vec4 q = GetClosestPointOnTriangle<MustIncludeD>(inA, inC, inD, set);
			float dist_sq = q.LengthSq();
			if (dist_sq < best_dist_sq)
			{
				best_dist_sq = dist_sq;
				closest_point = q;
				closest_set = (set & 0b0001) + ((set & 0b0110) << 1);
			}
		}

		// Repeat test for face adb
		if (origin_out_of_planes.GetZ()) // OriginOutsideOfPlane(inA, inD, inB, inC)
		{
			// Keep original vertex order, it doesn't matter if the triangle is facing inward or outward
			// and it improves consistency for GJK which will always add a new vertex D and keep the closest
			// feature from the previous iteration in ABC
			uint32 set;
			Vec4 q = GetClosestPointOnTriangle<MustIncludeD>(inA, inB, inD, set);
			float dist_sq = q.LengthSq();
			if (dist_sq < best_dist_sq)
			{
				best_dist_sq = dist_sq;
				closest_point = q;
				closest_set = (set & 0b0011) + ((set & 0b0100) << 1);
			}
		}

		// Repeat test for face bdc
		if (origin_out_of_planes.GetW()) // OriginOutsideOfPlane(inB, inD, inC, inA)
		{
			// Keep original vertex order, it doesn't matter if the triangle is facing inward or outward
			// and it improves consistency for GJK which will always add a new vertex D and keep the closest
			// feature from the previous iteration in ABC
			uint32 set;
			Vec4 q = GetClosestPointOnTriangle<MustIncludeD>(inB, inC, inD, set);
			float dist_sq = q.LengthSq();
			if (dist_sq < best_dist_sq)
			{
				closest_point = q;
				closest_set = set << 1;
			}
		}

		outSet = closest_set;
		return closest_point;
	}

	/// Returns a bitmask for each of the 5 tetrahedral faces of the pentachoron if the origin is outside it.
	/// Uses barycentric coordinates via 4x4 Gram matrix / Cramer's rule.
	/// Bit 0: face ABCD (opposite E), Bit 1: face BCDE (opposite A),
	/// Bit 2: face ACDE (opposite B), Bit 3: face ABDE (opposite C), Bit 4: face ABCE (opposite D)
	inline uint32 OriginOutsideOfPentachoronPlanes(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, Vec4Arg inD, Vec4Arg inE)
	{
		// Edge vectors from A
		Vec4 e1 = inB - inA, e2 = inC - inA, e3 = inD - inA, e4 = inE - inA;

		// 4x4 Gram matrix G: G_ij = e_i . e_j (symmetric, 10 unique entries)
		float d00 = e1.Dot(e1), d01 = e1.Dot(e2), d02 = e1.Dot(e3), d03 = e1.Dot(e4);
		float d11 = e2.Dot(e2), d12 = e2.Dot(e3), d13 = e2.Dot(e4);
		float d22 = e3.Dot(e3), d23 = e3.Dot(e4);
		float d33 = e4.Dot(e4);

		// RHS: b_i = -A . e_i
		float b0 = -inA.Dot(e1), b1 = -inA.Dot(e2), b2 = -inA.Dot(e3), b3 = -inA.Dot(e4);

		// 2x2 sub-determinants from rows 2,3 of G (6 values)
		float s23_01 = d02 * d13 - d12 * d03;
		float s23_02 = d02 * d23 - d22 * d03;
		float s23_03 = d02 * d33 - d23 * d03;
		float s23_12 = d12 * d23 - d22 * d13;
		float s23_13 = d12 * d33 - d23 * d13;
		float s23_23 = d22 * d33 - d23 * d23;

		// 3x3 cofactors of G — first row (C0j) and second row (C1j)
		float C00 =  d11 * s23_23 - d12 * s23_13 + d13 * s23_12;
		float C01 = -(d01 * s23_23 - d12 * s23_03 + d13 * s23_02);
		float C02 =  d01 * s23_13 - d11 * s23_03 + d13 * s23_01;
		float C03 = -(d01 * s23_12 - d11 * s23_02 + d12 * s23_01);

		float C11 =  d00 * s23_23 - d02 * s23_03 + d03 * s23_02;
		float C12 = -(d00 * s23_13 - d01 * s23_03 + d03 * s23_01);
		float C13 =  d00 * s23_12 - d01 * s23_02 + d02 * s23_01;

		// 2x2 sub-determinants from rows 1,3 of G (5 values)
		float s13_01 = d01 * d13 - d11 * d03;
		float s13_02 = d01 * d23 - d12 * d03;
		float s13_03 = d01 * d33 - d13 * d03;
		float s13_12 = d11 * d23 - d12 * d13;
		float s13_13 = d11 * d33 - d13 * d13;

		// Third row cofactors
		float C22 =  d00 * s13_13 - d01 * s13_03 + d03 * s13_01;
		float C23 = -(d00 * s13_12 - d01 * s13_02 + d02 * s13_01);

		// 2x2 sub-determinants from rows 1,2 of G (3 values)
		float s12_01 = d01 * d12 - d11 * d02;
		float s12_02 = d01 * d22 - d12 * d02;
		float s12_12 = d11 * d22 - d12 * d12;

		// Fourth row cofactor
		float C33 = d00 * s12_12 - d01 * s12_02 + d02 * s12_01;

		// Gram determinant (expansion along first row)
		float det_G = d00 * C00 + d01 * C01 + d02 * C02 + d03 * C03;

		// Cramer's rule: barycentric coordinate numerators = adj(G) * b
		float num_B = b0 * C00 + b1 * C01 + b2 * C02 + b3 * C03;
		float num_C = b0 * C01 + b1 * C11 + b2 * C12 + b3 * C13;
		float num_D = b0 * C02 + b1 * C12 + b2 * C22 + b3 * C23;
		float num_E = b0 * C03 + b1 * C13 + b2 * C23 + b3 * C33;
		float num_A = det_G - num_B - num_C - num_D - num_E;

		// Origin is outside the face opposite vertex V iff lambda_V < 0, and det_G >= 0 (Gram), so
		// iff num_V < 0. TODO(4D): the slop must be relative to det_G - see OriginOutsideOfPlane. An
		// absolute FLT_EPSILON on num_V * det_G reports every face as "outside" once the simplex is
		// small (both scale as a high power of the simplex size), which stops the pentachoron from
		// ever reporting that it encloses the origin. A degenerate simplex (det_G == 0) still reports
		// all faces outside, which is the safe answer.
		// Degenerate simplex: report every face outside (see OriginOutsideOfTetrahedronPlanes)
		if (det_G <= 1.0e-6f * d00 * d11 * d22 * d33)
			return 0x1f;

		float slop = FLT_EPSILON * det_G;
		uint32 result = 0;
		if (num_E < slop) result |= 0x01;	// face ABCD, opposite E
		if (num_A < slop) result |= 0x02;	// face BCDE, opposite A
		if (num_B < slop) result |= 0x04;	// face ACDE, opposite B
		if (num_C < slop) result |= 0x08;	// face ABDE, opposite C
		if (num_D < slop) result |= 0x10;	// face ABCE, opposite D
		return result;
	}

	/// Get the closest point between pentachoron (inA, inB, inC, inD, inE) and the origin.
	/// outSet specifies which feature was closest: bit 0 = A, bit 1 = B, bit 2 = C, bit 3 = D, bit 4 = E.
	/// If MustIncludeE is true, the function assumes that E is part of the closest feature and does less work.
	template <bool MustIncludeE = false>
	inline Vec4	GetClosestPointOnPentachoron(Vec4Arg inA, Vec4Arg inB, Vec4Arg inC, Vec4Arg inD, Vec4Arg inE, uint32 &outSet)
	{
		// Start out assuming point inside pentachoron, so closest to itself
		uint32 closest_set = 0b11111;
		Vec4 closest_point = Vec4::sZero();
		float best_dist_sq = FLT_MAX;

		// Determine for each of the 5 tetrahedral faces if the origin is outside
		uint32 origin_out_of_planes = OriginOutsideOfPentachoronPlanes(inA, inB, inC, inD, inE);

		// Face ABCD (opposite E) — bit 0
		if (origin_out_of_planes & 0x01)
		{
			if constexpr (MustIncludeE)
			{
				// E must be in closest feature, so ABCD cannot be closest.
				// Return A as a fallback (same pattern as tetrahedron MustIncludeD).
				closest_set = 0b00001;
				closest_point = inA;
			}
			else
			{
				closest_point = GetClosestPointOnTetrahedron<false>(inA, inB, inC, inD, closest_set);
			}
			best_dist_sq = closest_point.LengthSq();
		}

		// Face BCDE (opposite A) — bit 1
		if (origin_out_of_planes & 0x02)
		{
			uint32 set;
			Vec4 q = GetClosestPointOnTetrahedron<MustIncludeE>(inB, inC, inD, inE, set);
			float dist_sq = q.LengthSq();
			if (dist_sq < best_dist_sq)
			{
				best_dist_sq = dist_sq;
				closest_point = q;
				// Tet bits: 0=B,1=C,2=D,3=E → pent bits: 1=B,2=C,3=D,4=E
				closest_set = set << 1;
			}
		}

		// Face ACDE (opposite B) — bit 2
		if (origin_out_of_planes & 0x04)
		{
			uint32 set;
			Vec4 q = GetClosestPointOnTetrahedron<MustIncludeE>(inA, inC, inD, inE, set);
			float dist_sq = q.LengthSq();
			if (dist_sq < best_dist_sq)
			{
				best_dist_sq = dist_sq;
				closest_point = q;
				// Tet bits: 0=A,1=C,2=D,3=E → pent bits: 0=A,2=C,3=D,4=E
				closest_set = (set & 0b0001) | ((set & 0b1110) << 1);
			}
		}

		// Face ABDE (opposite C) — bit 3
		if (origin_out_of_planes & 0x08)
		{
			uint32 set;
			Vec4 q = GetClosestPointOnTetrahedron<MustIncludeE>(inA, inB, inD, inE, set);
			float dist_sq = q.LengthSq();
			if (dist_sq < best_dist_sq)
			{
				best_dist_sq = dist_sq;
				closest_point = q;
				// Tet bits: 0=A,1=B,2=D,3=E → pent bits: 0=A,1=B,3=D,4=E
				closest_set = (set & 0b0011) | ((set & 0b1100) << 1);
			}
		}

		// Face ABCE (opposite D) — bit 4
		if (origin_out_of_planes & 0x10)
		{
			uint32 set;
			Vec4 q = GetClosestPointOnTetrahedron<MustIncludeE>(inA, inB, inC, inE, set);
			float dist_sq = q.LengthSq();
			if (dist_sq < best_dist_sq)
			{
				closest_point = q;
				// Tet bits: 0=A,1=B,2=C,3=E → pent bits: 0=A,1=B,2=C,4=E
				closest_set = (set & 0b0111) | ((set & 0b1000) << 1);
			}
		}

		outSet = closest_set;
		return closest_point;
	}
};

JPH_PRECISE_MATH_OFF

JPH_NAMESPACE_END
