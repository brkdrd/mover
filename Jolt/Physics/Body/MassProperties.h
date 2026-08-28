// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/ObjectStream/SerializableObject.h>
#include <Jolt/Math/Bivec.h>

JPH_NAMESPACE_BEGIN

class StreamIn;
class StreamOut;

/// Describes the mass and inertia properties of a 4D body. Used during body construction only.
///
/// In 4D rigid body dynamics the inertia tensor is a 6x6 linear map bivector -> bivector.
/// We store the underlying 4x4 symmetric second moment tensor
/// \f$M_{ij} = \int \rho(\mathbf{x}) \, x_i x_j \, dV\f$
/// from which the 6 bivector principal moments are derived (in the principal frame) as:
///   \f$I_{ij} = M_{ii} + M_{jj}\f$
/// for each of the 6 rotation planes (ij \in {12, 13, 14, 23, 24, 34}).
///
/// This representation is compact (10 unique values vs 21 for the full 6x6 inertia)
/// and generalizes cleanly under rotation (M' = R M R^T), translation (parallel axis theorem)
/// and scale (M'_ij = s_i s_j M_ij).
///
/// Note: mInertia is stored in a Mat44 but is NOT an orthogonal rotation matrix.
/// It is a 4x4 symmetric matrix (4 Lane4 columns). Only Mat44's linear-algebra operations
/// (add, scale, transpose, inverse, multiply) are used on it — not sRotation or similar helpers.
class JPH_EXPORT MassProperties
{
	JPH_DECLARE_SERIALIZABLE_NON_VIRTUAL(JPH_EXPORT, MassProperties)

public:
	/// Eigendecomposition of the 4x4 symmetric second moment tensor M:
	/// finds an orthogonal R (SO(4)) and a diagonal D such that \f$R \: D \: R^{T} = M\f$.
	/// @param outRotation The 4x4 rotation matrix R (columns are the principal axes)
	/// @param outDiagonal The 4 principal second moments \f$(M_1, M_2, M_3, M_4)\f$, sorted largest first
	/// @return True if successful, false if eigendecomposition failed to converge
	bool					DecomposePrincipalMomentsOfInertia(Mat44 &outRotation, Vec4 &outDiagonal) const;

	/// Compute the 6 bivector principal moments from a 4D second moment diagonal.
	/// In the principal frame the bivector inertia is diagonal: \f$I_{ij} = M_{ii} + M_{jj}\f$.
	static Bivec			sBivecInertiaFromSecondMomentDiagonal(Vec4Arg inDiagonal);

	/// Set the mass and second moment tensor of a solid 4D box (tesseract) with edge sizes inBoxSize and uniform density.
	/// For a solid tesseract: mass = \f$\rho \prod_i S_i\f$, and M is diagonal with \f$M_{ii} = m S_i^2 / 12\f$.
	void					SetMassAndInertiaOfSolidBox(Vec4Arg inBoxSize, float inDensity);

	/// Scale the mass to the given value and scale the inertia proportionally.
	void					ScaleToMass(float inMass);

	/// Calculate the edge sizes of a solid 4D box that would produce the given second moment diagonal.
	/// Inverse of SetMassAndInertiaOfSolidBox: \f$S_i = \sqrt{12 M_i / m}\f$.
	static Vec4				sGetEquivalentSolidBoxSize(float inMass, Vec4Arg inInertiaDiagonal);

	/// Rotate the second moment tensor by a 4x4 orthogonal rotation matrix: \f$M' = R \, M \, R^T\f$
	void					Rotate(Mat44Arg inRotation);

	/// Shift the body by inTranslation (origin → origin + t).
	/// Parallel axis theorem for the second moment tensor (when the original COM is at the origin):
	/// \f$M'_{ij} = M_{ij} + m \, t_i \, t_j\f$
	void					Translate(Vec4Arg inTranslation);

	/// Scale the body by a 4D scale factor.
	/// The second moment tensor scales as \f$M'_{ij} = s_i \, s_j \, M_{ij}\f$,
	/// and mass scales linearly with 4D hypervolume: \f$m' = |s_1 s_2 s_3 s_4| \, m\f$.
	/// Individual scale components may be negative (to flip axes) — mass uses absolute value.
	void					Scale(Vec4Arg inScale);

	/// Saves the state of this object in binary form to inStream.
	void					SaveBinaryState(StreamOut &inStream) const;

	/// Restore the state of this object from inStream.
	void					RestoreBinaryState(StreamIn &inStream);

	/// Mass of the shape (kg)
	float					mMass = 0.0f;

	/// Second moment tensor \f$M_{ij} = \int \rho(\mathbf{x}) \, x_i x_j \, dV\f$ (kg m^2), 4x4 symmetric.
	/// Stored in a Mat44 for SIMD-friendly layout — NOT an orthogonal matrix.
	Mat44					mInertia = Mat44::sZero();
};

JPH_NAMESPACE_END
