// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Collision/Shape/CapsuleShape.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>
#include <Jolt/Physics/Collision/Shape/ScaleHelpers.h>
#include <Jolt/Physics/Collision/RayCast.h>
#include <Jolt/Physics/Collision/CastResult.h>
#include <Jolt/Physics/Collision/CollidePointResult.h>
#include <Jolt/Physics/Collision/TransformedShape.h>
#include <Jolt/Physics/Collision/CollideSoftBodyVertexIterator.h>
#include <Jolt/Geometry/RayCapsule.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>
#include <Jolt/Core/StreamIn.h>
#include <Jolt/Core/StreamOut.h>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_DEBUG_RENDERER

JPH_NAMESPACE_BEGIN

JPH_IMPLEMENT_SERIALIZABLE_VIRTUAL(CapsuleShapeSettings)
{
	JPH_ADD_BASE_CLASS(CapsuleShapeSettings, ConvexShapeSettings)

	JPH_ADD_ATTRIBUTE(CapsuleShapeSettings, mRadius)
	JPH_ADD_ATTRIBUTE(CapsuleShapeSettings, mHalfHeightOfCylinder)
}

ShapeSettings::ShapeResult CapsuleShapeSettings::Create() const
{
	if (mCachedResult.IsEmpty())
	{
		Ref<Shape> shape;
		if (IsValid() && IsSphere())
		{
			// If the capsule has no height, use a sphere instead
			shape = new SphereShape(mRadius, mMaterial);
			mCachedResult.Set(shape);
		}
		else
			shape = new CapsuleShape(*this, mCachedResult);
	}
	return mCachedResult;
}

CapsuleShape::CapsuleShape(const CapsuleShapeSettings &inSettings, ShapeResult &outResult) :
	ConvexShape(EShapeSubType::Capsule, inSettings, outResult),
	mRadius(inSettings.mRadius),
	mHalfHeightOfCylinder(inSettings.mHalfHeightOfCylinder)
{
	if (inSettings.mHalfHeightOfCylinder <= 0.0f)
	{
		outResult.SetError("Invalid height");
		return;
	}

	if (inSettings.mRadius <= 0.0f)
	{
		outResult.SetError("Invalid radius");
		return;
	}

	outResult.Set(this);
}

class CapsuleShape::CapsuleNoConvex final : public Support
{
public:
					CapsuleNoConvex(Vec4Arg inHalfHeightOfCylinder, float inConvexRadius) :
		mHalfHeightOfCylinder(inHalfHeightOfCylinder),
		mConvexRadius(inConvexRadius)
	{
		static_assert(sizeof(CapsuleNoConvex) <= sizeof(SupportBuffer), "Buffer size too small");
		JPH_ASSERT(IsAligned(this, alignof(CapsuleNoConvex)));
	}

	virtual Vec4	GetSupport(Vec4Arg inDirection) const override
	{
		if (inDirection.GetY() > 0)
			return mHalfHeightOfCylinder;
		else
			return -mHalfHeightOfCylinder;
	}

	virtual float	GetConvexRadius() const override
	{
		return mConvexRadius;
	}

private:
	Vec4			mHalfHeightOfCylinder;
	float			mConvexRadius;
};

class CapsuleShape::CapsuleWithConvex final : public Support
{
public:
					CapsuleWithConvex(Vec4Arg inHalfHeightOfCylinder, float inRadius) :
		mHalfHeightOfCylinder(inHalfHeightOfCylinder),
		mRadius(inRadius)
	{
		static_assert(sizeof(CapsuleWithConvex) <= sizeof(SupportBuffer), "Buffer size too small");
		JPH_ASSERT(IsAligned(this, alignof(CapsuleWithConvex)));
	}

	virtual Vec4	GetSupport(Vec4Arg inDirection) const override
	{
		float len = inDirection.Length();
		Vec4 radius = len > 0.0f? Vec4(inDirection * (mRadius / len)) : Vec4::sZero();

		if (inDirection.GetY() > 0)
			return radius + mHalfHeightOfCylinder;
		else
			return radius - mHalfHeightOfCylinder;
	}

	virtual float	GetConvexRadius() const override
	{
		return 0.0f;
	}

private:
	Vec4			mHalfHeightOfCylinder;
	float			mRadius;
};

const ConvexShape::Support *CapsuleShape::GetSupportFunction(ESupportMode inMode, SupportBuffer &inBuffer, Vec4Arg inScale) const
{
	JPH_ASSERT(IsValidScale(inScale));

	// Get scaled capsule
	Vec4 abs_scale = inScale.Abs();
	float scale = abs_scale.GetX();
	Vec4 scaled_half_height_of_cylinder = Vec4(0, scale * mHalfHeightOfCylinder, 0, 0);
	float scaled_radius = scale * mRadius;

	switch (inMode)
	{
	case ESupportMode::IncludeConvexRadius:
		return new (&inBuffer) CapsuleWithConvex(scaled_half_height_of_cylinder, scaled_radius);

	case ESupportMode::ExcludeConvexRadius:
	case ESupportMode::Default:
		return new (&inBuffer) CapsuleNoConvex(scaled_half_height_of_cylinder, scaled_radius);
	}

	JPH_ASSERT(false);
	return nullptr;
}

void CapsuleShape::GetSupportingFace(const SubShapeID &inSubShapeID, Vec4Arg inDirection, Vec4Arg inScale, RMat44Arg inCenterOfMassTransform, SupportingFace &outVertices) const
{
	JPH_ASSERT(inSubShapeID.IsEmpty(), "Invalid subshape ID");
	JPH_ASSERT(IsValidScale(inScale));

	// Get direction in the plane perpendicular to the capsule axis (X, Z, W)
	Vec4 direction = inDirection;
	direction.SetY(0.0f);

	// Check zero vector, in this case we're hitting from top/bottom so there's no supporting face
	float len = direction.Length();
	if (len == 0.0f)
		return;

	// Get scaled capsule
	Vec4 abs_scale = inScale.Abs();
	float scale = abs_scale.GetX();
	Vec4 scaled_half_height_of_cylinder = Vec4(0, scale * mHalfHeightOfCylinder, 0, 0);
	float scaled_radius = scale * mRadius;

	// Get support point for top and bottom sphere in the opposite of 'direction' (including convex radius)
	Vec4 support = Vec4((scaled_radius / len) * direction);
	Vec4 support_top = scaled_half_height_of_cylinder - support;
	Vec4 support_bottom = -scaled_half_height_of_cylinder - support;

	// Get projection on inDirection
	// Note that inDirection is not normalized, so we need to divide by inDirection.Length() to get the actual projection
	// We've multiplied both sides of the if below with inDirection.Length()
	float proj_top = support_top.Dot(inDirection);
	float proj_bottom = support_bottom.Dot(inDirection);

	// If projection is roughly equal then return line, otherwise we return nothing as there's only 1 point
	if (abs(proj_top - proj_bottom) < cCapsuleProjectionSlop * inDirection.Length())
	{
		outVertices.push_back(Vec4(inCenterOfMassTransform * support_top));
		outVertices.push_back(Vec4(inCenterOfMassTransform * support_bottom));
	}
}

MassProperties CapsuleShape::GetMassProperties() const
{
	MassProperties p;

	float density = GetDensity();
	float r = mRadius;
	float H = mHalfHeightOfCylinder;
	float r2 = r * r, r3 = r2 * r, r4 = r2 * r2, r5 = r4 * r, r6 = r4 * r2;
	float H2 = H * H, H3 = H2 * H;
	const float PI = JPH_PI;
	const float PI2 = JPH_PI * JPH_PI;

	// A 4D capsule = a 4D cylinder (3-ball cross section of radius r swept along Y over [-H, H])
	// capped by two half-4-balls of radius r centered at (0, +-H, 0, 0).
	// mInertia stores the raw second-moment tensor M_ij = integral(rho * x_i * x_j dV), matching BoxShape's convention.
	// By symmetry M is diagonal with M_xx = M_zz = M_ww (perpendicular) and M_yy (along the axis).

	// Cylinder: 3-ball volume = 4/3 pi r^3, swept over length 2H
	float mass_cyl = density * (8.0f / 3.0f) * PI * r3 * H;
	// Caps: two hemispheres = one full 4-ball, volume = pi^2/2 r^4
	float mass_caps = density * 0.5f * PI2 * r4;
	p.mMass = mass_cyl + mass_caps;

	// Cylinder second moments:
	//   M_yy = rho * (integral_{-H}^{H} y^2 dy) * (3-ball area)   = rho * (2H^3/3) * (4/3 pi r^3) = rho * 8 pi/9 r^3 H^3
	//   M_xx = rho * (integral_{3-ball} x^2 dA) * (2H)            = rho * (4 pi r^5/15) * (2H)     = rho * 8 pi/15 r^5 H
	float Myy_cyl = density * (8.0f * PI / 9.0f) * r3 * H3;
	float Mxx_cyl = density * (8.0f * PI / 15.0f) * r5 * H;

	// Cap second moments (per hemisphere, in its own centered frame):
	//   V_h    = pi^2 r^4 / 4                (half of the 4-ball volume)
	//   S1     = integral y dV = 4 pi/15 r^5 (first moment about the cut plane)
	//   S2     = integral y^2 dV = pi^2 r^6 / 24
	//   Sperp  = integral x^2 dV = pi^2 r^6 / 24
	// Each cap is offset to y = +-H, so by the parallel-axis theorem for second moments:
	//   M_yy_cap = rho * (H^2 V_h + 2H S1 + S2), and the two caps sum (both give the same value).
	float V_h = 0.25f * PI2 * r4;
	float S1 = (4.0f * PI / 15.0f) * r5;
	float S2 = PI2 * r6 / 24.0f;
	float Sperp = PI2 * r6 / 24.0f;
	float Myy_caps = 2.0f * density * (H2 * V_h + 2.0f * H * S1 + S2);
	float Mxx_caps = 2.0f * density * Sperp;

	float Mxx = Mxx_cyl + Mxx_caps; // = M_zz = M_ww
	float Myy = Myy_cyl + Myy_caps;

	// TODO(4D): validate this second-moment derivation against a brute-force numerical integral
	// (per the project's "never trust, always test the algebra" rule). It reduces correctly to the
	// 4-ball second moment (mass * r^2 / 6 per axis) as H -> 0.
	Mat44 m = Mat44::sZero();
	m(0, 0) = Mxx;
	m(1, 1) = Myy;
	m(2, 2) = Mxx;
	m(3, 3) = Mxx;
	p.mInertia = m;

	return p;
}

Vec4 CapsuleShape::GetSurfaceNormal(const SubShapeID &inSubShapeID, Vec4Arg inLocalSurfacePosition) const
{
	JPH_ASSERT(inSubShapeID.IsEmpty(), "Invalid subshape ID");

	if (inLocalSurfacePosition.GetY() > mHalfHeightOfCylinder)
		return Vec4((inLocalSurfacePosition - Vec4(0, mHalfHeightOfCylinder, 0, 0)).Normalized());
	else if (inLocalSurfacePosition.GetY() < -mHalfHeightOfCylinder)
		return Vec4((inLocalSurfacePosition - Vec4(0, -mHalfHeightOfCylinder, 0, 0)).Normalized());
	else
	{
		// Perpendicular (X, Z, W) part
		Vec4 horiz(inLocalSurfacePosition.GetX(), 0.0f, inLocalSurfacePosition.GetZ(), inLocalSurfacePosition.GetW());
		float len = horiz.Length();
		return len > 0.0f? Vec4(horiz.Normalized()) : Vec4::sAxisX();
	}
}

AABox CapsuleShape::GetLocalBounds() const
{
	Vec4 extent = Vec4::sReplicate(mRadius) + Vec4(0, mHalfHeightOfCylinder, 0, 0);
	return AABox(-extent, extent);
}

AABox CapsuleShape::GetWorldSpaceBounds(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale) const
{
	JPH_ASSERT(IsValidScale(inScale));

	Vec4 abs_scale = inScale.Abs();
	float scale = abs_scale.GetX();
	Vec4 extent = Vec4::sReplicate(scale * mRadius);
	Vec4 height = Vec4(0, scale * mHalfHeightOfCylinder, 0, 0);
	Vec4 p1 = Vec4(inCenterOfMassTransform * -height);
	Vec4 p2 = Vec4(inCenterOfMassTransform * height);
	return AABox(Vec4::sMin(p1, p2) - extent, Vec4::sMax(p1, p2) + extent);
}

#ifdef JPH_DEBUG_RENDERER
void CapsuleShape::Draw(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inUseMaterialColors, bool inDrawWireframe) const
{
	DebugRenderer::EDrawMode draw_mode = inDrawWireframe? DebugRenderer::EDrawMode::Wireframe : DebugRenderer::EDrawMode::Solid;
	inRenderer->DrawCapsule(inCenterOfMassTransform * Mat44::sScale(inScale.Abs().GetX()), mHalfHeightOfCylinder, mRadius, inUseMaterialColors? GetMaterial()->GetDebugColor() : inColor, DebugRenderer::ECastShadow::On, draw_mode);
}
#endif // JPH_DEBUG_RENDERER

bool CapsuleShape::CastRay(const RayCast &inRay, const SubShapeIDCreator &inSubShapeIDCreator, RayCastResult &ioHit) const
{
	// Test ray against capsule
	float fraction = RayCapsule(inRay.mOrigin, inRay.mDirection, mHalfHeightOfCylinder, mRadius);
	if (fraction < ioHit.mFraction)
	{
		ioHit.mFraction = fraction;
		ioHit.mSubShapeID2 = inSubShapeIDCreator.GetID();
		return true;
	}
	return false;
}

void CapsuleShape::CollidePoint(Vec4Arg inPoint, const SubShapeIDCreator &inSubShapeIDCreator, CollidePointCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	float radius_sq = Square(mRadius);

	// Get vertical distance to the top/bottom sphere centers
	float delta_y = abs(inPoint.GetY()) - mHalfHeightOfCylinder;

	// Get squared distance in the plane perpendicular to the capsule axis (X, Z, W)
	float perp_sq = Square(inPoint.GetX()) + Square(inPoint.GetZ()) + Square(inPoint.GetW());

	// Check if the point is in one of the two spheres
	bool in_sphere = perp_sq + Square(delta_y) <= radius_sq;

	// Check if the point is in the cylinder in the middle
	bool in_cylinder = delta_y <= 0.0f && perp_sq <= radius_sq;

	if (in_sphere || in_cylinder)
		ioCollector.AddHit({ TransformedShape::sGetBodyID(ioCollector.GetContext()), inSubShapeIDCreator.GetID() });
}

void CapsuleShape::CollideSoftBodyVertices(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const CollideSoftBodyVertexIterator &inVertices, uint inNumVertices, int inCollidingShapeIndex) const
{
	// TODO(4D): CollideSoftBodyVertexIterator still exposes 3D positions. Stubbed until the soft-body
	// subsystem is migrated to 4D (SoftBody is cut for v1).
	(void)inCenterOfMassTransform; (void)inScale; (void)inVertices; (void)inNumVertices; (void)inCollidingShapeIndex;
}

void CapsuleShape::GetTrianglesStart(GetTrianglesContext &ioContext, const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale) const
{
	// TODO(4D): The GetTriangles vertex pipeline is still 3D. Stubbed until the 4D (tetrahedra-based)
	// boundary tessellation is defined.
	(void)ioContext; (void)inBox; (void)inPositionCOM; (void)inRotation; (void)inScale;
}

int CapsuleShape::GetTrianglesNext(GetTrianglesContext &ioContext, int inMaxTrianglesRequested, Float4 *outTriangleVertices, const PhysicsMaterial **outMaterials) const
{
	// TODO(4D): See GetTrianglesStart. Returns 0 to indicate no tessellation data is currently produced.
	(void)ioContext; (void)inMaxTrianglesRequested; (void)outTriangleVertices; (void)outMaterials;
	return 0;
}

void CapsuleShape::SaveBinaryState(StreamOut &inStream) const
{
	ConvexShape::SaveBinaryState(inStream);

	inStream.Write(mRadius);
	inStream.Write(mHalfHeightOfCylinder);
}

void CapsuleShape::RestoreBinaryState(StreamIn &inStream)
{
	ConvexShape::RestoreBinaryState(inStream);

	inStream.Read(mRadius);
	inStream.Read(mHalfHeightOfCylinder);
}

bool CapsuleShape::IsValidScale(Vec4Arg inScale) const
{
	return ConvexShape::IsValidScale(inScale) && ScaleHelpers::IsUniformScale(inScale.Abs());
}

Vec4 CapsuleShape::MakeScaleValid(Vec4Arg inScale) const
{
	Vec4 scale = ScaleHelpers::MakeNonZeroScale(inScale);

	return scale.GetSign() * ScaleHelpers::MakeUniformScale(scale.Abs());
}

void CapsuleShape::sRegister()
{
	ShapeFunctions &f = ShapeFunctions::sGet(EShapeSubType::Capsule);
	f.mConstruct = []() -> Shape * { return new CapsuleShape; };
	f.mColor = Color::sGreen;
}

JPH_NAMESPACE_END
