// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include <Jolt/Jolt.h>

#include <Jolt/Physics/Collision/Shape/SphereShape.h>
#include <Jolt/Physics/Collision/Shape/ScaleHelpers.h>
#include <Jolt/Physics/Collision/RayCast.h>
#include <Jolt/Physics/Collision/CastResult.h>
#include <Jolt/Physics/Collision/CollidePointResult.h>
#include <Jolt/Physics/Collision/TransformedShape.h>
#include <Jolt/Physics/Collision/CollideSoftBodyVertexIterator.h>
#include <Jolt/Geometry/RaySphere.h>
#include <Jolt/Geometry/Plane.h>
#include <Jolt/Core/StreamIn.h>
#include <Jolt/Core/StreamOut.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>
#ifdef JPH_DEBUG_RENDERER
	#include <Jolt/Renderer/DebugRenderer.h>
#endif // JPH_DEBUG_RENDERER

JPH_NAMESPACE_BEGIN

JPH_IMPLEMENT_SERIALIZABLE_VIRTUAL(SphereShapeSettings)
{
	JPH_ADD_BASE_CLASS(SphereShapeSettings, ConvexShapeSettings)

	JPH_ADD_ATTRIBUTE(SphereShapeSettings, mRadius)
}

ShapeSettings::ShapeResult SphereShapeSettings::Create() const
{
	if (mCachedResult.IsEmpty())
		Ref<Shape> shape = new SphereShape(*this, mCachedResult);
	return mCachedResult;
}

SphereShape::SphereShape(const SphereShapeSettings &inSettings, ShapeResult &outResult) :
	ConvexShape(EShapeSubType::Sphere, inSettings, outResult),
	mRadius(inSettings.mRadius)
{
	if (inSettings.mRadius <= 0.0f)
	{
		outResult.SetError("Invalid radius");
		return;
	}

	outResult.Set(this);
}

float SphereShape::GetScaledRadius(Vec4Arg inScale) const
{
	JPH_ASSERT(IsValidScale(inScale));

	Vec4 abs_scale = inScale.Abs();
	return abs_scale.GetX() * mRadius;
}

AABox SphereShape::GetLocalBounds() const
{
	Vec4 half_extent = Vec4::sReplicate(mRadius);
	return AABox(-half_extent, half_extent);
}

AABox SphereShape::GetWorldSpaceBounds(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale) const
{
	float scaled_radius = GetScaledRadius(inScale);
	Vec4 half_extent = Vec4::sReplicate(scaled_radius);
	AABox bounds(-half_extent, half_extent);
	bounds.Translate(inCenterOfMassTransform.GetTranslation());
	return bounds;
}

class SphereShape::SphereNoConvex final : public Support
{
public:
	explicit		SphereNoConvex(float inRadius) :
		mRadius(inRadius)
	{
		static_assert(sizeof(SphereNoConvex) <= sizeof(SupportBuffer), "Buffer size too small");
		JPH_ASSERT(IsAligned(this, alignof(SphereNoConvex)));
	}

	virtual Vec4	GetSupport(Vec4Arg inDirection) const override
	{
		return Vec4::sZero();
	}

	virtual float	GetConvexRadius() const override
	{
		return mRadius;
	}

private:
	float			mRadius;
};

class SphereShape::SphereWithConvex final : public Support
{
public:
	explicit		SphereWithConvex(float inRadius) :
		mRadius(inRadius)
	{
		static_assert(sizeof(SphereWithConvex) <= sizeof(SupportBuffer), "Buffer size too small");
		JPH_ASSERT(IsAligned(this, alignof(SphereWithConvex)));
	}

	virtual Vec4	GetSupport(Vec4Arg inDirection) const override
	{
		float len = inDirection.Length();
		return len > 0.0f? (mRadius / len) * inDirection : Vec4::sZero();
	}

	virtual float	GetConvexRadius() const override
	{
		return 0.0f;
	}

private:
	float			mRadius;
};

const ConvexShape::Support *SphereShape::GetSupportFunction(ESupportMode inMode, SupportBuffer &inBuffer, Vec4Arg inScale) const
{
	float scaled_radius = GetScaledRadius(inScale);

	switch (inMode)
	{
	case ESupportMode::IncludeConvexRadius:
		return new (&inBuffer) SphereWithConvex(scaled_radius);

	case ESupportMode::ExcludeConvexRadius:
	case ESupportMode::Default:
		return new (&inBuffer) SphereNoConvex(scaled_radius);
	}

	JPH_ASSERT(false);
	return nullptr;
}

MassProperties SphereShape::GetMassProperties() const
{
	MassProperties p;

	// Calculate mass: 4-ball hypervolume = pi^2/2 * r^4, multiplied by density
	float r2 = mRadius * mRadius;
	float r4 = r2 * r2;
	p.mMass = (0.5f * JPH_PI * JPH_PI) * r4 * GetDensity();

	// MassProperties stores the second moment M_ij = integral rho*x_i*x_j dV,
	// not planar inertia. A solid 4-ball has M_ii = mass*r^2/6; the
	// bivector moment for plane ij is M_ii + M_jj = mass*r^2/3.
	float inertia = (1.0f / 6.0f) * p.mMass * r2;
	Mat44 m = Mat44::sZero();
	m(0, 0) = inertia;
	m(1, 1) = inertia;
	m(2, 2) = inertia;
	m(3, 3) = inertia;
	p.mInertia = m;

	return p;
}

Vec4 SphereShape::GetSurfaceNormal(const SubShapeID &inSubShapeID, Vec4Arg inLocalSurfacePosition) const
{
	JPH_ASSERT(inSubShapeID.IsEmpty(), "Invalid subshape ID");

	float len = inLocalSurfacePosition.Length();
	return len != 0.0f? inLocalSurfacePosition / len : Vec4::sAxisY();
}

void SphereShape::GetSubmergedVolume(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const Plane &inSurface, float &outTotalVolume, float &outSubmergedVolume, Vec4 &outCenterOfBuoyancy, RVec4Arg inBaseOffset) const
{
	float scaled_radius = GetScaledRadius(inScale);
	float r2 = scaled_radius * scaled_radius;
	// 4-ball hypervolume = pi^2/2 * r^4
	outTotalVolume = (0.5f * JPH_PI * JPH_PI) * r2 * r2;

	float distance_to_surface = inSurface.SignedDistance(Vec4(inCenterOfMassTransform.GetTranslation()));
	if (distance_to_surface >= scaled_radius)
	{
		// Above surface
		outSubmergedVolume = 0.0f;
		outCenterOfBuoyancy = Vec4::sZero();
	}
	else if (distance_to_surface <= -scaled_radius)
	{
		// Under surface
		outSubmergedVolume = outTotalVolume;
		outCenterOfBuoyancy = Vec4(inCenterOfMassTransform.GetTranslation());
	}
	else
	{
		// Intersecting surface.
		// TODO(4D): Exact submerged-hypervolume of a 4-ball cut by a hyperplane involves a 4D
		// spherical cap formula (incomplete beta function). Placeholder: linear interpolation
		// between 0 and full volume by signed distance across the diameter. This is coarse but
		// keeps buoyancy approximately continuous until a closed-form is added.
		float t = 0.5f - distance_to_surface / (2.0f * scaled_radius);
		outSubmergedVolume = t * outTotalVolume;

		// Center-of-buoyancy placeholder: the sphere center offset into the submerged half by a fraction
		// of the radius proportional to the cap thickness.
		float z = 0.5f * (scaled_radius - distance_to_surface);
		outCenterOfBuoyancy = Vec4(inCenterOfMassTransform.GetTranslation()) - z * inSurface.GetNormal();

	#ifdef JPH_DEBUG_RENDERER
		(void)inBaseOffset;
	#endif // JPH_DEBUG_RENDERER
	}

#ifdef JPH_DEBUG_RENDERER
	// TODO(4D): Submerged-volume visualization uses 3D DebugRenderer APIs (DrawPie/DrawWireSphere) that
	// operate on RVec3. Deferred until the renderer is migrated to 4D.
	(void)inBaseOffset;
#endif // JPH_DEBUG_RENDERER
}

#ifdef JPH_DEBUG_RENDERER
void SphereShape::Draw(DebugRenderer *inRenderer, RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, ColorArg inColor, bool inUseMaterialColors, bool inDrawWireframe) const
{
	// TODO(4D): DebugRenderer's DrawUnitSphere expects a 3D transform. Deferred until renderer 4D migration.
	(void)inRenderer; (void)inCenterOfMassTransform; (void)inScale; (void)inColor; (void)inUseMaterialColors; (void)inDrawWireframe;
}
#endif // JPH_DEBUG_RENDERER

bool SphereShape::CastRay(const RayCast &inRay, const SubShapeIDCreator &inSubShapeIDCreator, RayCastResult &ioHit) const
{
	float fraction = RaySphere(inRay.mOrigin, inRay.mDirection, Vec4::sZero(), mRadius);
	if (fraction < ioHit.mFraction)
	{
		ioHit.mFraction = fraction;
		ioHit.mSubShapeID2 = inSubShapeIDCreator.GetID();
		return true;
	}
	return false;
}

void SphereShape::CastRay(const RayCast &inRay, const RayCastSettings &inRayCastSettings, const SubShapeIDCreator &inSubShapeIDCreator, CastRayCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	float min_fraction, max_fraction;
	int num_results = RaySphere(inRay.mOrigin, inRay.mDirection, Vec4::sZero(), mRadius, min_fraction, max_fraction);
	if (num_results > 0 // Ray should intersect
		&& max_fraction >= 0.0f // End of ray should be inside sphere
		&& min_fraction < ioCollector.GetEarlyOutFraction()) // Start of ray should be before early out fraction
	{
		// Better hit than the current hit
		RayCastResult hit;
		hit.mBodyID = TransformedShape::sGetBodyID(ioCollector.GetContext());
		hit.mSubShapeID2 = inSubShapeIDCreator.GetID();

		// Check front side hit
		if (inRayCastSettings.mTreatConvexAsSolid || min_fraction > 0.0f)
		{
			hit.mFraction = max(0.0f, min_fraction);
			ioCollector.AddHit(hit);
		}

		// Check back side hit
		if (inRayCastSettings.mBackFaceModeConvex == EBackFaceMode::CollideWithBackFaces
			&& num_results > 1 // Ray should have 2 intersections
			&& max_fraction < ioCollector.GetEarlyOutFraction()) // End of ray should be before early out fraction
		{
			hit.mFraction = max_fraction;
			ioCollector.AddHit(hit);
		}
	}
}

void SphereShape::CollidePoint(Vec4Arg inPoint, const SubShapeIDCreator &inSubShapeIDCreator, CollidePointCollector &ioCollector, const ShapeFilter &inShapeFilter) const
{
	// Test shape filter
	if (!inShapeFilter.ShouldCollide(this, inSubShapeIDCreator.GetID()))
		return;

	if (inPoint.LengthSq() <= Square(mRadius))
		ioCollector.AddHit({ TransformedShape::sGetBodyID(ioCollector.GetContext()), inSubShapeIDCreator.GetID() });
}

void SphereShape::CollideSoftBodyVertices(RMat44Arg inCenterOfMassTransform, Vec4Arg inScale, const CollideSoftBodyVertexIterator &inVertices, uint inNumVertices, int inCollidingShapeIndex) const
{
	// TODO(4D): CollideSoftBodyVertexIterator still exposes Vec3 positions. Once the soft-body
	// subsystem is migrated to 4D, rewrite this loop to operate on Vec4 directly.
	(void)inCenterOfMassTransform; (void)inScale; (void)inVertices; (void)inNumVertices; (void)inCollidingShapeIndex;
}

void SphereShape::GetTrianglesStart(GetTrianglesContext &ioContext, const AABox &inBox, Vec4Arg inPositionCOM, RotorArg inRotation, Vec4Arg inScale) const
{
	// TODO(4D): GetTrianglesContextVertexList still uses the 3D vertex pipeline. Replaced with a no-op placeholder
	// until the 4D tessellation (tetrahedra-based boundary) is defined.
	(void)ioContext; (void)inBox; (void)inPositionCOM; (void)inRotation; (void)inScale;
}

int SphereShape::GetTrianglesNext(GetTrianglesContext &ioContext, int inMaxTrianglesRequested, Float4 *outTriangleVertices, const PhysicsMaterial **outMaterials) const
{
	// TODO(4D): See GetTrianglesStart. Returns 0 to indicate no tessellation data is currently produced.
	(void)ioContext; (void)inMaxTrianglesRequested; (void)outTriangleVertices; (void)outMaterials;
	return 0;
}

void SphereShape::SaveBinaryState(StreamOut &inStream) const
{
	ConvexShape::SaveBinaryState(inStream);

	inStream.Write(mRadius);
}

void SphereShape::RestoreBinaryState(StreamIn &inStream)
{
	ConvexShape::RestoreBinaryState(inStream);

	inStream.Read(mRadius);
}

bool SphereShape::IsValidScale(Vec4Arg inScale) const
{
	return ConvexShape::IsValidScale(inScale) && ScaleHelpers::IsUniformScale(inScale.Abs());
}

Vec4 SphereShape::MakeScaleValid(Vec4Arg inScale) const
{
	Vec4 scale = ScaleHelpers::MakeNonZeroScale(inScale);

	return scale.GetSign() * ScaleHelpers::MakeUniformScale(scale.Abs());
}

void SphereShape::sRegister()
{
	ShapeFunctions &f = ShapeFunctions::sGet(EShapeSubType::Sphere);
	f.mConstruct = []() -> Shape * { return new SphereShape; };
	f.mColor = Color::sGreen;
}

JPH_NAMESPACE_END
