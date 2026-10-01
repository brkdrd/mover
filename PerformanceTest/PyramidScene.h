// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2023 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

// Jolt includes
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>

// Local includes
#include "PerformanceTestScene.h"
#include "Layers.h"

// A scene that creates a pyramid of boxes to create a very large island
// TODO(4D): a 4D pyramid has three horizontal axes, so it is a 4-deep loop and its body count grows
// as the sum of cubes rather than squares. Upstream's height of 15 would give
// (15*16/2)^2 = 14400 bodies where 3D gets 1240; the height is lowered to 8 ((8*9/2)^2 = 1296) so
// the scene stays the same size as the 3D one it is meant to be compared against.
class PyramidScene : public PerformanceTestScene
{
public:
	virtual const char *	GetName() const override
	{
		return "Pyramid";
	}

	virtual void			StartTest(PhysicsSystem &inPhysicsSystem, EMotionQuality inMotionQuality) override
	{
		BodyInterface &bi = inPhysicsSystem.GetBodyInterface();

		// Floor. Wide in W as well - in 4D the floor's boundary is a 3D volume.
		bi.CreateAndAddBody(BodyCreationSettings(new BoxShape(Vec4(50.0f, 1.0f, 50.0f, 50.0f), 0.0f), RVec4(Vec4(0.0f, -1.0f, 0.0f, 0.0f)), Rotor::sIdentity(), EMotionType::Static, Layers::NON_MOVING), EActivation::DontActivate);

		const float cBoxSize = 2.0f;
		const float cBoxSeparation = 0.5f;
		const float cHalfBoxSize = 0.5f * cBoxSize;
		const int cPyramidHeight = 8;

		RefConst<Shape> box_shape = new BoxShape(Vec4::sReplicate(cHalfBoxSize), 0.0f); // No convex radius to force more collisions

		// Pyramid
		for (int i = 0; i < cPyramidHeight; ++i)
			for (int j = i / 2; j < cPyramidHeight - (i + 1) / 2; ++j)
				for (int k = i / 2; k < cPyramidHeight - (i + 1) / 2; ++k)
					for (int l = i / 2; l < cPyramidHeight - (i + 1) / 2; ++l)
					{
						float offset = (i & 1)? cHalfBoxSize : 0.0f;
						RVec4 position(Real(-cPyramidHeight + cBoxSize * j + offset), Real(1.0f + (cBoxSize + cBoxSeparation) * i), Real(-cPyramidHeight + cBoxSize * k + offset), Real(-cPyramidHeight + cBoxSize * l + offset));
						BodyCreationSettings settings(box_shape, position, Rotor::sIdentity(), EMotionType::Dynamic, Layers::MOVING);
						settings.mAllowSleeping = false; // No sleeping to force the large island to stay awake
						bi.CreateAndAddBody(settings, EActivation::Activate);
					}
	}
};
