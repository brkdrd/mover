// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Geometry/AABox.h>

JPH_NAMESPACE_BEGIN

class MortonCode
{
public:
	/// Converts a floating point value in the range [0, 1] to an 8-bit fixed point integer.
	/// Then expands the 8-bit integer into 32 bits by inserting 3 zeros after each bit (for 4D interleaving).
	static uint32 sExpandBits(float inV)
	{
		JPH_ASSERT(inV >= 0.0f && inV <= 1.0f);
		uint32 v = uint32(inV * 255.0f + 0.5f);
		JPH_ASSERT(v < 256);
		// Spread 8 bits across 32 bits with 3-bit gaps
		v = (v | (v << 16)) & 0x030000FFu;
		v = (v | (v <<  8)) & 0x0300F00Fu;
		v = (v | (v <<  4)) & 0x030C30C3u;
		v = (v | (v <<  2)) & 0x09249249u;
		return v;
	}

	/// Calculate the 4D morton code for inVector, given that all vectors lie in inVectorBounds
	static uint32 sGetMortonCode(Vec4Arg inVector, const AABox &inVectorBounds)
	{
		// Convert to 8 bit fixed point
		Vec4 scaled = (inVector - inVectorBounds.mMin) / inVectorBounds.GetSize();
		uint x = sExpandBits(scaled.GetX());
		uint y = sExpandBits(scaled.GetY());
		uint z = sExpandBits(scaled.GetZ());
		uint w = sExpandBits(scaled.GetW());
		return (x << 3) + (y << 2) + (z << 1) + w;
	}
};

JPH_NAMESPACE_END
