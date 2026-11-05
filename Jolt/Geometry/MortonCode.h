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
		v = (v | (v << 12)) & 0x000F000Fu;
		v = (v | (v <<  6)) & 0x03030303u;
		v = (v | (v <<  3)) & 0x11111111u;
		return v;
	}

	/// Calculate a 32-bit 4D Morton key: 8 bits/axis, nearest quantization to [0, 255].
	/// Bit 4*b+3/2/1/0 stores X/Y/Z/W bit b respectively. Distinct quantized
	/// coordinates are injective; points in the same quantization cell intentionally share a key.
	/// Bounds must contain the point; zero-extent axes quantize to zero.
	static uint32 sGetMortonCode(Vec4Arg inVector, const AABox &inVectorBounds)
	{
		// Convert to 8 bit fixed point
		Vec4 size = inVectorBounds.GetSize();
		Vec4 safe_size = Vec4::sSelect(size, Vec4::sOne(), Vec4::sEquals(size, Vec4::sZero()));
		Vec4 scaled = Vec4::sClamp((inVector - inVectorBounds.mMin) / safe_size, Vec4::sZero(), Vec4::sOne());
		uint x = sExpandBits(scaled.GetX());
		uint y = sExpandBits(scaled.GetY());
		uint z = sExpandBits(scaled.GetZ());
		uint w = sExpandBits(scaled.GetW());
		return (x << 3) | (y << 2) | (z << 1) | w;
	}
};

JPH_NAMESPACE_END
