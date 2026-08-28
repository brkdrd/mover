// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2023 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

/// Enum used in BodyCreationSettings and MotionProperties to indicate which degrees of freedom a body has.
/// In 4D: 4 translation DOFs + 6 rotation DOFs (one per coordinate plane) = 10 bits.
enum class EAllowedDOFs : uint16
{
	None				= 0,											///< No degrees of freedom are allowed. Note that this is not valid and will crash. Use a static body instead.
	TranslationX		= 1 << 0,										///< Body can move in world space X axis
	TranslationY		= 1 << 1,										///< Body can move in world space Y axis
	TranslationZ		= 1 << 2,										///< Body can move in world space Z axis
	TranslationW		= 1 << 3,										///< Body can move in world space W axis
	RotationXY			= 1 << 4,										///< Body can rotate in the XY plane (e12 bivector component)
	RotationXZ			= 1 << 5,										///< Body can rotate in the XZ plane (e13 bivector component)
	RotationXW			= 1 << 6,										///< Body can rotate in the XW plane (e14 bivector component)
	RotationYZ			= 1 << 7,										///< Body can rotate in the YZ plane (e23 bivector component)
	RotationYW			= 1 << 8,										///< Body can rotate in the YW plane (e24 bivector component)
	RotationZW			= 1 << 9,										///< Body can rotate in the ZW plane (e34 bivector component)
	All					= 0b1111111111,									///< All degrees of freedom are allowed
	AllTranslation		= TranslationX | TranslationY | TranslationZ | TranslationW,
	AllRotation			= RotationXY | RotationXZ | RotationXW | RotationYZ | RotationYW | RotationZW,
};

/// Bitwise OR operator for EAllowedDOFs
constexpr EAllowedDOFs operator | (EAllowedDOFs inLHS, EAllowedDOFs inRHS)
{
	return EAllowedDOFs(uint16(inLHS) | uint16(inRHS));
}

/// Bitwise AND operator for EAllowedDOFs
constexpr EAllowedDOFs operator & (EAllowedDOFs inLHS, EAllowedDOFs inRHS)
{
	return EAllowedDOFs(uint16(inLHS) & uint16(inRHS));
}

/// Bitwise XOR operator for EAllowedDOFs
constexpr EAllowedDOFs operator ^ (EAllowedDOFs inLHS, EAllowedDOFs inRHS)
{
	return EAllowedDOFs(uint16(inLHS) ^ uint16(inRHS));
}

/// Bitwise NOT operator for EAllowedDOFs
constexpr EAllowedDOFs operator ~ (EAllowedDOFs inAllowedDOFs)
{
	return EAllowedDOFs(~uint16(inAllowedDOFs));
}

/// Bitwise OR assignment operator for EAllowedDOFs
constexpr EAllowedDOFs & operator |= (EAllowedDOFs &ioLHS, EAllowedDOFs inRHS)
{
	ioLHS = ioLHS | inRHS;
	return ioLHS;
}

/// Bitwise AND assignment operator for EAllowedDOFs
constexpr EAllowedDOFs & operator &= (EAllowedDOFs &ioLHS, EAllowedDOFs inRHS)
{
	ioLHS = ioLHS & inRHS;
	return ioLHS;
}

/// Bitwise XOR assignment operator for EAllowedDOFs
constexpr EAllowedDOFs & operator ^= (EAllowedDOFs &ioLHS, EAllowedDOFs inRHS)
{
	ioLHS = ioLHS ^ inRHS;
	return ioLHS;
}

JPH_NAMESPACE_END
