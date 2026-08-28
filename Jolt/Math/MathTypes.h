// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

JPH_NAMESPACE_BEGIN

class Vec3;
class DVec4;
class Lane4;
class Vec4;
class Vec8;
class Bivec;
class Rotor;
class UVec4;
class BVec16;
class Quat;
class Mat44;
class RMat44;

// Types to use for passing arguments to functions
using Vec3Arg = const Vec3;
#ifdef JPH_USE_AVX
	using DVec4Arg = const DVec4;
#else
	using DVec4Arg = const DVec4 &;
#endif
using Lane4Arg = const Lane4;
using Vec4Arg = const Vec4;
using Vec8Arg = const Vec8;
using UVec4Arg = const UVec4;
using BVec16Arg = const BVec16;
using QuatArg = const Quat;
using BivecArg = const Bivec &;
using RotorArg = const Rotor &;
using Mat44Arg = const Mat44 &;
using RMat44Arg = const RMat44 &;

JPH_NAMESPACE_END
