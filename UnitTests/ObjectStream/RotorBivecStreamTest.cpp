// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Core/Factory.h>
#include <Jolt/Math/Rotor.h>
#include <Jolt/Math/Bivec.h>
#include <Jolt/ObjectStream/ObjectStreamOut.h>
#include <Jolt/ObjectStream/ObjectStreamIn.h>
#include <Jolt/ObjectStream/SerializableObject.h>
#include <Jolt/ObjectStream/TypeDeclarations.h>

// Round-trips the 4D rotation primitives (Vec4 / Rotor / Bivec) through the text and binary
// object streams (Phase B step 5). Rotor and Bivec were newly registered as ObjectStream
// primitives; this pins that they serialize losslessly.
JPH_SUPPRESS_WARNINGS_STD_BEGIN
#include <sstream>
JPH_SUPPRESS_WARNINGS_STD_END

namespace RotorBivecStreamTestNS {

class RotorSer
{
	JPH_DECLARE_SERIALIZABLE_NON_VIRTUAL(JPH_NO_EXPORT, RotorSer)
public:
	Vec4	mVec4 = Vec4::sZero();
	DVec4 mDVec4 = DVec4::sZero();
	Array<DVec4> mPositions;
	Rotor	mRotor = Rotor::sIdentity();
	Bivec	mBivec = Bivec::sZero();
};

JPH_IMPLEMENT_SERIALIZABLE_NON_VIRTUAL(RotorSer)
{
	JPH_ADD_ATTRIBUTE(RotorSer, mVec4)
	JPH_ADD_ATTRIBUTE(RotorSer, mDVec4)
	JPH_ADD_ATTRIBUTE(RotorSer, mPositions)
	JPH_ADD_ATTRIBUTE(RotorSer, mRotor)
	JPH_ADD_ATTRIBUTE(RotorSer, mBivec)
}

} // namespace RotorBivecStreamTestNS

using RotorBivecStreamTestNS::RotorSer;

TEST_SUITE("RotorBivecStreamTest")
{
	static void sRoundTrip(ObjectStreamOut::EStreamType inType)
	{
		Factory::sInstance->Register(JPH_RTTI(RotorSer));

		RotorSer in;
		in.mVec4 = Vec4(1.5f, -2.5f, 3.5f, -4.5f);
		in.mDVec4 = DVec4(0x1.0000000000001p40, -0x1.0000000000003p39, 7.25, -0x1.0000000000005p41);
		in.mPositions = { in.mDVec4, DVec4(11, -13, 17, -19) };
		in.mRotor = Rotor::sRotation(Vec4::sAxisX(), Vec4::sAxisY(), 0.7f) * Rotor::sRotation(Vec4::sAxisZ(), Vec4::sAxisW(), -1.3f);
		in.mBivec = Bivec(0.1f, -0.2f, 0.3f, -0.4f, 0.5f, -0.6f);

		stringstream stream;
		CHECK(ObjectStreamOut::sWriteObject(stream, inType, in));
		CHECK(stream.str().substr(0, 8) == (inType == ObjectStreamOut::EStreamType::Text? "TOS 2.00" : "BOS 2.00"));

		RotorSer *out = nullptr;
		CHECK(ObjectStreamIn::sReadObject(stream, out));
		CHECK(out != nullptr);
		if (out != nullptr)
		{
			CHECK(out->mVec4.IsClose(in.mVec4, 1.0e-10f));
			for (uint axis = 0; axis < 4; ++axis)
				CHECK(out->mDVec4[axis] == in.mDVec4[axis]);
			CHECK(out->mPositions.size() == in.mPositions.size());
			if (out->mPositions.size() != in.mPositions.size())
			{
				delete out;
				return;
			}
			for (size_t i = 0; i < in.mPositions.size(); ++i)
				for (uint axis = 0; axis < 4; ++axis)
					CHECK(out->mPositions[i][axis] == in.mPositions[i][axis]);
			CHECK(out->mRotor.GetVec8().IsClose(in.mRotor.GetVec8(), 1.0e-10f));
			CHECK(out->mBivec.GetVec8().IsClose(in.mBivec.GetVec8(), 1.0e-10f));
			delete out;
		}
	}

	TEST_CASE("FourCoordinatePrecisionConversion")
	{
		Factory::sInstance->Register(JPH_RTTI(RotorSer));
		RotorSer input;
		input.mVec4 = Vec4(1, -2, 3, -11);
		input.mDVec4 = DVec4(5, -6, 7, -13);
		for (bool widen : { false, true })
		{
			stringstream stream;
			CHECK(ObjectStreamOut::sWriteObject(stream, ObjectStreamOut::EStreamType::Text, input));
			std::string bytes = stream.str();
			const char *source_type = widen? "dvec4" : "vec4";
			size_t position = bytes.find(source_type);
			CHECK(position != std::string::npos);
			if (position == std::string::npos)
				return;
			// Alter only the scalar attribute's declared source precision, not its four values.
			bytes.replace(position, strlen(source_type), widen? "vec4" : "dvec4");
			stringstream converted(bytes);
			RotorSer *output = nullptr;
			CHECK(ObjectStreamIn::sReadObject(converted, output));
			CHECK(output != nullptr);
			if (output == nullptr)
				return;
			for (uint axis = 0; axis < 4; ++axis)
			{
				CHECK(output->mVec4[axis] == input.mVec4[axis]);
				CHECK(output->mDVec4[axis] == input.mDVec4[axis]);
			}
			delete output;
		}
	}

	TEST_CASE("LegacyThreeDoubleObjectStreamIsRejected")
	{
		Factory::sInstance->Register(JPH_RTTI(RotorSer));
		for (auto type : { ObjectStreamOut::EStreamType::Text, ObjectStreamOut::EStreamType::Binary })
		{
			stringstream stream;
			CHECK(ObjectStreamOut::sWriteObject(stream, type, RotorSer()));
			std::string bytes = stream.str();
			bytes.replace(0, 8, type == ObjectStreamOut::EStreamType::Text? "TOS 1.00" : "BOS 1.00");
			stringstream old_stream(bytes);
			RotorSer *out = nullptr;
			CHECK_FALSE(ObjectStreamIn::sReadObject(old_stream, out));
			CHECK(out == nullptr);
			delete out;
		}
	}

	TEST_CASE("TestRotorBivecStreamText")
	{
		sRoundTrip(ObjectStreamOut::EStreamType::Text);
	}

	TEST_CASE("TestRotorBivecStreamBinary")
	{
		sRoundTrip(ObjectStreamOut::EStreamType::Binary);
	}
}
