// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT
#include "UnitTestFramework.h"
#include <Jolt/Core/StreamWrapper.h>
#include <sstream>

TEST_SUITE("DoublePositionStreamTests")
{
	TEST_CASE("ScalarWireFormatContainsFourDoubles")
	{
		const DVec4 input(0x1.0000000000001p40, -0x1.0000000000003p39, 7.25, -0x1.0000000000005p41);
		std::stringstream stream;
		StreamOutWrapper out(stream);
		out.Write(input);
		CHECK(stream.str().size() == 4 * sizeof(double));
		// Inspect the wire independently of StreamIn; two broken endpoints must not agree.
		double wire[4] = {};
		const std::string bytes = stream.str();
		memcpy(wire, bytes.data(), min(bytes.size(), sizeof(wire)));
		CHECK(wire[0] == input.GetX());
		CHECK(wire[1] == input.GetY());
		CHECK(wire[2] == input.GetZ());
		CHECK(wire[3] == input.GetW());
		DVec4 result = DVec4::sZero();
		StreamInWrapper in(stream);
		in.Read(result);
		CHECK_FALSE(in.IsFailed());
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(result[axis] == input[axis]);
	}

	TEST_CASE("ReadFourDoubleWireFormatAndFollowingSentinel")
	{
		const double wire[4] = { 1.25, -2.5, 3.75, -9.125 };
		const uint32 sentinel = 0x1234abcd;
		std::stringstream stream;
		stream.write(reinterpret_cast<const char *>(wire), sizeof(wire));
		stream.write(reinterpret_cast<const char *>(&sentinel), sizeof(sentinel));
		StreamInWrapper in(stream);
		DVec4 result = DVec4::sZero();
		in.Read(result);
		uint32 next = 0;
		in.Read(next);
		CHECK_FALSE(in.IsFailed());
		for (uint axis = 0; axis < 4; ++axis)
			CHECK(result[axis] == wire[axis]);
		CHECK(next == sentinel);
	}

	TEST_CASE("ArrayWireFormatContainsAllFourCoordinates")
	{
		Array<DVec4> input = { DVec4(1, 2, 3, 17), DVec4(4, 5, 6, -23) };
		std::stringstream stream;
		StreamOutWrapper out(stream);
		out.Write(input);
		CHECK(stream.str().size() == sizeof(uint32) + 8 * sizeof(double));
		StreamInWrapper in(stream);
		Array<DVec4> result;
		in.Read(result);
		CHECK_FALSE(in.IsFailed());
		CHECK(result.size() == input.size());
		if (in.IsFailed() || result.size() != input.size())
			return;
		for (size_t i = 0; i < input.size(); ++i)
			for (uint axis = 0; axis < 4; ++axis)
				CHECK(result[i][axis] == input[i][axis]);
	}

	TEST_CASE("TruncatedFourthDoubleIsFailure")
	{
		const double wire[3] = { 1, 2, 3 };
		std::stringstream stream;
		stream.write(reinterpret_cast<const char *>(wire), sizeof(wire));
		StreamInWrapper in(stream);
		DVec4 result = DVec4::sZero();
		in.Read(result);
		CHECK(in.IsFailed());
		CHECK(in.IsEOF());
	}
}
