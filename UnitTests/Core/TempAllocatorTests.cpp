// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-License-Identifier: MIT

#include "UnitTestFramework.h"
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Math/Bivec.h>
#include <Jolt/Math/Rotor.h>
#include <array>
#include <limits>

TEST_SUITE("TempAllocatorTests")
{
	static void sCheckAlignment(const void *inAddress)
	{
		CHECK(inAddress != nullptr);
		CHECK(IsAligned(inAddress, JPH_RVECTOR_ALIGNMENT));
		CHECK(IsAligned(inAddress, alignof(RVec4)));
		CHECK(IsAligned(inAddress, alignof(Vec8)));
		CHECK(IsAligned(inAddress, alignof(Bivec)));
		CHECK(IsAligned(inAddress, alignof(Rotor)));
	}

	static void sCheckAllocations(TempAllocator &inAllocator)
	{
		// In particular, a 16-byte chunk must not misalign the next bivector.
		constexpr uint sizes[] = { 1, 16, 17, 31, 32, 33, 48, 63, 65, 127, 257 };
		std::array<void *, std::size(sizes)> addresses;
		CHECK(inAllocator.Allocate(0) == nullptr);
		inAllocator.Free(nullptr, 0);
		for (size_t i = 0; i < std::size(sizes); ++i)
		{
			CAPTURE(i);
			addresses[i] = inAllocator.Allocate(sizes[i]);
			sCheckAlignment(addresses[i]);
			memset(addresses[i], int(i + 1), sizes[i]);
		}
		for (size_t i = std::size(sizes); i-- > 0;)
		{
			CAPTURE(i);
			const uint8 *bytes = static_cast<const uint8 *>(addresses[i]);
			for (uint j = 0; j < sizes[i]; ++j)
				CHECK(bytes[j] == uint8(i + 1));
			inAllocator.Free(addresses[i], sizes[i]);
		}
	}

	TEST_CASE("FixedBufferAlignsEveryVariedSizeAllocationAndRewinds")
	{
		TempAllocatorImpl allocator(4096);
		sCheckAllocations(allocator);
		CHECK(allocator.IsEmpty());
		CHECK(allocator.GetUsage() == 0);
		// Repeat to verify Free reverses exactly the same padding as Allocate.
		sCheckAllocations(allocator);
		CHECK(allocator.IsEmpty());
	}

	TEST_CASE("MallocAlignsEveryVariedSizeAllocation")
	{
		TempAllocatorMalloc allocator;
		sCheckAllocations(allocator);
	}

	TEST_CASE("FallbackAlignsBothFixedAndOverflowAllocationsAndRewinds")
	{
		// The first small allocations fit; the remaining allocations fall back.
		TempAllocatorImplWithMallocFallback allocator(96);
		sCheckAllocations(allocator);
		sCheckAllocations(allocator);
	}

	TEST_CASE("CapacityChecksDoNotWrapNearUintMax")
	{
		TempAllocatorImpl allocator(128);
		const uint max_size = std::numeric_limits<uint>::max();
		for (uint delta = 0; delta < 64; ++delta)
			CHECK_FALSE(allocator.CanAllocate(max_size - delta));
		CHECK(allocator.GetUsage() == 0);
		void *block = allocator.Allocate(1);
		for (uint delta = 0; delta < 64; ++delta)
			CHECK_FALSE(allocator.CanAllocate(max_size - delta));
		allocator.Free(block, 1);
		CHECK(allocator.IsEmpty());
	}

	TEST_CASE("ExactCapacityZeroSizeAndLIFOReuse")
	{
		// An unaligned capacity must not make the last padded allocation fit.
		constexpr uint alignment = uint(std::max({ size_t(JPH_RVECTOR_ALIGNMENT), alignof(RVec4), alignof(Vec8), alignof(Bivec), alignof(Rotor) }));
		TempAllocatorImpl allocator(3 * alignment - 1);
		void *first = allocator.Allocate(1);
		void *second = allocator.Allocate(alignment);
		CHECK(allocator.GetUsage() == 2 * alignment);
		CHECK(allocator.CanAllocate(0));
		CHECK_FALSE(allocator.CanAllocate(1));
		CHECK(allocator.Allocate(0) == nullptr);
		allocator.Free(nullptr, 0);
		CHECK(allocator.GetUsage() == 2 * alignment);
		allocator.Free(second, alignment);
		CHECK(allocator.Allocate(alignment) == second);
		allocator.Free(second, alignment);
		allocator.Free(first, 1);
		CHECK(allocator.IsEmpty());
	}
}
