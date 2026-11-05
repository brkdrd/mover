// Jolt Physics Library (https://github.com/jrouwe/JoltPhysics)
// SPDX-FileCopyrightText: 2021 Jorrit Rouwe
// SPDX-License-Identifier: MIT

#pragma once

#include <Jolt/Core/NonCopyable.h>
#include <Jolt/Math/Bivec.h>
#include <Jolt/Math/Real.h>
#include <Jolt/Math/Rotor.h>

JPH_NAMESPACE_BEGIN

/// Allocator for temporary allocations.
/// This allocator works as a stack: The blocks must always be freed in the reverse order as they are allocated.
/// Note that allocations and frees can take place from different threads, but the order is guaranteed though
/// job dependencies, so it is not needed to use any form of locking.
class JPH_EXPORT TempAllocator : public NonCopyable
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Required alignment for all temporary allocations, including custom TempAllocator implementations.
	/// 4D angular quantities can require wider SIMD alignment than single-precision spatial vectors.
	static constexpr size_t			cAlignment = std::max({ size_t(JPH_RVECTOR_ALIGNMENT), alignof(RVec4), alignof(Vec8), alignof(Bivec), alignof(Rotor) });

	/// If the platform default allocator cannot satisfy the temporary allocation contract
	static constexpr bool			needs_aligned_allocate = cAlignment > JPH_DEFAULT_ALLOCATE_ALIGNMENT;

	/// Destructor
	virtual							~TempAllocator() = default;

	/// Allocates inSize bytes of memory, returned memory address must be cAlignment byte aligned (also for custom allocators).
	/// A zero-size allocation returns nullptr.
	virtual void *					Allocate(uint inSize) = 0;

	/// Frees inSize bytes of memory located at inAddress
	virtual void					Free(void *inAddress, uint inSize) = 0;
};

/// Default implementation of the temp allocator that allocates a large block through malloc upfront
class JPH_EXPORT TempAllocatorImpl final : public TempAllocator
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructs the allocator with a maximum allocatable size of inSize
	explicit						TempAllocatorImpl(size_t inSize) : mSize(inSize)
	{
		if constexpr (needs_aligned_allocate)
			mBase = static_cast<uint8 *>(AlignedAllocate(inSize, cAlignment));
		else
			mBase = static_cast<uint8 *>(JPH::Allocate(inSize));
	}

	/// Destructor, frees the block
	virtual							~TempAllocatorImpl() override
	{
		JPH_ASSERT(mTop == 0);
		if constexpr (needs_aligned_allocate)
			AlignedFree(mBase);
		else
			JPH::Free(mBase);
	}

	// See: TempAllocator
	virtual void *					Allocate(uint inSize) override
	{
		if (inSize == 0)
		{
			return nullptr;
		}
		else
		{
			if (!CanAllocate(inSize))
			{
				Trace("TempAllocator: Out of memory trying to allocate %u bytes", inSize);
				std::abort();
			}
			void *address = mBase + mTop;
			mTop += AlignUp(size_t(inSize), cAlignment);
			return address;
		}
	}

	// See: TempAllocator
	virtual void					Free(void *inAddress, uint inSize) override
	{
		if (inAddress == nullptr)
		{
			JPH_ASSERT(inSize == 0);
		}
		else
		{
			const size_t aligned_size = AlignUp(size_t(inSize), cAlignment);
			if (aligned_size < inSize || aligned_size > mTop || mBase + (mTop - aligned_size) != inAddress)
			{
				Trace("TempAllocator: Freeing in the wrong order");
				std::abort();
			}
			mTop -= aligned_size;
		}
	}

	/// Check if no allocations have been made
	bool							IsEmpty() const
	{
		return mTop == 0;
	}

	/// Get the total size of the fixed buffer
	size_t							GetSize() const
	{
		return mSize;
	}

	/// Get current usage in bytes of the buffer
	size_t							GetUsage() const
	{
		return mTop;
	}

	/// Check if an allocation of inSize can be made in this fixed buffer allocator
	bool							CanAllocate(uint inSize) const
	{
		// Widen before rounding (uint can wrap near UINT_MAX). The first check also
		// rejects rounding overflow on platforms where size_t is no wider than uint.
		const size_t aligned_size = AlignUp(size_t(inSize), cAlignment);
		return aligned_size >= inSize && aligned_size <= mSize - mTop;
	}

	/// Check if memory block at inAddress is owned by this allocator
	bool							OwnsMemory(const void *inAddress) const
	{
		return inAddress >= mBase && inAddress < mBase + mSize;
	}

private:
	uint8 *							mBase;							///< Base address of the memory block
	size_t							mSize;							///< Size of the memory block
	size_t							mTop = 0;						///< End of currently allocated area
};

/// Implementation of the TempAllocator that just falls back to malloc/free
/// Note: This can be quite slow when running in the debugger as large memory blocks need to be initialized with 0xcd
class JPH_EXPORT TempAllocatorMalloc final : public TempAllocator
{
public:
	JPH_OVERRIDE_NEW_DELETE

	// See: TempAllocator
	virtual void *					Allocate(uint inSize) override
	{
		if (inSize > 0)
		{
			if constexpr (needs_aligned_allocate)
				return AlignedAllocate(inSize, cAlignment);
			else
				return JPH::Allocate(inSize);
		}
		else
			return nullptr;
	}

	// See: TempAllocator
	virtual void					Free(void *inAddress, [[maybe_unused]] uint inSize) override
	{
		if (inAddress != nullptr)
		{
			if constexpr (needs_aligned_allocate)
				AlignedFree(inAddress);
			else
				JPH::Free(inAddress);
		}
	}
};

/// Implementation of the TempAllocator that tries to allocate from a large preallocated block, but falls back to malloc when it is exhausted
class JPH_EXPORT TempAllocatorImplWithMallocFallback final : public TempAllocator
{
public:
	JPH_OVERRIDE_NEW_DELETE

	/// Constructs the allocator with an initial fixed block if inSize
	explicit						TempAllocatorImplWithMallocFallback(uint inSize) :
		mAllocator(inSize)
	{
	}

	// See: TempAllocator
	virtual void *					Allocate(uint inSize) override
	{
		if (mAllocator.CanAllocate(inSize))
			return mAllocator.Allocate(inSize);
		else
			return mFallbackAllocator.Allocate(inSize);
	}

	// See: TempAllocator
	virtual void					Free(void *inAddress, uint inSize) override
	{
		if (inAddress == nullptr)
		{
			JPH_ASSERT(inSize == 0);
		}
		else
		{
			if (mAllocator.OwnsMemory(inAddress))
				mAllocator.Free(inAddress, inSize);
			else
				mFallbackAllocator.Free(inAddress, inSize);
		}
	}

private:
	TempAllocatorImpl				mAllocator;
	TempAllocatorMalloc				mFallbackAllocator;
};

JPH_NAMESPACE_END
