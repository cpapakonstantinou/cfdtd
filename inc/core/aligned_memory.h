#ifndef __ALIGNED_MEMORY_H__
#define __ALIGNED_MEMORY_H__

// Copyright 2025  Constantine Papakonstantinou
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files, the "Software", to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

/**
 * \file core/aligned_memory.h header only support for smart memory aligned pointers
 * \author cpapakonstantinou
 * \date 2021
 **/
#include <array>
#include <cstddef>
#include <memory>
#include <stdlib.h>
#include <sys/mman.h>
#include <unordered_map>

/** \brief  Bytes held in aligned arrays, by group: what `fdtd mesh` reports and `memcfg` estimates.
 *
 *  Every array of the solver is allocated by palign() and freed by pdelete(), so the tally is exact
 *  for arrays. A `Scope` names the group of the allocations made while it is alive.
 */
namespace fdtd::core::memory
{
enum Group { fields, material_maps, material_models, cpml, source, monitors, other, group_count };
inline constexpr std::array<const char*, group_count> GROUP_NAMES = {
	"fields", "material maps", "material models", "CPML", "source", "monitors", "other"};

inline std::array<std::size_t, group_count> live{}; ///< bytes held now, per group
inline Group current = other;
inline std::size_t held = 0; ///< bytes held now in all groups, and their maximum so far
inline std::size_t peak = 0; ///< bytes held now in all groups, and their maximum so far
struct Block
{
	Group group;
	std::size_t bytes;
	bool mapped; ///< from mmap, lazy mode, to be unmapped
};
inline std::unordered_map<void*, Block> blocks; ///< block -> its group and size

/** Lazy mode, `fdtd memory`: arrays are address space only. An array costs memory once it is written,
 *  so the set-up of a rank can run, and be counted, on a machine that could not hold its arrays.
 *  The pages read as zero, so the arrays' zero fill is skipped (zero_fill()).
 */
inline bool lazy = false;

/// Zero `count` elements at `first`; nothing to do in lazy mode, fresh pages are zero.
template<class T>
inline void zero_fill(T* first, std::size_t count)
{
	if (!lazy)
		for (std::size_t i = 0; i < count; ++i)
			first[i] = T{};
}

/// Allocations made during the lifetime of a Scope belong to `g`.
struct Scope
{
	Group previous;
	explicit Scope(Group g) : previous(current) { current = g; }
	~Scope() { current = previous; }
	Scope(const Scope&) = delete;
	Scope& operator=(const Scope&) = delete;
};

} // namespace fdtd::core::memory

/**
 * \brief   allocated a memory aligned array
 * \param   align   the byte alignment
 * \param   size    the memory size
 * \returns the memory aligned pointer
 *
 * \note    aligned_alloc() requires the byte count to be a multiple of the
 *          alignment (C11 7.22.3.1); the request is rounded up so that strict
 *          implementations, e.g. macOS, do not return NULL. Zero-sized requests
 *          allocate one alignment unit instead of failing.
 */
template<class T>
static T* palign(size_t align, size_t size)
{
	size_t bytes = size * sizeof(T);
	bytes = ((bytes + align - 1) / align) * align;
	if (bytes == 0)
		bytes = align;

	T* ptr = nullptr;
	const bool mapped = fdtd::core::memory::lazy;
	if (mapped)
	{
		void* p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
		ptr = p == MAP_FAILED ? nullptr : static_cast<T*>(p);
	}
	else
		ptr = static_cast<T*>(aligned_alloc(align, bytes));
	if(ptr==nullptr)
		 throw std::bad_alloc();
	fdtd::core::memory::live[fdtd::core::memory::current] += bytes;
	fdtd::core::memory::held += bytes;
	if (fdtd::core::memory::held > fdtd::core::memory::peak)
		fdtd::core::memory::peak = fdtd::core::memory::held;
	fdtd::core::memory::blocks[ptr] = {fdtd::core::memory::current, bytes, mapped};

	return ptr;
}

/**
 * \brief   frees a memory aligned array allocated with palign
 * \param   ptr   the pointer to free
 */
static inline void pdelete(void * ptr)
{
	if (const auto it = fdtd::core::memory::blocks.find(ptr); it != fdtd::core::memory::blocks.end())
	{
		const fdtd::core::memory::Block block = it->second;
		fdtd::core::memory::live[block.group] -= block.bytes;
		fdtd::core::memory::held -= block.bytes;
		fdtd::core::memory::blocks.erase(it);
		if (block.mapped)
		{
			munmap(ptr, block.bytes);
			return;
		}
	}
	free(ptr);
}

/**
 *  type alias unique_ptrs since the deleter is part of the type
 *  note for shared_ptrs the deleter is only part of the constructor, not the type
 */
///alias for memory aligned unique pointer types
template<class T> using unique_ptr_aligned = std::unique_ptr<T, decltype(&pdelete)>;
///alias for memory aligned unique pointer array types
template<class T> using unique_ptr_aaray = unique_ptr_aligned<T[]>;

/**
 * \brief 	declares a memory aligned unique_ptr
 * \param 	align 	the byte alignment
 * \param 	size 	the memory size
 * \returns a unique_ptr pointing at a memory aligned address
 */
template<class T>
unique_ptr_aligned<T> make_unique_aligned(size_t align, size_t size)
{
	return unique_ptr_aligned<T>(palign<T>(align, size), &pdelete);
}

/**
 * \brief 	declares a memory aligned unique_ptr array
 * \param 	align 	the byte alignment
 * \param 	size 	the memory size
 * \returns a unique_ptr pointing at a memory aligned address
 */
template<class T>
unique_ptr_aaray<T> make_unique_aarray(size_t align, size_t size)
{
	return unique_ptr_aaray<T>(palign<T>(align, size), &pdelete);
}

/**
 * \brief 	declares a memory aligned shared_ptr
 * \param 	align 	the byte alignment
 * \param 	size 	the memory size
 * \returns a shared_ptr pointing at a memory aligned address
 */
template<class T>
std::shared_ptr<T> make_shared_aligned(size_t align, size_t size)
{
	return std::shared_ptr<T>(palign<T>(align, size), &pdelete);
}

/**
 * \brief 	declares a memory aligned shared_ptr array
 * \param 	align 	the byte alignment
 * \param 	size 	the memory size
 * \returns a shared_ptr array pointing at a memory aligned address
 */
template<class T>
std::shared_ptr<T[]> make_shared_aarray(size_t align, size_t size)
{
	return std::shared_ptr<T[]>(palign<T>(align, size), &pdelete);
}

#endif//__ALIGNED_MEMORY_H__
