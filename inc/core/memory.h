/** \file   core/memory.h
 *  \brief  Bytes held in arrays, by group: what `fdtd mesh` reports and `fdtd memory` counts.
 *
 *  Each array charges the blocks it allocates when it is created; the charge is returned when the
 *  last copy releases the memory. A Scope names the group of the charges made while it is alive.
 */
#ifndef FDTD_CORE_MEMORY_H
#define FDTD_CORE_MEMORY_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <unordered_set>

namespace fdtd::core::memory
{
enum Group { fields, material_maps, material_models, cpml, source, monitors, other, group_count };
inline constexpr std::array<const char*, group_count> GROUP_NAMES = {
	"fields", "material_maps", "material_models", "cpml", "source", "monitors", "other"};

inline std::array<std::size_t, group_count> live{}; ///< bytes held now, per group
inline Group current = other; ///< group of the charges made now
inline std::size_t held = 0; ///< bytes held now in all groups
inline std::size_t peak = 0; ///< largest `held` so far

/** Lazy mode, `fdtd memory`: an array costs memory only once written, so the set-up of a rank can
 *  run, and be counted, on a machine that could not hold its arrays. Large blocks then come from
 *  fresh pages that read as zero, and zero_fill() leaves them untouched.
 */
inline bool lazy = false;
inline constexpr std::size_t LAZY_FRESH_BYTES = std::size_t{128} << 10; ///< blocks this large are fresh pages in lazy mode

/// size of one aligned block as an allocator reserves it: rounded up to the alignment, at least one unit
constexpr std::size_t block(std::size_t bytes, std::size_t align)
{
	return std::max(align, (bytes + align - 1) / align * align);
}

/// zero `count` elements; in lazy mode a large block is left as the fresh pages it came from
template <class T>
void zero_fill(T* first, std::size_t count)
{
	if (!lazy || count * sizeof(T) < LAZY_FRESH_BYTES)
		std::fill_n(first, count, T{});
}

/// bytes charged to the current group for the lifetime of the object
class Charge
{
public:
	explicit Charge(std::size_t bytes);
	~Charge();
	Charge(const Charge&) = delete;
	Charge& operator=(const Charge&) = delete;

	Group group() const { return _group; }
	std::size_t bytes() const { return _bytes; }

	/// the charges not yet returned: what an array still holds
	static std::unordered_set<const Charge*>& outstanding();

private:
	Group _group;
	std::size_t _bytes;
};

/// charges made during the lifetime of a Scope belong to `g`
class Scope
{
public:
	explicit Scope(Group g) : _previous(current) { current = g; }
	~Scope() { current = _previous; }
	Scope(const Scope&) = delete;
	Scope& operator=(const Scope&) = delete;

private:
	Group _previous;
};

/// enter lazy mode: large blocks come from fresh pages, see `lazy`
void enter_lazy_mode();

} // namespace fdtd::core::memory

#endif // FDTD_CORE_MEMORY_H
