/** \file   core/memory.cc
 *  \brief  Lazy mode of the memory tally.
 */
#include "core/memory.h"

#include <malloc.h>

using fdtd::core::memory::Charge;

Charge::Charge(std::size_t bytes) : _group(current), _bytes(bytes)
{
	live[_group] += _bytes;
	held += _bytes;
	peak = std::max(peak, held);
	outstanding().insert(this);
}

Charge::~Charge()
{
	live[_group] -= _bytes;
	held -= _bytes;
	outstanding().erase(this);
}

std::unordered_set<const Charge*>& Charge::outstanding()
{
	static std::unordered_set<const Charge*> set;
	return set;
}

void fdtd::core::memory::enter_lazy_mode()
{
	// a fixed threshold: glibc otherwise raises it after large frees and serves later large blocks
	// from reused, written memory
	mallopt(M_MMAP_THRESHOLD, static_cast<int>(LAZY_FRESH_BYTES));
	lazy = true;
}
