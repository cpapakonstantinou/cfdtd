/** \file   fields/fields.cc
 *  \brief  The electromagnetic fields of this rank, inc/fields/fields.h.
 */
#include "fields/fields.h"

#include <new>

namespace fdtd::fields
{
Fields &fields()
{
	static Fields instance;
	return instance;
}

bool Fields::allocate(int nx, int ny, int nz, const Needs &needs)
{
	try
	{
		for (int c = 0; c < 3; ++c)
		{
			for (int set = 0; set < (needs.harmonic ? 2 : 1); ++set)
			{
				E[set][c] = field3<double>(nx, ny, nz);
				H[set][c] = field3<double>(nx, ny, nz);
				E_prev[set][c] = field3<double>(nx, ny, nz);
			}
			if (needs.e_prev_2)
				E_prev_2[c] = field3<double>(nx, ny, nz);
			if (needs.h_prev)
				H_prev[c] = field3<double>(nx, ny, nz);
			if (needs.e_scratch)
				E_scratch[c] = field3<double>(nx, ny, nz);
		}
	}
	catch (const std::bad_alloc &)
	{
		return false;
	}
	return true;
}

void Fields::release()
{
	for (int c = 0; c < 3; ++c)
	{
		for (int set = 0; set < 2; ++set)
			E[set][c].reset(), H[set][c].reset(), E_prev[set][c].reset();
		E_prev_2[c].reset(), H_prev[c].reset(), E_scratch[c].reset();
	}
}
} // namespace fdtd::fields
