#include "grid/state.h"
#include "grid/subdomain.h"
// no shared state

namespace fdtd::grid
{
Subdomain &this_rank()
{
	static Subdomain instance;
	return instance;
}

MPI_Comm member_comm(bool member)
{
	int mine = member ? 1 : 0;
	int all = 0;
	MPI_Allreduce(&mine, &all, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
	if (all)
		return MPI_COMM_WORLD;
	MPI_Comm comm;
	MPI_Comm_split(MPI_COMM_WORLD, member ? 0 : MPI_UNDEFINED, fdtd::grid::world_rank, &comm);
	return comm;
}
} // namespace fdtd::grid
