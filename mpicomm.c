#include <mpi.h>

#include "mpicomm.h"

neighbor neighbor_init(int pid) {
	return (neighbor) { pid, -1, -1, list_init(), NULL, NULL };
}

void neighbor_add_cell(neighbor *n, tree *cell) {
	list_add(&(n->part_lists), &(cell->particles));
	n->numsends = list_length(n->part_lists);
}

MPI_Request neighbor_send_cell_count(neighbor n) {
	MPI_Request request;
	// tell neighbors how many cells you will send them
	MPI_Isend(&n.numsends, 1, MPI_INT, n.pid, 1, MPI_COMM_WORLD, &request);
	// recv same info back from them
	MPI_Irecv(&n.numrecvs, 1, MPI_INT, n.pid, 1, MPI_COMM_WORLD, &request);
	return request;
}

MPI_Request* neighbor_send_cells(neighbor n) {
	MPI_Request *requests = (MPI_Request*) calloc( n.numsends*sizeof(MPI_Request) );
	// allocate the buffers for sending/recving particles
	n.sendbufs = (particle **)malloc(n.numsends * sizeof(particle*));
	n.recvbufs = (particle **)malloc(n.numrecvs * sizeof(particle*));
	// now loop over each cell and send the particles themselves
	list_reset_iter(&n.part_lists);
	while(list_has_next(n.part_lists)) {
		// pointer to the particle list (not an array)
		List *curSendList = (List*) list_get_next(&n.part_lists);
		int i;
		for (i = 0; i < n.numsends; ++i) {
			mpi_list_send(curSendList, n, i);
		}

		for (i = 0; i < n.numrecvs; ++i) {
			requests[i] = mpi_list_recv(n, i);
		}
	}
	return requests;
}
