#include <mpi.h>

#include "mpicomm.h"

neighbor neighbor_init(int pid) {
	return (neighbor) { pid, -1, -1, list_init(), NULL, NULL };
}

MPI_Request* neighbor_send_cell_count(neighbor n) {
	MPI_Request *request;
	// tell neighbors how many cells you will send them
	MPI_Isend(&n.numsends, 1, MPI_INT, n.pid, 1, MPI_COMM_WORLD, request);
	// recv same info back from them
	MPI_Irecv(&n.numrecvs, 1, MPI_INT, n.pid, 1, MPI_COMM_WORLD, request);
	return request;
}

MPI_Request* neighbor_send_cell(neighbor n) {
	MPI_Request *request;
	int numsends, numrecvs;
	numsends = list_length(n.part_lists);
	// allocate the buffers for sending/recving particles
	n.sendbufs = (particle **)malloc(numsends * sizeof(particle*));
	n.recvbufs = (particle **)malloc(numrecv * sizeof(particle *));
	// now loop over each cell and send the particles themselves
	list_reset_iter(&n.part_lists);
	while(list_has_next(n.part_lists)) {
		// pointer to the particle list (not an array)
		List *curSendList = (List*) list_get_next(&n.part_lists);
		int i;
		for (i = 0; i < numsend; ++i) {
			// allocate sendbuffer. mpi_list_send will populate and send it
			n.sendbufs[i] = (particle *)malloc(nParts * sizeof(particle));
			mpi_list_send(curSendList, n.pid, n.sendbufs[i]);
		}

		for (i = 0; i < numrecv; ++i) {
			// don't allocate recvbuffer, mpi_list_recv will do that
			request = mpi_list_recv(n.pid, n.recvbuffs[i]);
		}
	}
	return request;
}
