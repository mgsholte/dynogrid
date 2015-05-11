#include <mpi.h>
#include <stdlib.h>

#include "mpicomm.h"
#include "mpi_dyno.h"

neighbor* neighbor_init(int pid) {
	neighbor *n = (neighbor*) malloc(sizeof(neighbor));
	*n = (neighbor) { pid, 0, 0, list_init(), NULL, NULL, NULL, NULL };
	return n;
}

void neighbor_free(neighbor *n) {
	// loop over the cells you received
	int iCell;
	// free recv buffers
	for (iCell = 0; iCell < n->ncellrecvs; ++iCell) {
		//TODO: don't think we need the if statement, just the free
		if (n->recvbufs[iCell]) {
			free(n->recvbufs[iCell]);
			n->recvbufs[iCell] = NULL;
		}
	}
	// free send buffers
	for (iCell = 0; iCell < n->ncellsends; ++iCell) {
		//TODO: don't think we need the if statement, just the free
		if (n->sendbufs[iCell]) {
			free(n->sendbufs[iCell]);
			n->sendbufs[iCell] = NULL;
		}
	}
	// free array storing the buffers
	free(n->sendbufs);
	free(n->recvbufs);
	// free the neighbor itself
	free(n);
}

void neighbor_add_cell(neighbor *n, tree *cell) {
	list_add(&(n->part_lists), &(cell->particles));
	//TODO: remove if unnecessary
	//n->ncellsends = list_length(n->part_lists);
	n->ncellsends += 1;
}

MPI_Request neighbor_send_cell_count(neighbor *n) {
	MPI_Request request;
	// tell neighbors how many cells you will send them
	MPI_Isend(&(n->ncellsends), 1, MPI_INT, n->pid, 1, MPI_COMM_WORLD, &request);
	// recv same info back from them
	MPI_Irecv(&(n->ncellrecvs), 1, MPI_INT, n->pid, 1, MPI_COMM_WORLD, &request);
	return request;
}

void neighbor_send_cells(neighbor *n) {
	//MPI_Request *requests = (MPI_Request*) calloc(n.ncellsends, sizeof(MPI_Request));
	// allocate the buffers for sending particles
	n->sendbufs = (particle **)malloc(n->ncellsends * sizeof(particle*));
	// need to allocate array of # of particles to send in each cell
	n->sendlens = (int*) malloc( n->ncellsends*sizeof(int) );
	// now loop over each cell and send the particles themselves
	int i = 0;
	list_reset_iter(&n->part_lists);
	while(list_has_next(n->part_lists)) {
		// pointer to the particle list (not an array)
		List *curSendList = (List*) list_get_next(&n->part_lists);
		
		mpi_list_send(*curSendList, n, i);
		//requests[i] = mpi_list_recv(n, i);
		++i;
	}
}

void neighbor_recv_cells(neighbor *n) {
	// allocate the buffers for recving particles
	n->recvbufs = (particle **)malloc(n->ncellrecvs * sizeof(particle*));
	// need to allocate array of # of particles to recv in each cell
	n->recvlens = (int*) malloc( n->ncellrecvs*sizeof(int) );
	int i = 0;
	list_reset_iter(&n->part_lists);
	while(list_has_next(n->part_lists)) {
		// pointer to the particle list (not an array)
		List *curSendList = (List*) list_get_next(&n->part_lists);
		
		mpi_list_recv(n, i);
		++i;
	}
}
