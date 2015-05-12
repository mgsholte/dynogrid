#include <stdlib.h>

#include "mpicomm.h"

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
		free(n->recvbufs[iCell]);
	}
	// free send buffers
	for (iCell = 0; iCell < n->ncellsends; ++iCell) {
		free(n->sendbufs[iCell]);
	}
	// free list len buffers
	free(n->recvlens);
	free(n->sendlens);
	// free array storing the buffers
	free(n->sendbufs);
	free(n->recvbufs);
	// free the particle list
	list_free(n->part_lists);
	// free the neighbor itself
	free(n);
}

void neighbor_add_cell(neighbor *n, tree *cell) {
	// only add cells that have particles to send
	if (list_length(cell->particles) != 0) {
		//TODO: could init sendlens here and set its values rather than doing that in send_cell_lengths
		list_add(&(n->part_lists), &(cell->particles));
		n->ncellsends += 1;
	}
}

// always return 2 requests, one for the send, the other for the recv
MPI_Request* neighbor_comm_cell_count(neighbor *n) {
	MPI_Request *reqs = (MPI_Request*) malloc(2 * sizeof(MPI_Request));
	// tell neighbors how many cells you will send them
	MPI_Isend(&(n->ncellsends), 1, MPI_INT, n->pid, TAG_N_CELLS, MPI_COMM_WORLD, &reqs[0]);
	// recv same info back from them
	MPI_Irecv(&(n->ncellrecvs), 1, MPI_INT, n->pid, TAG_N_CELLS, MPI_COMM_WORLD, &reqs[1]);
	return reqs;
}

//TODO: name -> neighbor_comm_cell_lengths since it is both send and recv
// always returns an array of 2 requests
MPI_Request* neighbor_comm_cell_lengths(neighbor *n) {
	MPI_Request *reqs = (MPI_Request*) malloc(2 * sizeof(MPI_Request));

	// need to allocate the array of # of particles to send in each cell
	n->sendlens = (int*) malloc(n->ncellsends * sizeof(int));
	// now fill the entries in the sendlens array with the lengths of the lists to send
	int iCell = 0;
	list_reset_iter(&n->part_lists);
	while(list_has_next(n->part_lists)) {
		// pointer to the particle list (not an array)
		List *curSendList = (List*) list_get_next(&n->part_lists);
		n->sendlens[iCell++] = list_length(*curSendList);
	}
	
	// do the non-blocking send
	MPI_Isend(
		n->sendlens, n->ncellsends, MPI_INT,  // send data
		n->pid, TAG_LIST_LENGTH, MPI_COMM_WORLD, // recver data
		&reqs[0]); // request handle
	
	// need to allocate the array of # of particles to recv in each cell
	n->recvlens = (int*) malloc(n->ncellrecvs * sizeof(int));

	// followed by its non-blocking recv counterpart
	MPI_Irecv(
		n->recvlens, n->ncellrecvs, MPI_INT, // recv data
		n->pid, TAG_LIST_LENGTH, MPI_COMM_WORLD, // sender data
		&reqs[1]); // request handle

	// return the request handles so we can wait for the communication to finish
	return reqs;
}

// 2 isends for every cell being sent to the neighbor
MPI_Request* neighbor_comm_cells(neighbor *n) {
	// allocate the array to hold the buffers for sending particles
	n->sendbufs = (particle **)malloc(n->ncellsends * sizeof(particle*));

	// now loop over each cell and send the particles themselves
	int iCell = 0;
	list_reset_iter(&n->part_lists);
	while(list_has_next(n->part_lists)) {
		// pointer to the particle list (not an array)
		List *curSendList = (List*) list_get_next(&n->part_lists);
		// allocate the buffer used for the send of this cell
		particle *sendbuf = (particle *)malloc(n->sendlens[iCell] * sizeof(particle));
		// populate the buffer with the list contents
		int iPart = 0;
		list_reset_iter(curSendList);
		while(list_has_next(*curSendList)) {
			sendbuf[iPart++] = *((particle*) list_get_next(curSendList));
			// we don't need to keep the particle on this proc anymore so remove it from the list (which also frees it)
			list_pop(curSendList);
		}
		// save handles to send/recv buffers to free later
		n->sendbufs[iCell++] = sendbuf;
	}

	MPI_Request *reqs = (MPI_Request*) malloc((n->ncellsends+n->ncellrecvs) * sizeof(MPI_Request));

	// do the non-blocking send for each cell going to your neighbor
	for (iCell = 0; iCell < n->ncellsends; ++iCell) {
		MPI_Isend(
			n->sendbufs[iCell], n->sendlens[iCell], mpi_particle,  // send data
			n->pid, TAG_PARTICLES, MPI_COMM_WORLD, // recver data
			&reqs[iCell]); // request handle
	}

	// allocate the array to hold the buffers for sending particles
	n->recvbufs = (particle **)malloc(n->ncellrecvs * sizeof(particle*));

	// do the non-blocking recv for each cell sent by your neighbor
	for (iCell = 0; iCell < n->ncellrecvs; ++iCell) {
		// allocate the buffer to receive the contents from your neighbor
		n->recvbufs[iCell] = (particle *)malloc(n->recvlens[iCell] * sizeof(particle));
		MPI_Irecv(
			n->recvbufs[iCell], n->recvlens[iCell], mpi_particle, // recv data
			n->pid, TAG_PARTICLES, MPI_COMM_WORLD, // sender data
			&reqs[iCell+n->ncellsends]); // request handle. offset b/c we add these after the send reqs
	}
		
	return reqs;
}
