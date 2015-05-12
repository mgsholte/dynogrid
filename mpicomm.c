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
		free(n->recvbufs[iCell]);
	}
	// free send buffers
	for (iCell = 0; iCell < n->ncellsends; ++iCell) {
		free(n->sendbufs[iCell]);
	}
	// free array storing the buffers
	if (n->ncellsends != 0) {
		free(n->sendbufs);
	}
	if (n->ncellrecvs != 0) {
		free(n->recvbufs);
	}
	// free the neighbor itself
	free(n);
}

void neighbor_add_cell(neighbor *n, tree *cell) {
	if (list_length(cell->particles)!=0){
		list_add(&(n->part_lists), &(cell->particles));
		//TODO: remove if unnecessary
		//n->ncellsends = list_length(n->part_lists);
		n->ncellsends += 1;
	}
}

// always return 2 requests, one for the send, the other for the recv
MPI_Request* neighbor_send_cell_count(neighbor *n) {
	MPI_Request *reqs = (MPI_Request*) malloc(2*sizeof(MPI_Request));
	// tell neighbors how many cells you will send them
	MPI_Isend(&(n->ncellsends), 1, MPI_INT, n->pid, TAG_N_CELLS, MPI_COMM_WORLD, &reqs[0]);
	// recv same info back from them
	MPI_Irecv(&(n->ncellrecvs), 1, MPI_INT, n->pid, TAG_N_CELLS, MPI_COMM_WORLD, &reqs[1]);
	return reqs;
}

//TODO: name -> neighbor_comm_cell_lengths since it is both send and recv
MPI_Request* neighbor_send_cell_lengths(neighbor *n) {
	MPI_Request *reqs = (MPI_Request*) malloc(2*sizeof(MPI_Request));
	if (n->ncellsends == 0) { // no cells to send so nothing to wait on
		reqs[0] = reqs[1] = MPI_REQUEST_NULL;
		return reqs;
	}

	reqs = (MPI_Request*) malloc(n->ncellsends * sizeof(MPI_Request));
	// need to allocate the array of # of particles to send in each cell
	n->sendlens = (int*) malloc(n->ncellsends * sizeof(int));
	// now fill the entries in the sendlens array with the lengths of the lists to send
	int i = 0;
	list_reset_iter(&n->part_lists);
	while(list_has_next(n->part_lists)) {
		// pointer to the particle list (not an array)
		List *curSendList = (List*) list_get_next(&n->part_lists);
		n->sendlens[i++] = list_length(*curSendList);
	}
	
	// do the non-blocking send
	MPI_Isend(
		n->sendlens, n->ncellsends, MPI_INT,  // send data
		n->pid, TAG_LIST_LENGTH, MPI_COMM_WORLD, // recver data
		&reqs[0]); // request handle

	// followed by its non-blocking recv counterpart
	MPI_Irecv(
		n->recvlens, n->ncellrecvs, MPI_INT, // recv data
		n->pid, TAG_LIST_LENGTH, MPI_COMM_WORLD, // sender data
		&reqs[1]); // request handle

	// return the request handles so we can wait for them to finish
	return reqs;
}

// 2 isends for every cell being sent to the neighbor
MPI_Request* neighbor_send_cells(neighbor *n) {
	MPI_Request *reqs;

	// no communication to do in this case
	if (n->ncellsends == 0 && n->ncellrecvs == 0) {
		return NULL;
	}
	
	int iCell = 0;
	if (n->ncellsends != 0) { 
		// allocate the array to hold the buffers for sending particles
		n->sendbufs = (particle **)malloc(n->ncellsends * sizeof(particle*));

		// now loop over each cell and send the particles themselves
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
	}

	reqs = (MPI_Request*) malloc((n->ncellsends+n->ncellrecvs) * sizeof(MPI_Request));

	// do the non-blocking send for each cell going to your neighbor
	for (iCell = 0; iCell < n->ncellsends; ++iCell) {
		MPI_Isend(
			n->sendbufs[iCell], n->sendlens[iCell], MPI_INT,  // send data
			n->pid, TAG_PARTICLES, MPI_COMM_WORLD, // recver data
			&reqs[iCell]); // request handle
	}

	// do the non-blocking recv for each cell sent by your neighbor
	for (iCell = 0; iCell < n->ncellrecvs; ++iCell) {
		// allocate the buffer to receive the contents from your neighbor
		n->recvbufs[iCell] = (particle *)malloc(n->recvlens[iCell] * sizeof(particle));
		MPI_Irecv(
			n->recvbufs[iCell], n->recvlens[iCell], MPI_INT, // recv data
			n->pid, TAG_PARTICLES, MPI_COMM_WORLD, // sender data
			&reqs[iCell+n->ncellsends]); // request handle. offset b/c we add these after the send reqs
	}
		
	return reqs;
}
