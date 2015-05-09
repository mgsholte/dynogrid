#ifndef MPICOMM_H
#define MPICOMM_H

#include "decs.h"
#include "list.h"

typedef struct {
	// the id of the neighbor processor
	// the # of cells to send/recv
	int pid, ncellsends, ncellrecvs;
	// a list of the particle lists to send to this neighbor
	List part_lists;
	// array of buffers for sending/recving each particle list
	particle **sendbufs, **recvbufs;
	// the lengths of each of the send/recv buffers
	int *sendlens, *recvlens;
} neighbor;

neighbor neighbor_init(int pid);

void neighbor_add_cell(neighbor *n, tree *cell);

MPI_Request neighbor_send_cell_count(neighbor n);
MPI_Request* neighbor_send_cells(neighbor n);

#endif //MPICOMM_H
