#ifndef MPICOMM_H
#define MPICOMM_H

#include "decs.h"
#include "list.h"
#include "tree.h"

typedef struct {
	// the id of the neighbor processor
	// the # of cells to send/recv aka the lens of the send/recv bufs
	int pid, ncellsends, ncellrecvs;
	// a list of the particle lists to send to this neighbor
	List part_lists;
	// array of buffers for sending/recving each particle list
	particle **sendbufs, **recvbufs;
	// the # of particles to send in a given cell aka sendlen[i] is the length of sendbuf[i], etc.
	int *sendlens, *recvlens;
} neighbor;

neighbor* neighbor_init(int pid);
void neighbor_free(neighbor *n);

void neighbor_add_cell(neighbor *n, tree *cell);

// return the number of total requests to wait on
int neighbor_send_cell_count(neighbor *n, MPI_Request *reqs);
MPI_Request* neighbor_send_cells(neighbor *n);
MPI_Request* neighbor_recv_cells(neighbor *n);

#endif //MPICOMM_H
