#ifndef MPICOMM_H
#define MPICOMM_H

#include "decs.h"
#include "list.h"
#include "tree.h"

typedef struct {
	// the id of the neighbor processor
	int pid;
	// the # of cells to send/recv aka the lens of the send/recv bufs
	// and, by extension, the lens of the sendlens/recvlens arrays
	int ncellsends, ncellrecvs;
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

// give to (and recv from) your neighbor a count of how many cells will be communicated in the next steps
MPI_Request* neighbor_send_cell_count(neighbor *n);

// give to (and recv from) your neighbor an array of how many particles are in each cell being communicated
MPI_Request* neighbor_send_cell_lengths(neighbor *n);

// give to (and recv from) your neigbor an array of all the arrays of particles in a cell -- one for each cell
MPI_Request* neighbor_send_cells(neighbor *n);

#endif //MPICOMM_H
