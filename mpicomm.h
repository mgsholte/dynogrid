#ifndef MPICOMM_H
#define MPICOMM_H

#include "decs.h"
#include "list.h"

typedef struct {
	// the id of the neighbor processor
	// the # of cells to send/recv
	int pid, numsends, numrecvs;
	// a list of the particle lists to send to this neighbor
	List part_lists;
	// array of buffers for sending/recving each particle list
	particle **sendbufs, **recvbufs;
} neighbor;

neighbor neighbor_init(int pid);

void neighbor_send(neighbor n);

#endif //MPICOMM_H
