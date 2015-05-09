#ifndef MPICOMM_H
#define MPICOMM_H

typedef struct {
	// the id of the neighbor processor
	int pid;
	// a list of the particle lists to send to this neighbor
	List part_lists;
	// array of buffers for sending/recving each particle list
	particle **sendbufs, **recvbufs;
} neighbor;

neighbor neighbor_init(int pid);

void neighbor_send(neighbor n);

#endif //MPICOMM_H
