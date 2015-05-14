#ifndef MPI_DYNO_H
#define MPI_DYNO_H

#include <mpi.h>

#include "decs.h"
#include "list.h"

MPI_Request* mpi_tree_send(List *tree_list, int to_pid, simple_tree** simple_trees_array, particle** all_particles_array, int** part_counts, char* dir6);

MPI_Request* mpi_tree_recv(int from_pid, simple_tree** simple_trees_array, particle** all_particles_array, int** part_counts, int* (buf_lens[]), char* dir6);

List* mpi_tree_unpack(simple_tree** simple_trees_array, particle** all_particles_array, int** part_counts, int* (buf_lens[]));

#endif //MPI_DYNO_H
