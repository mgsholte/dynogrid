#ifndef MPI_DYNO_H
#define MPI_DYNO_H

#include <mpi.h>

#include "decs.h"
#include "list.h"
#include "mpicomm.h"

int init_mpi_vec3();
int init_mpi_particle();
int init_mpi_grid_point();

int init_mpi_customs();
int free_mpi_customs();

MPI_Request mpi_list_send(List part_list, neighbor n, int iCell);
MPI_Request mpi_list_recv(neighbor n, int iCell);

#endif //MPI_DYNO_H
