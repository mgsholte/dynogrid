#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "decs.h"
#include "list.h"


MPI_Request mpi_list_send(List part_list, int to_pid, particle* part_array){
	/*	pack particles from List into array (freeing each particle as we pack it...
		this will free the particles, while keeping the list structure in place for
		reuse next time	*/
	int error;
	list_reset_iter(&part_list);
	int i = 0;

	while(list_has_next(part_list)){
		if(i > 4*part_per_cell){
			printf("ERROR! Exceeded buffer in mpi_list_send function. Exiting.\ni is: %d\n", i);
			MPI_Finalize();
   			exit(-1);
		}//end if
		part_array[i] = *((particle*) list_get_next(&part_list));
		list_pop(&part_list);
		i++;
	}//end while

	/* MPI command to actually send the array of particles (non-blocking) */
	error = MPI_Isend(part_array, i, MPI_Particle, to_pid, 1, MPI_COMM_WORLD, request);
	return req;
}//end mpi_list_pass


MPI_Request mpi_list_recv(List part_list, int from_pid, particle part_array){
	int error, i, count tag;
	MPI_Request *request req;

	error = MPI_Irecv((void*) part_array, count, MPI_Particle, from_pid, tag, MPI_COMM_WORLD, req);
	// error = MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[])
	
	// list_reset_iter(part_list);
	// i = 0;

	// //TODO: we actually need to MPI_Waitall before we can do this part below...
	// for( i = 0; i < count; i++){
	// 	list_add(part_list) = &(part_array[i]);
	// }//end for
	return req;
}//end mpi_list_pass