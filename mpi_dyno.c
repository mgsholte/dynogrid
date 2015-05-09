#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mpi_dyno.h"
#include "list.h"

//---------------------------------------------------------------------
//FUNCTIONS TO INITIALIZE MPI CUSTOM DATA TYPES:

/*	custom MPI Datatype for our vec3 struct */
int init_mpi_vec3(){
	int err;
	//declare the 4 fields required to create a custom MPI Datatype:
	int count; //number of fields in our struct
	int block_lengths[3] = {1,1,1}; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[3]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[1] = {MPI_DOUBLE}; //the different data types included in the struct
	// MPI_Datatype mpi_vec3; //the new custom MPI Datatype
	
	//set count, blocks, and types:
	count = 3;

	//get the size of an MPI_DOUBLE datatype:
	MPI_Aint size_of_mpi_double;
	err = MPI_Type_extent(MPI_DOUBLE, &size_of_mpi_double);

	//set offsets[]:	
	offsets[0] = (MPI_Aint) (0);
	offsets[1] = size_of_mpi_double;
	offsets[2] = 2*size_of_mpi_double;

	err = MPI_Type_create_struct(count, block_lengths, offsets, types, &mpi_vec3);
	err = MPI_Type_commit(&mpi_vec3);
	return err;
}//end init_mpi_vec3()


/*	custom MPI Datatype for our particle struct */
int init_mpi_particle(){
	int err;
	//declare the 4 fields required to create a custom MPI Datatype:
	int count; //number of fields in our struct
	int block_lengths[5] = {1,1,1,1,1}; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[5]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[2] = {mpi_vec3, MPI_DOUBLE}; //the different data types included in the struct
	// MPI_Datatype mpi_grid_point; //the new custom MPI Datatype
	
	//set count, blocks, and types:
	count = 5;

	//get the size of an mpi_vec3 datatype (our custom made datatype):
	MPI_Aint size_of_mpi_vec3, size_of_mpi_double;
	err = MPI_Type_extent(mpi_vec3, &size_of_mpi_vec3);
	err = MPI_Type_extent(MPI_DOUBLE, &size_of_mpi_double);

	//set offsets[]:	
	offsets[0] = (MPI_Aint) (0);
	offsets[1] = size_of_mpi_vec3;
	offsets[2] = 2*size_of_mpi_vec3;
	offsets[3] = 2*size_of_mpi_vec3 + size_of_mpi_double;
	offsets[4] = 2*size_of_mpi_vec3 + 2*size_of_mpi_double;

	err = MPI_Type_create_struct(count, block_lengths, offsets, types, &mpi_particle);
	err = MPI_Type_commit(&mpi_particle);
	return err;
}//end init_mpi_particle()


/*	custom MPI Datatype for our grid_point struct */
int init_mpi_grid_point(){
	int err;
	//declare the 4 fields required to create a custom MPI Datatype:
	int count; //number of fields in our struct
	int block_lengths[2] = {1,1}; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[2]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[1] = {mpi_vec3}; //the different data types included in the struct
	// MPI_Datatype mpi_grid_point; //the new custom MPI Datatype
	
	//set count, blocks, and types:
	count = 2;

	//get the size of an mpi_vec3 datatype (our custom made datatype):
	MPI_Aint size_of_mpi_vec3;
	err = MPI_Type_extent(mpi_vec3, &size_of_mpi_vec3);

	//set offsets[]:	
	offsets[0] = (MPI_Aint) (0);
	offsets[1] = size_of_mpi_vec3;

	err = MPI_Type_create_struct(count, block_lengths, offsets, types, &mpi_grid_point);
	err = MPI_Type_commit(&mpi_grid_point);
	return err;
}//end init_mpi_grid_point()


/* Initialize ALL mpi custom datatypes */
int init_mpi_customs(){
	int err;
	// //allocate space on the heap to store our custom datatype handles:
	// mpi_vec3 = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));
	// mpi_particle = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));
	// mpi_grid_point = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));
	// mpi_tree = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));
	// mpi_tree_node = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));

	//initializes our MPI custom data types:
	err = init_mpi_vec3();
	err = init_mpi_particle();	
	err = init_mpi_grid_point();	

	return err;
}//end init_mpi_custom()


/*	deallocate our custom MPI datatypes and free the global pointers
	to our custom datatype handles */
int free_mpi_customs(){
	int err0, err1, err2, err3, err4;
	err0 = MPI_Type_free (&mpi_vec3);
	err1 = MPI_Type_free (&mpi_particle);
	err2 = MPI_Type_free (&mpi_grid_point);
	err3 = MPI_Type_free (&mpi_tree);
	err4 = MPI_Type_free (&mpi_tree_node);

	// //free global pointers to custom MPI Datatype handles:
	// free(mpi_vec3);
	// free(mpi_particle);
	// free(mpi_grid_point);
	// free(mpi_tree);
	// free(mpi_tree_node);

	//check for errors:
	if(err0 != MPI_SUCCESS){
		return err0;
	}
	if(err1 != MPI_SUCCESS){
		return err1;
	}
	if(err2 != MPI_SUCCESS){
		return err2;
	}
	if(err3 != MPI_SUCCESS){
		return err3;
	}
	if(err4 != MPI_SUCCESS){
		return err4;
	}
	return MPI_SUCCESS;
}//end free_mpi_customs()

//---------------------------------------------------------------------
//MPI SEND AND RECV SUBROUTINES:

/*	send a particle list */
// 'n' is the neighboring proc where we are sending the particle list 'part_list' of cell 'iCell'
MPI_Request mpi_list_send(List part_list, neighbor n, int iCell) {
	/*	pack particles from List into array (freeing each particle as we pack it...
		this will free the particles, while keeping the list structure in place for
		reuse next time	*/
	MPI_Request req;
	int nParts = list_length(part_list);
	particle *buf = (particle *)malloc(nParts * sizeof(particle));
	
	/* TODO: we don't need to send this?
	// tell the receiving proc how many particles it will be receiving
	MPI_Isend(&nParts, 1, MPI_INT, to_pid, 1, MPI_COMM_WORLD, req);
	*/

	int i = 0;
	list_reset_iter(&part_list);
	// fill the sendbuffer with particles to send
	while(list_has_next(part_list)) {
		if(i > nParts) { //TODO: will never happen if list_length works properly
			printf("ERROR! Exceeded buffer in mpi_list_send function. Exiting.\ni is: %d\n", i);
			MPI_Finalize();
   			exit(-1);
		}//end if
		buf[i++] = *((particle*) list_get_next(&part_list));
		list_pop(&part_list);
	}//end while
	// save the handle to the sendbuffer so we can free it later
	n.sendbufs[iCell] = buf;
	n.sendlens[iCell] = nParts;

	/* MPI command to actually send the array of particles (non-blocking) */
	MPI_Isend(buf, nParts, mpi_particle, n.pid, 1, MPI_COMM_WORLD, &req);
	return req;
}//end mpi_list_send

/*	recv a particle list */
// iCell is the index of the cell being received in neighbors recvbuffer
MPI_Request mpi_list_recv(neighbor n, int iCell) {
	int nParts, tag = 0;
	MPI_Request req;

	// recv the # of particles it will be receiving from the sending proc
	MPI_Recv(&nParts, 1, MPI_INT, n.pid, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	/* TODO: could these replace the above irecv?
	// determine # of particles that were sent so we can alloc recvbuffer
	MPI_Status status;
	MPI_Improbe(n.pid, tag, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, mpi_particle, &nParts);
	*/

	// allocate buffer to receive these particles
	n.recvbufs[iCell] = (particle *)malloc(nParts * sizeof(particle));
	// receive the particles themselves
	MPI_Irecv((void*) n.recvbufs[iCell], nParts, mpi_particle, n.pid, tag, MPI_COMM_WORLD, &req);
	
	return req;
}//end mpi_list_pass
