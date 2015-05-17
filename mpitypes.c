#include <stdio.h>
#include <mpi.h>

#include "mpitypes.h"

//---------------------------------------------------------------------
//FUNCTIONS TO INITIALIZE MPI CUSTOM DATA TYPES:

/*	custom MPI Datatype for our vec3 struct */
int init_mpi_vec3(){
	int err;
	//declare the 4 fields required to create a custom MPI Datatype:
	int count; //number of fields in our struct
	int block_lengths[3] = {1,1,1}; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[3]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[3] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE}; //the different data types included in the struct
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
	MPI_Datatype types[5] = {mpi_vec3, mpi_vec3, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; //the different data types included in the struct
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
	MPI_Datatype types[2] = {mpi_vec3, mpi_vec3}; //the different data types included in the struct
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


/*	custom MPI Datatype for our tree struct */
int init_mpi_tree(){
	int err;
	//declare the 4 fields required to create a custom MPI Datatype:
	int count; //number of fields in our struct
	int block_lengths[4] = {1,1,1,9}; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[4]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[4] = {mpi_vec3, MPI_INT, MPI_INT, MPI_INT}; //the different data types included in the struct
	// MPI_Datatype mpi_tree; //the new custom MPI Datatype
	
	//set count, blocks, and types:
	count = 4;

	//get the size of an mpi_vec3 datatype (our custom made datatype):
	MPI_Aint size_of_mpi_vec3;
	MPI_Aint size_of_mpi_int;
	err = MPI_Type_extent(mpi_vec3, &size_of_mpi_vec3);
	err = MPI_Type_extent(MPI_INT, &size_of_mpi_int);

	//set offsets[]:	
	offsets[0] = (MPI_Aint) (0);
	offsets[1] = size_of_mpi_vec3;
	offsets[2] = size_of_mpi_vec3 + size_of_mpi_int;
	offsets[3] = size_of_mpi_vec3 + 2*size_of_mpi_int;

	err = MPI_Type_create_struct(count, block_lengths, offsets, types, &mpi_tree);
	err = MPI_Type_commit(&mpi_tree);
	return err;
}//end init_mpi_tree()


/* Initialize ALL mpi custom datatypes */
int init_mpi_customs(){
	int err;
	// //allocate space on the heap to store our custom datatype handles:
	// mpi_vec3 = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));
	// mpi_particle = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));
	// mpi_grid_point = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));
	// mpi_tree = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));
	// mpi_tree_node = (MPI_Datatype*) malloc(sizeof(MPI_Datatype));

	if (pid == 0) {
		printf("initing mpi_types\n");
	}
	//initializes our MPI custom data types:
	err = init_mpi_vec3();
	err = init_mpi_particle();	
	err = init_mpi_grid_point();	
	err = init_mpi_tree();	
	if (pid == 0) {
		printf("finished initing mpi_types\n");
	}

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
	// if(err4 != MPI_SUCCESS){
	// 	return err4;
	// }
	return MPI_SUCCESS;
}//end free_mpi_customs()
