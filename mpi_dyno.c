#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "decs.h"
#include "list.h"

//---------------------------------------------------------------------
//FUNCTIONS TO INITIALIZE MPI CUSTOM DATA TYPES:

/*	custom MPI Datatype for our vec3 struct */
int init_mpi_vec3(){
	int err;
	//declare the 4 fields required to create a custom MPI Datatype:
	int count; //number of fields in our struct
	int block_lengths[3]; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[3]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[1]; //the different data types included in the struct
	// MPI_Datatype mpi_vec3; //the new custom MPI Datatype
	
	//set count, blocks, and types:
	count = 3;
	block_lengths = {1,1,1};
	types = {MPI_DOUBLE};

	//get the size of an MPI_DOUBLE datatype:
	MPI_Aint size_of_mpi_double;
	err = MPI_Type_extent(MPI_DOUBLE, &size_of_mpi_double);

	//set offsets[]:	
	offsets[0] = static_cast<MPI_Aint>(0);
	offsets[1] = size_of_mpi_double;
	offsets[2] = 2*size_of_mpi_double;

	err = MPI_Type_create_struct(count, block_lengths, offsets, types, mpi_vec3);
	err = MPI_Type_commit(mpi_vec3);
	return err;
}//end init_mpi_vec3()


/*	custom MPI Datatype for our particle struct */
int init_mpi_particle(){
	int err;
	//declare the 4 fields required to create a custom MPI Datatype:
	int count; //number of fields in our struct
	int block_lengths[5]; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[5]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[2]; //the different data types included in the struct
	// MPI_Datatype mpi_grid_point; //the new custom MPI Datatype
	
	//set count, blocks, and types:
	count = 5;
	block_lengths = {1,1,1,1,1};
	types = {mpi_vec3, MPI_DOUBLE};

	//get the size of an mpi_vec3 datatype (our custom made datatype):
	MPI_Aint size_of_mpi_vec3; size_of_mpi_double;
	err = MPI_Type_extent(mpi_vec3, &size_of_mpi_vec3);
	err = MPI_Type_extent(MPI_DOUBLE, &size_of_mpi_double);

	//set offsets[]:	
	offsets[0] = static_cast<MPI_Aint>(0);
	offsets[1] = size_of_mpi_vec3;
	offsets[2] = 2*size_of_mpi_vec3;
	offsets[3] = 2*size_of_mpi_vec3 + size_of_mpi_double;
	offsets[4] = 2*size_of_mpi_vec3 + 2*size_of_mpi_double;

	err = MPI_Type_create_struct(count, block_lengths, offsets, types, mpi_particle);
	err = MPI_Type_commit(mpi_particle);
	return err;
}//end init_mpi_particle()


/*	custom MPI Datatype for our grid_point struct */
int init_mpi_grid_point(){
	int err;
	//declare the 4 fields required to create a custom MPI Datatype:
	int count; //number of fields in our struct
	int block_lengths[2]; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[2]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[1]; //the different data types included in the struct
	// MPI_Datatype mpi_grid_point; //the new custom MPI Datatype
	
	//set count, blocks, and types:
	count = 2;
	block_lengths = {1,1};
	types = {mpi_vec3};

	//get the size of an mpi_vec3 datatype (our custom made datatype):
	MPI_Aint size_of_mpi_vec3;
	err = MPI_Type_extent(mpi_vec3, &size_of_mpi_vec3);

	//set offsets[]:	
	offsets[0] = static_cast<MPI_Aint>(0);
	offsets[1] = size_of_mpi_vec3;

	err = MPI_Type_create_struct(count, block_lengths, offsets, types, mpi_grid_point);
	err = MPI_Type_commit(mpi_grid_point);
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
	err0 = MPI_Type_free (mpi_vec3);
	err1 = MPI_Type_free (mpi_particle);
	err2 = MPI_Type_free (mpi_grid_point);
	err3 = MPI_Type_free (mpi_tree);
	err4 = MPI_Type_free (mpi_tree_node);

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
MPI_Request mpi_list_send(List *p_part_list, int to_pid, particle *part_array) {
	/*	pack particles from List into array (freeing each particle as we pack it...
		this will free the particles, while keeping the list structure in place for
		reuse next time	*/
	int error;
	MPI_Request *req;
	list_reset_iter(p_part_list);
	int i = 0;
	int nParts = list_length(*p_part_list);
	// tell the receiving proc how many particles it will be receiving
	MPI_Isend(&nParts, 1, MPI_INT, to_pid, 1, MPI_COMM_WORLD, req);
	// fill the sendbuffer with particles to send
	while(list_has_next(*p_part_list)) {
		if(i > nParts) {
			printf("ERROR! Exceeded buffer in mpi_list_send function. Exiting.\ni is: %d\n", i);
			MPI_Finalize();
   			exit(-1);
		}//end if
		part_array[i++] = *((particle*) list_get_next(p_part_list));
		list_pop(&part_list);
	}//end while

	/* MPI command to actually send the array of particles (non-blocking) */
	error = MPI_Isend(part_array, i, mpi_particle, to_pid, 1, MPI_COMM_WORLD, req);
	return req;
}//end mpi_list_send

/*	recv a particle list */
MPI_Request* mpi_list_recv(int from_pid, particle* part_array){
	int error, i, tag;
	MPI_Request *req;
	int nParts;
	// recv the # of particles it will be receiving from the sending proc
	MPI_Irecv(&nParts, 1, MPI_INT, from_pid, 1, MPI_COMM_WORLD, req);
	// allocate buffer to receive these particles
	part_array = (particle *)malloc(nParts * sizeof(particle));
	// perform the receive
	error = MPI_Irecv((void*) part_array, nParts, MPI_Particle, from_pid, tag, MPI_COMM_WORLD, req);
	// error = MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[])
	
	// list_reset_iter(part_list);
	// i = 0;

	// //TODO: we actually need to MPI_Waitall before we can do this part below...
	// for( i = 0; i < count; i++){
	// 	list_add(part_list) = &(part_array[i]);
	// }//end for
	return req;
}//end mpi_list_pass

/*	send a tree (and all of it's decendants) */
// MPI_Request mpi_tree_send(tree tree, int to_pid, TreeNode* tree_node_array){
// // error = MPI_Isend(tree_array, i, mpi_???, to_pid, 1, MPI_COMM_WORLD, request);
// 	return req;
// }//end mpi_tree_send



/*	recv a tree (and all of it's decendants) */
// MPI_Request mpi_tree_recv(tree tree, int to_pid, TreeNode* tree_node_array){
// // error = MPI_Irecv((void*) tree_array, i, mpi_???, to_pid, 1, MPI_COMM_WORLD, request);
// 	return req;
// }//end mpi_tree_send


//-----------------------------------------------------------------------------------









// int test_mpi_data_type(){
	 
// 	//MPI custom data type: Struct Derived Datatype: C Example
// 	#define NELEM 25 

// 	int main(argc,argv) 
// 	int argc; 
// 	char *argv[];  { 
// 	int numtasks, rank, source=0, dest, tag=1, i; 

// 	typedef struct { 
// 	  float x, y, z; 
// 	  float velocity; 
// 	  int  n, type; 
// 	  }          Particle; 
// 	Particle     p[NELEM], particles[NELEM]; 
// 	MPI_Datatype particletype, oldtypes[2];  
// 	int          blockcounts[2]; 

// 	/* MPI_Aint type used to be consistent with syntax of */ 
// 	/* MPI_Type_extent routine */ 
// 	MPI_Aint    offsets[2], extent; 

// 	MPI_Status stat; 

// 	MPI_Init(&argc,&argv); 
// 	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
// 	MPI_Comm_size(MPI_COMM_WORLD, &numtasks); 

// 	 /* Setup description of the 4 MPI_FLOAT fields x, y, z, velocity */ 
// 	offsets[0] = 0; 
// 	oldtypes[0] = MPI_FLOAT; 
// 	blockcounts[0] = 4; 
// 	/* Setup description of the 2 MPI_INT fields n, type */ 
// 	/* Need to first figure offset by getting size of MPI_FLOAT */ 
// 	MPI_Type_extent(MPI_FLOAT, &extent); 
// 	offsets[1] = 4 * extent; 
// 	oldtypes[1] = MPI_INT; 
// 	blockcounts[1] = 2; 

// 	/* Now define structured type and commit it */ 

// 	MPI_Type_struct(2, blockcounts, offsets, oldtypes, &particletype); 
// 	MPI_Type_commit(&particletype); 

// 	/* Initialize the particle array and then send it to each task */ 
// 	if (rank == 0) { 
// 	  for (i=0; i < NELEM; i++) { 
// 	     particles[i].x = i * 1.0; 
// 	     particles[i].y = i * -1.0; 
// 	     particles[i].z = i * 1.0;  
// 	     particles[i].velocity = 0.25; 
// 	     particles[i].n = i; 
// 	     particles[i].type = i % 2;  
// 	     } 
	  
// 	 for (i=0; i < numtasks; i++)  
// 	     MPI_Send(particles, NELEM, particletype, i, tag, MPI_COMM_WORLD); 
// 	  } 
	  
// 	MPI_Recv(p, NELEM, particletype, source, tag, MPI_COMM_WORLD, &stat); 

// 	/* Print a sample of what was received */ 
// 	printf("rank= %d   %3.2f %3.2f %3.2f %3.2f %d %d\n", rank,p[3].x, 
// 	     p[3].y,p[3].z,p[3].velocity,p[3].n,p[3].type); 
	    
// 	MPI_Finalize(); 
 
// }//end test_mpi_data_type()
