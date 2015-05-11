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
	int block_lengths[3] = {1,1,9}; //the number of items in each block in our struct (e.g. arrays would have blockcounts of len(array))
	MPI_Aint offsets[3]; //the offset of the start of each block in the struct, relative to the start of the struct (i.e. offset[0] = 0)
	MPI_Datatype types[3] = {mpi_vec3, MPI_INT, MPI_INT}; //the different data types included in the struct
	// MPI_Datatype mpi_tree; //the new custom MPI Datatype
	
	//set count, blocks, and types:
	count = 3;

	//get the size of an mpi_vec3 datatype (our custom made datatype):
	MPI_Aint size_of_mpi_vec3;
	MPI_Aint size_of_mpi_int;
	err = MPI_Type_extent(mpi_vec3, &size_of_mpi_vec3);
	err = MPI_Type_extent(MPI_INT, &size_of_mpi_int);

	//set offsets[]:	
	offsets[0] = (MPI_Aint) (0);
	offsets[1] = size_of_mpi_vec3;
	offsets[2] = size_of_mpi_vec3 + size_of_mpi_int;

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

//---------------------------------------------------------------------
//MPI SEND AND RECV SUBROUTINES:

/*	send a particle list */
// 'n' is the neighboring proc where we are sending the particle list 'part_list' of cell 'iCell'
MPI_Request* mpi_list_send(List part_list, neighbor *n, int iCell) {
	/*	pack particles from List into array (freeing each particle as we pack it...
		this will free the particles, while keeping the list structure in place for
		reuse next time	*/
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
	n->sendbufs[iCell] = buf;
	n->sendlens[iCell] = nParts;

	MPI_Request *reqs = (MPI_Request*) malloc(2*sizeof(MPI_Request));
	/* MPI command to send the number of particles in the array (non-blocking) */
	MPI_Isend(&n->sendlens[iCell], 1, MPI_INT, n->pid, TAG_LIST_LENGTH, MPI_COMM_WORLD, &reqs[0]);

	/* MPI command to actually send the array of particles (non-blocking) */
	MPI_Isend(buf, nParts, mpi_particle, n->pid, TAG_PARTICLES, MPI_COMM_WORLD, &reqs[1]);

	return reqs;
}//end mpi_list_send

/*	recv a particle list */
// iCell is the index of the cell being received in neighbors recvbuffer
MPI_Request mpi_list_recv(neighbor *n, int iCell) {
	int nParts;
	MPI_Request req;

	// recv the # of particles it will be receiving from the sending proc
	MPI_Recv(&nParts, 1, MPI_INT, n->pid, TAG_LIST_LENGTH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	/* TODO: could these replace the above irecv?
	// determine # of particles that were sent so we can alloc recvbuffer
	MPI_Status status;
	MPI_Improbe(n.pid, tag, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, mpi_particle, &nParts);
	*/

	// allocate buffer to receive these particles
	// and store it and its length in 
	n->recvbufs[iCell] = (particle *)malloc(nParts * sizeof(particle));
	n->recvlens[iCell] = nParts;

	// receive the particles themselves
	MPI_Irecv((void*) n->recvbufs[iCell], nParts, mpi_particle, n->pid, TAG_PARTICLES, MPI_COMM_WORLD, &req);
	
	return req;
}//end mpi_list_pass

/*	send a list of trees and each tree's particles */
MPI_Request* mpi_tree_send(List tree_list, int to_pid, simple_tree** simple_trees_array, particle** all_particles_array, int** part_counts, char* dir6){
	int error;

	//Allocate an array to hold 3 MPI_Request items to be returned to caller:
	MPI_Request *reqs = (MPI_Request*) malloc(sizeof(MPI_Request)*4);

	int trees_len = list_length(tree_list);
	*simple_trees_array = (simple_tree*) malloc(sizeof(simple_tree) * trees_len);
	*part_counts = (int*) malloc(sizeof(int) * trees_len);
	int all_particles_count;
	
	//convert the list of tree pointers to an array of simple_tree's and pack into an array:
	tree temp_tree;
	int i = 0;
	list_reset_iter(&tree_list);
	while(list_has_next(tree_list)){
		temp_tree = *((tree*) list_get_next(&tree_list));
		all_particles_count += list_length(temp_tree.particles);
		simple_trees_array[i]->loc = temp_tree.loc;
		simple_trees_array[i]->owner = temp_tree.owner;
		int row,col; //iterating nummbers, no specific meaning
		for(row = 0; row < 3; row++){
			for(col = 0; col < 3; col++){
				simple_trees_array[i]->neighbor_owners[row*3 + col] = temp_tree.neighbor_owners[row][col];
			}//end inner of
		}//end outer for
		i++;
	}//end while

	//aggregate and pack all of the particles in the trees_list into one big array of particles:
	*all_particles_array = (particle*) malloc(sizeof(particle) * all_particles_count);
	list_reset_iter(&tree_list);
	i = 0;
	int j, k = 0;
	while(list_has_next(tree_list)){
		List part_list = ((tree*)list_get_next(&tree_list))->particles;
		list_reset_iter(&part_list);
		j = 0;
		while(list_has_next(part_list)){
			(*all_particles_array)[i] = *((particle*)list_get_next(&part_list));
			list_pop(&part_list);
			i++;
			j++;	
		}//end inner while
		(*part_counts)[k] = j;
		list_pop(&tree_list);
		k++;
	}//end while	

	/* MPI command to actually send the array of trees (non-blocking) */
	error = MPI_Isend(simple_trees_array, trees_len, mpi_tree, to_pid, 0, MPI_COMM_WORLD, &(reqs[0]));
	error = MPI_Isend(all_particles_array, all_particles_count, mpi_particle, to_pid, 1, MPI_COMM_WORLD, &(reqs[1]));
	error = MPI_Isend(*part_counts, trees_len, MPI_INT, to_pid, 2, MPI_COMM_WORLD, &(reqs[2]));
	error = MPI_Isend(dir6, 1, MPI_CHAR, to_pid, 3, MPI_COMM_WORLD, &(reqs[3]));

	return reqs;
}//end mpi_tree_send



/*	recv an array of trees and each tree's particles and each trees particle_count */
MPI_Request* mpi_tree_recv(int from_pid, simple_tree** simple_trees_array, particle** all_particles_array, int** part_counts, int* (buf_lens[]), char* dir6) {
	
	MPI_Status status;
	MPI_Request *reqs = (MPI_Request*) malloc(sizeof(MPI_Request)*4);

//IRECV THE ARRAY OF SIMPLE_TREES:
	int trees_len;
	//Allocate an array to hold 2 MPI_Request items to be returned to caller:
    // Probe for an incoming message:
    MPI_Probe(from_pid, 0, MPI_COMM_WORLD, &status);

    // When probe returns, the status object has the size and other
    // attributes of the incoming message. Get the message size
    MPI_Get_count(&status, MPI_INT, &trees_len);
    (*buf_lens)[0] = trees_len;

    // Allocate a buffer to hold the incoming simple_trees:
	*simple_trees_array = (simple_tree*) malloc(sizeof(simple_tree) * trees_len);

    // Now receive the message with the allocated buffer (non-blocking, so unpacking must be done in separate function)
	MPI_Irecv(simple_trees_array, trees_len, mpi_tree, from_pid, 0, MPI_COMM_WORLD, &(reqs[0]));
	
//IRECV THE ARRAY OF PARTICLES:
    int all_particles_count;
    MPI_Probe(from_pid, 1, MPI_COMM_WORLD, &status);

    // When probe returns, the status object has the size and other
    // attributes of the incoming message. Get the message size
    MPI_Get_count(&status, MPI_INT, &all_particles_count);
	(*buf_lens)[1] = all_particles_count;

    // Allocate a buffer to hold the incoming simple_trees:
	*all_particles_array = (particle*) malloc(sizeof(particle) * all_particles_count);

    // Now receive the message with the allocated buffer (non-blocking, so unpacking must be done in separate function)
	MPI_Irecv(all_particles_array, all_particles_count, mpi_particle, from_pid, 1, MPI_COMM_WORLD, &(reqs[1]));
	
//IRECV THE INT ARRAY OF PARTICLE_OFFSETS FOR EACH TREE IN THE LIST OF TREES BEING PASSED:
	(*buf_lens)[2] = trees_len;

    // Allocate a buffer to hold the incoming simple_trees:
	*part_counts = (int*) malloc(sizeof(int) * trees_len);

    // Now receive the message with the allocated buffer (non-blocking, so unpacking must be done in separate function)
	MPI_Irecv(*part_counts, trees_len, MPI_INT, from_pid, 2, MPI_COMM_WORLD, &(reqs[2]));
	
//IRECV THE SINGLE CHAR FOR DIR6:
    // Allocate a buffer to hold the incoming simple_trees:
	// dir6 = (char*) malloc(sizeof(char)); // don't think we need to do this, but keeping around just in case

    // Now receive the message with the allocated buffer (non-blocking, so unpacking must be done in separate function)
	MPI_Irecv(dir6, 1, MPI_CHAR, from_pid, 3, MPI_COMM_WORLD, &(reqs[3]));
	
	return reqs;
}//end mpi_tree_recv


List mpi_tree_unpack(simple_tree** simple_trees_array, particle** all_particles_array, int** part_counts, int* (buf_lens[])){
	int error;
	List trees_list = list_init();

	int trees_len = (*buf_lens)[0];

	int tree_num,part_num,tot_parts;
	// unpack array of simple_trees back into regular trees and from an array into a list:
	for(tree_num = 0; tree_num < trees_len; tree_num++) {
		tree* tree_ptr = (tree*) malloc(sizeof(tree));
		tree_ptr->root = (TreeNode*) malloc(sizeof(TreeNode));
		tree_ptr->loc = ((*simple_trees_array)[tree_num]).loc;
		tree_ptr->owner = ((*simple_trees_array)[tree_num]).owner;
		
		int row,col; //iterating nummbers, no specific meaning
		for(row = 0; row < 3; row++){
			for(col = 0; col < 3; col++){
				tree_ptr->neighbor_owners[row][col] = ((*simple_trees_array)[tree_num]).neighbor_owners[row*3 + col];
			}//end inner of
		}//end outer for

		tree_ptr->particles = list_init();
		tree_ptr->new_particles = list_init();
		list_add(&trees_list, tree_ptr);
		for(part_num = 0; part_num < (*part_counts)[tree_num]; part_num++){
			list_add(&tree_ptr->particles, (particle*) all_particles_array[tot_parts]);
			tot_parts++;
		}//end inner for
	}//end for
	//free(*simple_trees_array);
	return trees_list;
}//end mpi_tree_unpack()


//-----------------------------------------------------------------------------------




	// for(i = 0; i < trees_len; i++){
	// 	(trees_array[i]).loc = (trees[i])->loc;
	// 	(trees_array[i]).owner = (trees[i])->owner; 
	// }//end for





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

