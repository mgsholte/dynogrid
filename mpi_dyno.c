#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mpi_dyno.h"

//---------------------------------------------------------------------
//MPI SEND AND RECV SUBROUTINES:

/*	send a list of trees and each tree's particles */
MPI_Request* mpi_tree_send(List *tree_list, int to_pid, simple_tree** simple_trees_array, particle** all_particles_array, int** part_counts, char* dir6){
	int error;

	//Allocate an array to hold 3 MPI_Request items to be returned to caller:
	MPI_Request *reqs = (MPI_Request*) malloc(sizeof(MPI_Request)*4);

	int trees_len = list_length(tree_list);
	*simple_trees_array = (simple_tree*) malloc(sizeof(simple_tree) * trees_len);
	*part_counts = (int*) malloc(sizeof(int) * trees_len);
	int all_particles_count = 0;
	
	//convert the list of tree pointers to an array of simple_tree's and pack into an array:
	int i = 0;
	list_reset_iter(tree_list);
	while(list_has_next(tree_list)){
		tree *temp_tree = (tree*) list_get_next(tree_list);
		all_particles_count += list_length(temp_tree->particles);
		simple_trees_array[i]->loc = temp_tree->loc;
		simple_trees_array[i]->owner = temp_tree->owner;
		int row,col; //iterating nummbers, no specific meaning
		for(row = 0; row < 3; row++){
			for(col = 0; col < 3; col++){
				simple_trees_array[i]->neighbor_owners[row*3 + col] = temp_tree->neighbor_owners[row][col];
			}//end inner of
		}//end outer for
		i++;
	}//end while

	//aggregate and pack all of the particles in the trees_list into one big array of particles:
	*all_particles_array = (particle*) malloc(sizeof(particle) * all_particles_count);
	list_reset_iter(tree_list);
	i = 0;
	int j, k = 0;
	while(list_has_next(tree_list)){
		List *part_list = ((tree*)list_get_next(tree_list))->particles;
		list_reset_iter(part_list);
		j = 0;
		while(list_has_next(part_list)){
			(*all_particles_array)[i] = *((particle*)list_get_next(part_list));
			list_pop(part_list);
			i++;
			j++;	
		}//end inner while
		(*part_counts)[k] = j;
		list_pop(tree_list);
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


List* mpi_tree_unpack(simple_tree** simple_trees_array, particle** all_particles_array, int** part_counts, int* (buf_lens[])){
	int error;
	List* trees_list = list_init();

	int trees_len = (*buf_lens)[0];

	int tot_parts = 0;
	int tree_num,part_num;
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
		list_add(trees_list, tree_ptr);
		for(part_num = 0; part_num < (*part_counts)[tree_num]; part_num++){
			list_add(tree_ptr->particles, (particle*) all_particles_array[tot_parts]);
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

