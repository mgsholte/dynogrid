#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mpi_dyno.h"

//---------------------------------------------------------------------
//MPI SEND AND RECV SUBROUTINES:
MPI_Request* mpi_tree_presend(int *nTreeSends, int *nParticleSends, int to_pid) {
	MPI_Request *reqs = (MPI_Request*) malloc(2*sizeof(MPI_Request));

	MPI_Isend(nTreeSends, 1, MPI_INT, to_pid, TAG_N_SIMPLE_TREES, MPI_COMM_WORLD, &(reqs[0]));
	MPI_Isend(nParticleSends, 1, MPI_INT, to_pid, TAG_N_PARTICLES, MPI_COMM_WORLD, &(reqs[1]));

	return reqs;
}

MPI_Request* mpi_tree_prerecv(int *nTreeRecvs, int *nParticleRecvs, int from_pid) {
	MPI_Request *reqs = (MPI_Request*) malloc(2*sizeof(MPI_Request));

	MPI_Irecv(nTreeRecvs, 1, MPI_INT, from_pid, TAG_N_SIMPLE_TREES, MPI_COMM_WORLD, &(reqs[0]));
	MPI_Irecv(nParticleRecvs, 1, MPI_INT, from_pid, TAG_N_PARTICLES, MPI_COMM_WORLD, &(reqs[1]));

	return reqs;
}
/*	send a list of trees and each tree's particles */
MPI_Request* mpi_tree_send(List* tree_list, int to_pid, simple_tree *simple_trees_array, particle *all_particles_array, int *part_counts, char* dir6){
	int error;

	int trees_len = list_length(tree_list);

	//Allocate an array to hold 3 MPI_Request items to be returned to caller:
	MPI_Request *reqs = (MPI_Request*) malloc(sizeof(MPI_Request)*4);

	int set;
	for (set = 0; set < trees_len; set ++){
		part_counts[set] = 0;
	}
	int all_particles_count = 0;
	
	//convert the list of tree pointers to an array of simple_tree's and pack into an array:
	int i = 0;
	list_reset_iter(tree_list);
	while(list_has_next(tree_list)){
		tree* temp_tree = (tree*) list_get_next(tree_list);
		all_particles_count += list_length(temp_tree->particles);
		simple_trees_array[i].loc = temp_tree->loc;
		simple_trees_array[i].owner = temp_tree->owner;
		simple_trees_array[i].nPoints = temp_tree->nPoints;
		int row,col; //iterating nummbers, no specific meaning
		for(row = 0; row < 3; row++){
			for(col = 0; col < 3; col++){
				simple_trees_array[i].neighbor_owners[row*3 + col] = temp_tree->neighbor_owners[row][col];
			}//end inner of
		}//end outer for
		i++;
	}//end while

	//aggregate and pack all of the particles in the trees_list into one big array of particles:
	list_reset_iter(tree_list);
	i = 0;
	int j, k = 0;
	while(list_has_next(tree_list)){
		List *part_list = ((tree*)list_get_next(tree_list))->particles;
		list_reset_iter(part_list);
		j = 0;
		while(list_has_next(part_list)){
			all_particles_array[i] = *((particle*)list_get_next(part_list));
			list_pop(part_list, true);
			i++;
			j++;	
		}//end inner while
		part_counts[k] = j;
		list_pop(tree_list, false);
		k++;
	}//end while	

	/* MPI command to actually send the array of trees (non-blocking) */
	error = MPI_Isend(simple_trees_array, trees_len, mpi_tree, to_pid, TAG_SIMPLE_TREES, MPI_COMM_WORLD, &(reqs[0]));
	error = MPI_Isend(all_particles_array, all_particles_count, mpi_particle, to_pid, TAG_ALL_PARTICLES, MPI_COMM_WORLD, &(reqs[1]));
	error = MPI_Isend(part_counts, trees_len, MPI_INT, to_pid, TAG_ALL_PARTICLE_COUNTS, MPI_COMM_WORLD, &(reqs[2]));
	error = MPI_Isend(dir6, 1, MPI_CHAR, to_pid, TAG_DIR6, MPI_COMM_WORLD, &(reqs[3]));

	return reqs;
}//end mpi_tree_send



/*	recv an array of trees and each tree's particles and each trees particle_count */
MPI_Request* mpi_tree_recv(int from_pid, simple_tree* simple_trees_array, int nTreeRecvs, particle* all_particles_array, int nParticleRecvs, int* part_counts, char* dir6) {
	
	MPI_Status status;
	MPI_Request *reqs = (MPI_Request*) malloc(sizeof(MPI_Request)*4);

//IRECV THE ARRAY OF SIMPLE_TREES:
	
    // Now receive the message with the allocated buffer (non-blocking, so unpacking must be done in separate function)
	MPI_Irecv(simple_trees_array, nTreeRecvs, mpi_tree, from_pid, TAG_SIMPLE_TREES, MPI_COMM_WORLD, &(reqs[0]));
	
//IRECV THE ARRAY OF PARTICLES:
    // Now receive the message with the allocated buffer (non-blocking, so unpacking must be done in separate function)
	MPI_Irecv(all_particles_array, nParticleRecvs, mpi_particle, from_pid, TAG_ALL_PARTICLES, MPI_COMM_WORLD, &(reqs[1]));
	
//IRECV THE INT ARRAY OF PARTICLE_OFFSETS FOR EACH TREE IN THE LIST OF TREES BEING PASSED:
	// (*buf_lens)[2] = nTreesRecv;//SHOULDN'T NEED THIS ANYMORE

    // Now receive the message with the allocated buffer (non-blocking, so unpacking must be done in separate function)
	MPI_Irecv(part_counts, nTreeRecvs, MPI_INT, from_pid, TAG_ALL_PARTICLE_COUNTS, MPI_COMM_WORLD, &(reqs[2]));
	
//IRECV THE SINGLE CHAR FOR DIR6:
    // Allocate a buffer to hold the incoming simple_trees:
	// dir6 = (char*) malloc(sizeof(char)); // don't think we need to do this, but keeping around just in case

    // Now receive the message with the allocated buffer (non-blocking, so unpacking must be done in separate function)
	MPI_Irecv(dir6, 1, MPI_CHAR, from_pid, TAG_DIR6, MPI_COMM_WORLD, &(reqs[3]));
	
	return reqs;
}//end mpi_tree_recv


List* mpi_tree_unpack(simple_tree *simple_trees_array, particle *all_particles_array, int *part_counts, int trees_len){
	int error;
	List* trees_list = list_init();

	int tot_parts = 0;
	int tree_num,part_num;
	// unpack array of simple_trees back into regular trees and from an array into a list:
	for(tree_num = 0; tree_num < trees_len; tree_num++) {
		simple_tree temp = simple_trees_array[tree_num];
		tree* tree_ptr = tree_init(temp.loc, temp.owner);
		tree_ptr->nPoints = temp.nPoints;
		// tree_ptr->owner = ((*simple_trees_array)[tree_num]).owner;
		
		int row,col; //iterating nummbers, no specific meaning
		for(row = 0; row < 3; row++){
			for(col = 0; col < 3; col++){
				tree_ptr->neighbor_owners[row][col] = temp.neighbor_owners[row*3 + col];
			}//end inner of
		}//end outer for

		list_add(trees_list, tree_ptr);
		for(part_num = 0; part_num < part_counts[tree_num]; part_num++){
			particle *tmp = malloc(sizeof(particle));
			*tmp = all_particles_array[tot_parts];
			list_add(tree_ptr->particles, tmp);
			tot_parts++;
		}//end inner for
	}//end for
	free(simple_trees_array);
	free(all_particles_array);
	free(part_counts);
	return trees_list;
}//end mpi_tree_unpack()
