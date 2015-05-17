#include "balance.h"
#include <stdio.h>


bool is_real(tree* curCell){
	if (curCell != NULL) {
		if (curCell->owner == pid){
			return true;
		}
	}
	return false;
}//end is_real()

bool is_ghost(tree* curCell){
	if (curCell != NULL) {
		if (curCell->owner != pid && curCell->owner != -1){
			return true;
		}
	}
	return false;
}//end is_ghost()

void determine_neighbor_matchings(List* ne_matchings[], char dim, tree**** grid){
	// loop over each ghost cell.
	// loop over entire grid to find the ghost cells, this is the easiest way to find them
	// this can't be integrated with above identical loop since the push must be completed before checking to see which pushed things need to be sent to neighbors
	int i,j,k, ne_num; //ne_num is for "neighbor_number"
	
	tree *prevCell = NULL;
	tree *curCell = NULL;
	tree *nextCell = NULL;

	//instantiate all pointers to NULL for checking to see if whe need to malloc them in loop:
	for(ne_num = 0; ne_num < nProcs; ne_num++){
		ne_matchings[ne_num] = NULL;
	}//end for

	/*	determine surface matchings (i.e. list of cells to pass per neighbor)
		in the x (same as i) direction */
	int neighbor, dir = 0;
	int layer = -1;
	surface* prev_surf_match = NULL;
	
	if(dim == 'x'){
		for (i = imin; i < imax; ++i) {
			for (j = jmin; j < jmax; ++j) {
				for (k = kmin; k < kmax; ++k) {
					curCell = grid[i][j][k];
					if(is_ghost(curCell)) {
						//determine the direction of the neighbor (for this particular cell!!!)
						if(i == 0){
							/* 	we can't access the i-1 element of the grid, but we know that
								the neighbor must be to our LEFT direction (hence the dir = "-1") */
							dir = -1; //dir for "direction"
							nextCell = grid[i+1][j][k];
						}else if(i == wi-1){
							/* 	we can't access the i+1 element of the grid, but we know that
								the neighbor must be to our RIGHT direction (hence the dir ="+1") */
							dir = 1;
							prevCell = grid[i-1][j][k];
						}else{
							//we know we can access the i-1 and i+1 element of the grid
							prevCell = grid[i-1][j][k];
							nextCell = grid[i+1][j][k];
						
							if(is_real(nextCell)){
								dir = -1;
							}else if(is_real(prevCell)){
								dir = 1;
							}//end if else if
						}//end if else for checking extremes
						//set the neighbor pid:
						neighbor = curCell->owner;
						//set the layer value for the surface struct (i.e. )
						layer = i;
						
						surface* surf_match;
						List* list_of_cells; 

						if(prev_surf_match != NULL){
							if(	prev_surf_match->neighbor == neighbor &&
								prev_surf_match->layer == layer &&
								prev_surf_match->dir == dir){
								//then we are dealing with the same surface matching
								surf_match = prev_surf_match;
								list_of_cells = surf_match->cells;
							}//end inner if
						}else{
							surf_match = (surface*) malloc(sizeof(surface));
							list_of_cells = list_init();
						}//end if else
						
						surf_match->neighbor = curCell->owner;
						surf_match->dir = dir;
						surf_match->layer = layer;
						
						if(dir < 0){
							list_add(list_of_cells, nextCell);
						}else if(dir > 0){
							list_add(list_of_cells, prevCell);
						}
						if((ne_matchings[surf_match->neighbor]) == NULL){
							//then the List has not yet been initialized so we have to malloc it and initialize it:
							ne_matchings[surf_match->neighbor] = list_init();
						}
						/* reminder: 
								ne_matchings ......	an array of
													List pointers
									where each list holds "surface"s (a new type of struct...see Balance.h)
									each "surface" has meta data and
									a List of trees */
						list_add(ne_matchings[surf_match->neighbor], (void*) surf_match);
						prev_surf_match = surf_match;
						/* reminder ne_matchings[match->neighbor] is a List* */
					}//end if(is_ghost())
				}//end k loop
			}//end j loop
		}//end i loop
	}//end x dimension

	//TODO: y dimension (should be almost identical to x except for changing order of loops)

	//TODO: z dimension (should be almost identical to x except for changing order of loops)

	//for each neighbor, go through its list of surface matchings and combine contiguous matchings together:
	//note, this means that the surface matchings get combined only if they have the same neighbor, direction, AND layer
	surface* cur_surface;
	surface* next_surface;
	for(ne_num = 0; ne_num < nProcs; ne_num++){
		if(ne_matchings[ne_num] == NULL){
			continue;
		}//end if
		list_reset_iter(ne_matchings[ne_num]);
		int it, ct = 0;
		while(ct < list_length(ne_matchings[ne_num])){
			for(it = 0; it < ct; it++){
				//cycle through the list to get to the correct starting point...
				//probably a better way to do this whole scheme, but not sure...
				list_get_next(ne_matchings[ne_num]);
			}//end for
			cur_surface = (surface*) (list_get_next(ne_matchings[ne_num]));
			while(list_has_next(ne_matchings[ne_num])){
				next_surface = (surface*) (list_get_next(ne_matchings[ne_num]));
				if(	(cur_surface->layer == next_surface->layer) &&
					(cur_surface->dir == next_surface->dir)){
					//then the surface matchings should belong to the same surface, thus should be combined:
					//note, we already know that they have the same neighbor at this point bc of indexing into the array...
					//so, combine next_surface into cur_surface and then free next_surface:
					list_combine(cur_surface->cells, next_surface->cells);
					list_pop(ne_matchings[ne_num], true); //remove (and free) next_surface from the list
				}//end if
			}//end inner while
			ct++;
			list_reset_iter(ne_matchings[ne_num]);
		}//end outer while
	}//end for

}//end determine_neighbor_matchings(List** ne_matchings, char dim)


void Balance(tree ****grid){
	tree *curCell = NULL;
	int partwork = 0;
	int cellwork = 0;
	int i,j,k;
	double work, avgwork, mostwork, propensity;
	// Calculate amount of work on this processor
	for (i = imin; i < imax; ++i) {
		for (j = jmin; j < jmax; ++j) {
			for (k = kmin; k < kmax; ++k) {
				curCell = grid[i][j][k];
				if (curCell != NULL) {
					// if (curCell->owner == pid) {
					partwork += list_length(curCell->particles);
					// }
					// TODO: Keep track of number oc decendants
					cellwork += 1;// curCell->descendants;
				}
			}
		}
	}
	work = .8*partwork + .2*cellwork;

	// Find out how much work there is in total:
	MPI_Allreduce(&work, &avgwork, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	avgwork = avgwork / (double)nProcs;
	// Find out the most work that a processor in the group has:
	MPI_Allreduce(&work, &mostwork, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	// List (*cells_to_send)[nProcs];
	// IS THIS CORRECT SYNTAX?? (below, not above)
	List* ne_matchings[nProcs];
	determine_neighbor_matchings(ne_matchings, 'x', grid);

	int left_pid=-1, right_pid=-1, it, err_ct=0;
	for (it = 0; it < nProcs; it ++){
		if(ne_matchings[it] != NULL){
			err_ct++;
			list_reset_iter(ne_matchings[it]);
			surface *surf;
			if (list_length(ne_matchings[it]) > 1)
				printf("\n BADDDDDDDDDD\n");
			//while (list_has_next(ne_matchings[it])){}
			surf = list_get_next(ne_matchings[it]);
			if (surf->dir < 0)
				left_pid = surf->neighbor;
			else
				right_pid = surf->neighbor;
		}
	}
	if (err_ct > 2 )
		printf("\n err_ct is %d", err_ct);

	MPI_Request reqs[2];
	List *null_list = list_init();
	double left_pro, right_pro;
	
	while (mostwork > 1.1*avgwork){
		// If you have a bit too much, don't worry
		propensity = work - 1.05*avgwork;
		//Tell neighbors propensity
		MPI_Isend(&(propensity), 1, MPI_INT, left_pid, TAG_PROP_LEFT, MPI_COMM_WORLD, &reqs[0]);
		MPI_Isend(&(propensity), 1, MPI_INT, right_pid, TAG_PROP_RIGHT, MPI_COMM_WORLD, &reqs[1]);
		// recv same info back from them
		MPI_Recv(&(left_pro), 1, MPI_INT, left_pid, TAG_PROP_RIGHT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&(right_pro), 1, MPI_INT, right_pid, TAG_PROP_LEFT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Request_free(&reqs[0]);
		MPI_Request_free(&reqs[1]);

		if (propensity > 0){
			//If both neighbors are givers, give out
			if (left_pro > 0 && right_pro > 0){
				double center = (pymin + dy*(jmax-jmin)*.5);
				if (center < y_max/2){
					list_reset_iter(ne_matchings[left_pid]);
					give_take_surface(grid, list_get_next(ne_matchings[left_pid]), left_pid, null_list, right_pid);
				}
				else {
					list_reset_iter(ne_matchings[right_pid]);
					give_take_surface(grid, null_list, left_pid, list_get_next(ne_matchings[right_pid]), right_pid);
				}
			}
			// Else give to lowest propensity
			else if (left_pro < right_pro){
				list_reset_iter(ne_matchings[left_pid]);
				give_take_surface(grid, list_get_next(ne_matchings[left_pid]), left_pid, null_list, right_pid);
			}
			else {
				list_reset_iter(ne_matchings[right_pid]);
				give_take_surface(grid, null_list, left_pid, list_get_next(ne_matchings[right_pid]), right_pid);
			}
			//calculate propensity left and right
			//Give left or right and give null other way
			//update work based on what was given
		}
		else{
			//Give null left and right
			give_take_surface(grid, null_list, left_pid, null_list, right_pid);
		}
		//Recieve from left and right and update work
		MPI_Allreduce(&work, &mostwork, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	surface *curr;
	list_free(null_list, false);
	for (it = 0; it < nProcs; it ++){
		if(ne_matchings[it] != NULL){
			while (list_has_next(ne_matchings[it])){
				curr = list_get_next(ne_matchings[it]);
				list_free(curr->cells, false);
				list_pop(ne_matchings[it], true);
			}
			list_free(ne_matchings[it], true);
		}
	}
	
}
