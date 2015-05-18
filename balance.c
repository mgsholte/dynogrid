#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "balance.h"
#include "grid.h"
#include "dynamics.h"
#include "mpi_dyno.h"


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
	
	if(dim == 'y'){
		surface *surf_match = NULL;
		for (i = imin+1; i < imax-1; ++i) {
			for (j = jmin; j < jmax; ++j) {
				for (k = kmin+1; k < kmax-1; ++k) {
					curCell = grid[i][j][k];
					if(is_ghost(curCell)) {
						//determine the direction of the neighbor (for this particular cell!!!)
						if(j == 0){
							/* 	we can't access the j-1 element of the grid, but we know that
								the neighbor must be to our LEFT direction (hence the dir = "-1") */
							dir = -1; //dir for "direction"
							nextCell = grid[i][j+1][k];
						}else if(j == wj-1){
							/* 	we can't access the j+1 element of the grid, but we know that
								the neighbor must be to our RIGHT direction (hence the dir ="+1") */
							dir = 1;
							prevCell = grid[i][j-1][k];
						}else{
							//we know we can access the j-1 and j+1 element of the grid
							prevCell = grid[i][j-1][k];
							nextCell = grid[i][j+1][k];
						
							if(is_real(nextCell)){
								dir = -1;
							}else if(is_real(prevCell)){
								dir = 1;
							}//end if else if
						}//end if else for checking extremes
						//set the neighbor pid:
						neighbor = curCell->owner;
						//set the layer value for the surface struct (i.e. )
						layer = j;
						
						if( prev_surf_match != NULL &&
							prev_surf_match->neighbor == neighbor &&
							prev_surf_match->layer == layer &&
							prev_surf_match->dir == dir){
								//then we are dealing with the same surface matching
								surf_match = prev_surf_match;
						}else{
							surf_match = malloc(sizeof(*surf_match));

							surf_match->cells = list_init();
							surf_match->neighbor = curCell->owner;
							surf_match->dir = dir;
							surf_match->layer = layer;

							if((ne_matchings[surf_match->neighbor]) == NULL){
								//then the List has not yet been initialized so we have to malloc it and initialize it:
								ne_matchings[surf_match->neighbor] = list_init();
							}
							list_add(ne_matchings[surf_match->neighbor], surf_match);
						}//end if else
						
						if(dir < 0){
							list_add(surf_match->cells, nextCell);
						}else if(dir > 0){
							list_add(surf_match->cells, prevCell);
						}
						/* reminder: 
								ne_matchings ......	an array of
													List pointers
									where each list holds "surface"s (a new type of struct...see Balance.h)
									each "surface" has meta data and
									a List of trees */
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
			//printf("After while length of ne_matchings[%d] is %d\ncount is %d \n", ne_num, list_length(ne_matchings[ne_num]), ct);
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
					cellwork += curCell->nPoints;
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
	determine_neighbor_matchings(ne_matchings, 'y', grid);

	int left_pid=-1, right_pid=-1, it, err_ct=0;
	for (it = 0; it < nProcs; it ++){
		if(ne_matchings[it] != NULL){
			err_ct++;
			list_reset_iter(ne_matchings[it]);
			surface *surf;
			if (list_length(ne_matchings[it]) > 1)
				printf("\n BADDDDDDDDDD length is %d\n", list_length(ne_matchings[it]));
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
	//printf("\n Hello, I am proc %d, and my left and right neighbors are %d and %d, errct is %d\n", pid, left_pid, right_pid, err_ct);

	MPI_Request reqs[2];
	List *null_list = list_init();
	double left_pro, right_pro;

	while (mostwork > 1.1*avgwork){
		//TODO: get rid of this.
		work = avgwork;
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
		//Need to actually update work. Maybe I can convince Max to do this so I don't have to loop over all cells
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

// determine which trees to send where, and send them. then take trees given to you, adjust your base_grid if necessary, and insert the trees
void give_take_surface(tree**** base_grid, List* list_u, int id_u, List* list_d, int id_d) {
	//direction-dependent pointers to variables
	
	int nNeighbors = 2; // for 1D balancing
	
	/* // now getting list_u and list_d passed in, so creating moving_trees is unnecessary
	
	int neighbors[numNeighbors] = proc neighbors
	
	// create moving_trees: an array of List*s pointing to tree*s. The use of tree*s is not made explicit here, may need to cast later
	// always List*, never just List, for C reasons
	List* moving_trees[nNeighbors]
	for (neighbors)
		moving_trees[neighbor_id] = list_init();
	end
	
	// fill moving_trees
	int neighbor_id
	for (tree*s at position i = imax-2)
		if (tree->owner == pid) //not a ghost
			neighbor_id = base_grid[i+1][j][k]->owner
			list_add(moving_trees[neighbor_id], tree*)
		end
	end
	*/
	
	// buffer send arrays. both get malloced in mpi_tree_send, the ** is so after the * gets malloced and therefore changes, I still have the pointer to it
	simple_tree *buff_send_trees[nNeighbors];
	particle *buff_send_parts[nNeighbors];
	int *buff_send_part_list_lengths[nNeighbors];
	char send_dir6[nNeighbors]; // the direction of the incoming surface, pointing from giver to taker

	// need to track MPI_request objects for a waitall later, in order to be non-blocking
	MPI_Request req_send_trees[nNeighbors];
	MPI_Request req_send_parts[nNeighbors];
	MPI_Request req_send_lengths[nNeighbors];
	MPI_Request *trees_parts_and_lengths;
	
	// send first up then down (up to 1 can be a real send)
	List *move_list;
	int neighbor_id;
	int i,j,k;
	int l,m,n;
	vec3 pos;
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { move_list = list_u; neighbor_id = id_u; send_dir6[i] = 'u'; }
		else if (i == 1) { move_list = list_d; neighbor_id = id_d; send_dir6[i] = 'd'; }
		if (neighbor_id == -1){
			req_send_trees[i] = MPI_REQUEST_NULL;
			req_send_parts[i] = MPI_REQUEST_NULL;
			req_send_lengths[i] = MPI_REQUEST_NULL;
			buff_send_trees[i] = malloc(0);
			buff_send_parts[i] = malloc(0);
			buff_send_part_list_lengths[i] = malloc(0);
			continue;
		}

		// send move_list, cells within will be converted into these buffers before send
		trees_parts_and_lengths = mpi_tree_send(move_list, neighbor_id, &buff_send_trees[i], &buff_send_parts[i], &buff_send_part_list_lengths[i], &send_dir6[i]);
		req_send_trees[i] = trees_parts_and_lengths[0];
		req_send_parts[i] = trees_parts_and_lengths[1];
		req_send_lengths[i] = trees_parts_and_lengths[2];
		free(trees_parts_and_lengths);
		
		// clean out trees being given away by turning them into ghosts, their neighbor ghosts into NULLs or init ghosts, and their neighbor init ghosts into NULLs
		// only works for 1D, whole surfaces at a time
		list_reset_iter(move_list);
		while (list_has_next(move_list)) {
			pos = ((tree*) list_get_next(move_list))->loc;
			l = imin + (int) round((pos.x - pxmin)/dx);
			m = jmin + (int) round((pos.y - pymin)/dy);
			n = kmin + (int) round((pos.z - pzmin)/dz);
			base_grid[l][m][n]->owner = neighbor_id;
			if (send_dir6[i] == 'u') {
				for (j = -1; j <= 1; ++j) {
					for (k = -1; k <= 1; ++k) {
						if ( base_grid[l+j][m-1][n+k] != NULL )
							tree_free( base_grid[l+j][m-1][n+k] );
						base_grid[l+j][m-1][n+k] = NULL;
					}
				}
			} else if (send_dir6[i] == 'd') {
				for (j = -1; j <= 1; ++j) {
					for (k = -1; k <= 1; ++k) {
						if ( base_grid[l+j][m+2][n+k] != NULL )
							tree_free( base_grid[l+j][m+2][n+k] );
						base_grid[l+j][m+2][n+k] = NULL;
						base_grid[l+j][m+1][n+k]->owner = -2;
					}
				}
			}
		}
		
		// if giving up a surface, change imin and pxmin or imax
		if (send_dir6[i] == 'u' && move_list->length != 0) {
			++imin;
			pxmin += dx;
		} else if (send_dir6[i] == 'd' && move_list->length != 0) {
			--imax;
		}
			
	}
	
		
	// buffer receive arrays
	simple_tree* buff_recv_trees[nNeighbors];
	particle* buff_recv_parts[nNeighbors];
	int* buff_recv_part_list_lengths[nNeighbors];
	char recv_dir6[nNeighbors];
	int lengths_of_trees[nNeighbors]; //int for the tree array, pointer so it can be filled in mpi_tree_recv
	
	MPI_Request req_recv_trees[nNeighbors];
	MPI_Request req_recv_parts[nNeighbors];
	MPI_Request req_recv_lengths[nNeighbors];
	
	// receive from neighbors: either buffers or nothing
	// receive first from up then from down (up to 2 can be real receives)
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { neighbor_id = id_u; }
		else if (i == 1) { neighbor_id = id_d; }
		
		if (neighbor_id == -1){
			req_recv_trees[i] = MPI_REQUEST_NULL;
			req_recv_parts[i] = MPI_REQUEST_NULL;
			req_recv_lengths[i] = MPI_REQUEST_NULL;
			buff_recv_trees[i] = malloc(0);
			buff_recv_parts[i] = malloc(0);
			buff_recv_part_list_lengths[i] = malloc(0);
			continue;
		}

		trees_parts_and_lengths = mpi_tree_recv(neighbor_id, &buff_recv_trees[i], &buff_recv_parts[i], &buff_recv_part_list_lengths[i], &lengths_of_trees[i], &recv_dir6[i]);
		req_recv_trees[i] = trees_parts_and_lengths[0];
		req_recv_parts[i] = trees_parts_and_lengths[1];
		req_recv_lengths[i] = trees_parts_and_lengths[2];
		free(trees_parts_and_lengths);
	}	
	
	// wait for receives
	MPI_Waitall(nNeighbors, req_recv_trees, MPI_STATUSES_IGNORE);
	MPI_Waitall(nNeighbors, req_recv_parts, MPI_STATUSES_IGNORE);
	MPI_Waitall(nNeighbors, req_recv_lengths, MPI_STATUSES_IGNORE);
	
	// once receives are done, can start unpacking buffers, putting new trees where they belong, while adjusting base_grid as necessary
	// skips empty lists, but always expects a list from each neighbor in the current direction of passing
	List* new_trees;
	tree* new_tree;
	double pxmax = pxmin+dx*(imax-imin); //pxmin is start of ghosts using global x position
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { neighbor_id = id_u; }
		else if (i == 1) { neighbor_id = id_d; }
		
		new_trees = mpi_tree_unpack(&buff_recv_trees[i], &buff_recv_parts[i], &buff_recv_part_list_lengths[i], &lengths_of_trees[i]);
		list_reset_iter(new_trees);
		
		// first new_tree is used to check if base_grid is big enough, we loop through the rest later
		if (list_has_next(new_trees)) {
			new_tree = (tree*) list_get_next(new_trees);
		} else {
			continue; //no new_tree's to see here
		}
		
		// possibly change imin/imax plus possibly resize whole base_grid
		// new_tree only provides spatial to work with, so compare that to pxmin
		if ((new_tree->loc.x+dx/2) < pxmin) {	//all new_tree's have same x. using dx/2 to prevent rounding issues
			if (imin == 0) {
				resize_allocation(base_grid); //resets wi to 2*(imax-imin) and centers around (imin+imax)/2
			} else if (imin < 0) {
				printf("ERROR: unexpected imin < 0");
				MPI_Finalize();
				return;
			}
			--imin;
			pxmin -= dx;
		} else if((new_tree->loc.x-dx/2) > pxmax) {
			if (imax == wi-1) {
				resize_allocation(base_grid);
			} else if (imax > wi-1) {
				printf("ERROR: unexpected imax > wi-1");
				MPI_Finalize();
				return;
			}
			++imax;
		}
		
		// we're good now: let's insert our new_tree's
		convert_ghost2real_and_reghost(base_grid, new_tree, recv_dir6[i]);
		while (list_has_next(new_trees)) {
			new_tree = (tree*) list_get_next(new_trees);
			convert_ghost2real_and_reghost(base_grid, new_tree, recv_dir6[i]);
		}
		
	}
	
	// wait for sends
	MPI_Waitall(nNeighbors, req_send_trees, MPI_STATUSES_IGNORE);
	MPI_Waitall(nNeighbors, req_send_parts, MPI_STATUSES_IGNORE);
	MPI_Waitall(nNeighbors, req_send_lengths, MPI_STATUSES_IGNORE);
	
	// free buffers etc.
	for (i = 0; i < nNeighbors; ++i) {
		free(buff_send_trees[i]);
		free(buff_send_parts[i]);
		free(buff_send_part_list_lengths[i]);

		free(buff_recv_trees[i]);
		free(buff_recv_parts[i]);
		free(buff_recv_part_list_lengths[i]);
	}

	
} //end give_take_surface


// re-allocates base_grid such that wi is 2*(imax-imin) and grid is centered around (imin+imax)/2 (so 50% padding each side)
void resize_allocation(tree**** base_grid) {
	wi = 2*(imax-imin);
	wj = 2*(jmax-jmin);
	wk = 2*(kmax-kmin);
	int di = wi/4-imin; //how far imin is being shifted to re-center
	imin = wi/4;
	imax = imin+wi/2;
	int dj = wj/4-jmin;
	jmin = wj/4;
	jmax = jmin+wj/2;
	int dk = wk/4-kmin;
	kmin = wk/4;
	kmax = kmin+wk/2;
	
	tree ****new_grid = (tree****) malloc( wi * sizeof(tree***) ); // allocate an array of pointers to rows-depthwise
	int i, j, k, n;
	for (i = 0; i < wi; ++i) {
		new_grid[i] = (tree***) malloc( wj * sizeof(tree**) );  // allocate the row
		for (j = 0; j < wj; ++j) {
			new_grid[i][j] = (tree**) malloc( wk * sizeof(tree*) );  // allocate the row
			for (k = 0; k < wk; ++k) {
				// include init ghosts
				if (i >= imin && i <= imax &&
					j >= jmin && j <= jmax &&
					k >= kmin && k <= kmax) {
					
					new_grid[i][j][k] = base_grid[i-di][j-dj][k-dk]; //translating old tree*s onto new grid, includes some NULLs
					
				} else {
					new_grid[i][j][k] = NULL;
				}
			}
		}
	}
	
	for (i = 0; i < wi; ++i) {
		for (j = 0; j < wj; ++j) {
			free(base_grid[i][j]);
		}
		free(base_grid[i]);
	}
	free(base_grid);
	
	base_grid = new_grid;
} //end resize_allocation

// dir6 is direction that tree was passed to get here, can be 'l','r','u','d','f','b' (for 1D just 'l' or 'r')
// init ghosts should be overwritten by new ghosts and recreated as necessary
//  - when moving left, up, or forward, up to 23 new init ghosts are needed; otherwise up to 14 new init ghosts are needed
void convert_ghost2real_and_reghost(tree**** base_grid, tree* new_tree, char dir6) {
	int i = imin + 1 + (int) round((new_tree->loc.x - pxmin)/dx); //round to force correct int value
	int j = jmin + 1 + (int) round((new_tree->loc.y - pymin)/dy);
	int k = kmin + 1 + (int) round((new_tree->loc.z - pzmin)/dz);
	
	// for new tree, need to: give particles, set owner
	// we leave children and point values alone, those should be fine
	base_grid[i][j][k]->particles = new_tree->particles;
	base_grid[i][j][k]->owner = pid;
	
	// use pointers to iterate over the 9 (or fewer) new ghost cells, respecting which plane they lie in
	int di,dj,dk;
	int *d1, *d2;
	if (dir6 == 'l') {
		di = 1;
		d1 = &dj;
		d2 = &dk;
	} else if (dir6 == 'r') {
		di = -1;
		d1 = &dj;
		d2 = &dk;
	} else if (dir6 == 'u') {
		d2 = &di;
		dj = +1;
		d1 = &dk;
	} else if (dir6 == 'd') {
		d2 = &di;
		dj = -1;
		d1 = &dk;
	} else if (dir6 == 'f') {
		d1 = &di;
		d2 = &dj;
		dk = +1;
	} else if (dir6 == 'b') {
		d1 = &di;
		d2 = &dj;
		dk = -1;
	}
	
	// create new ghosts
	for (*d1 = -1; *d1 <= 1; ++*d1) {
		for (*d2 = -1; *d2 <= 1; ++*d2) {
			if (base_grid[i+di][j+dj][k+dk] == NULL) {
				
				// tree_init, includes setting points[0] and owner
				base_grid[i+di][j+dj][k+dk] = 
				tree_init(
				get_loc(i+di,j+dj,k+dk), 
				new_tree->
				neighbor_owners[1+*d1][1+*d2]);
				
			// if init ghost then convert into ghost instead of re-mallocing
			} else if (base_grid[i+di][j+dj][k+dk]->owner == -2) {
			
				// just change owner. will add points and call laser in next for loop
				base_grid[i+di][j+dj][k+dk]->owner = new_tree->neighbor_owners[1+*d1][1+*d2];
				
			}
			// if neither NULL nor init ghost, this cell is already real or ghost and should be left alone
		}
	}
	// create new init ghosts. then make other 7 points point to neighbor points, but only for ghosts made in previous loop
	// when dir6 is left, up, or forward, you're moving into init ghosts and so can use pre-existing ones in your layer. however you need to create a new layer of them.
	// when dir6 is right, down, or back, there are no init ghosts where you're going so you have to make a few new ones in your layer. however there is no new layer to make.
	// both cases are covered here by only putting new init ghosts to the right, down, and back of ghosts created in previous loop
	int n;
	for (*d1 = -1; *d1 <= 1; ++*d1) {
		for (*d2 = -1; *d2 <= 1; ++*d2) {
			// "if tree is a ghost created in the previous for loop"
			// recently made ghosts can be identified by points[1] == NULL (which is specifically set in tree_init)
			// second condition is to avoid choosing init ghosts as well
			if (base_grid[i+di][j+dj][k+dk]->root->points[1] == NULL && base_grid[i+di][j+dj][k+dk]->owner != -2) {
			
				// "for neighbors to the right, down, and back"
				// create new init ghosts as necessary. n should be thought of as binary for getX etc. (n for "neighbors")
				for (n = 1; n < 8; ++n) {
					if (base_grid[i+di+getX(n)][j+dj+getY(n)][k+dk+getZ(n)] == NULL) {
				
						// tree_init, includes setting points[0] and owner = -2
						base_grid[i+di+getX(n)][j+dj+getY(n)][k+dk+getZ(n)] = tree_init(get_loc(i+di+getX(n),j+dj+getY(n),k+dk+getZ(n)), -2);
					}
				}
				
				// make other 7 points point to neighbor points
				for (n = 1; n < 8; ++n) {
					base_grid[i+di][j+dj][k+dk]->root->points[n] = base_grid[i+di+getX(n)][j+dj+getY(n)][k+dk+getZ(n)]->root->points[0];
				}
				
				// laser
				tree_apply_fcn(base_grid[i+di][j+dj][k+dk], &laser);
			}
		}
	}
} //convert_ghost2real_and_reghost
