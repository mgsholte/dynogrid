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
		if (curCell->owner != pid && curCell->owner != -1 && curCell->owner != -2){
			return true;
		}
	}
	return false;
}//end is_ghost()

// heavily simplified for 1D
void determine_neighbor_matchings(List* ne_matchings[], char dim, tree**** grid){
	// loop over each ghost cell.
	// loop over entire grid to find the ghost cells, this is the easiest way to find them
	// this can't be integrated with above identical loop since the push must be completed before checking to see which pushed things need to be sent to neighbors
	int i,j,k, ne_num, jghost; //ne_num is for "neighbor_number"
	
	tree *curCell = NULL;

	//instantiate all pointers to NULL for checking to see if whe need to malloc them in loop:
	for(ne_num = 0; ne_num < nProcs; ne_num++){
		ne_matchings[ne_num] = NULL;
	}//end for

	/*	determine surface matchings (i.e. list of cells to pass per neighbor)
		in the x (same as i) direction */
	int neighbor, dir = 0;
	int layer = -1;
	surface* prev_surf_match = NULL;
	
		//for (j = jmin+1; j <= jmax-2; j += (jmax-jmin-3)) {
	if(dim == 'y'){
		// if processor is 1 real cell wide, kinda hack it to work. both surfaces are the same but it can only give one way so that's no problem
		if (jmax-jmin-3 == 0) {
			surface *surf_match = malloc(sizeof(surface));
			surf_match->cells = list_init();
			surf_match->dir = 1;
			jghost = jmax-1;

			surf_match->neighbor = grid[imin+1][jghost][kmin+1]->owner;

			if (surf_match->neighbor == -1) {
				return;
			}

			ne_matchings[surf_match->neighbor] = list_init();
			list_add(ne_matchings[surf_match->neighbor], surf_match);
		
			for (i = imin+1; i < imax-1; ++i) {
				for (k = kmin+1; k < kmax-1; ++k) {
					curCell = grid[i][j][k];
					list_add(surf_match->cells, curCell);
				}
			}
			
			return;
		}		
		// first set up the surface	
		surface *surf_match = malloc(sizeof(surface));

		// do interface with left neighbor
		j = jmin+1;
		surf_match->dir = -1;
		jghost = j + surf_match->dir;

		surf_match->cells = list_init();
		surf_match->neighbor = grid[imin+1][jghost][kmin+1]->owner;

		if (surf_match->neighbor != -1) {
			ne_matchings[surf_match->neighbor] = list_init();
			list_add(ne_matchings[surf_match->neighbor], surf_match);

			
			// loop over the entire surface at this j value
			for (i = imin+1; i < imax-1; ++i) {
				for (k = kmin+1; k < kmax-1; ++k) {
					curCell = grid[i][j][k];
					list_add(surf_match->cells, curCell);
					
				}//end k loop
			}//end i loop
		}
		
		// now do the right interface
		surf_match = malloc(sizeof(surface));

		j = jmax-2;
		surf_match->dir = 1;
		jghost = j + surf_match->dir;

		surf_match->cells = list_init();
		surf_match->neighbor = grid[imax-2][jghost][kmax-2]->owner;

		if (surf_match->neighbor != -1) {
			ne_matchings[surf_match->neighbor] = list_init();
			list_add(ne_matchings[surf_match->neighbor], surf_match);

			for (i = imin+1; i < imax-1; ++i) {
				for (k = kmin+1; k < kmax-1; ++k) {
					curCell = grid[i][j][k];
					list_add(surf_match->cells, curCell);
					
				}//end k loop
			}//end i loop
		}
	}//end y dimension

}//end determine_neighbor_matchings(List** ne_matchings, char dim)

double get_work(tree ****grid){
	tree * curCell = NULL;
	int partwork=0, cellwork=0;
	int i,j,k;
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
	return .8*partwork + .2*cellwork;
}

void Balance(tree ****grid){
	double work, avgwork, mostwork, propensity;

	work = get_work(grid);

	// Find out how much work there is in total:
	MPI_Allreduce(&work, &avgwork, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	avgwork = avgwork / (double)nProcs;
	// Find out the most work that a processor in the group has:
	MPI_Allreduce(&work, &mostwork, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	
	// List (*cells_to_send)[nProcs];
	List* ne_matchings[nProcs];
	determine_neighbor_matchings(ne_matchings, 'y', grid);

	// figure out the neighbor pids
	int left_pid=-1, right_pid=-1, it, ne_ct=0;
	for (it = 0; it < nProcs; it++){
		if(ne_matchings[it] != NULL){
			ne_ct++;
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
	if (ne_ct > 2 )
		printf("\n ne_ct is %d", ne_ct);
	//printf("\n Hello, I am proc %d, and my left and right neighbors are %d and %d, errct is %d\n", pid, left_pid, right_pid, ne_ct);

	MPI_Request reqs[4];
	List *null_list = list_init();
	double left_pro, right_pro;

	while (mostwork > 1.18*avgwork){
		determine_neighbor_matchings(ne_matchings, 'y', grid);
		// If you have a bit too much, don't worry
		propensity = work - 1.05*avgwork;
		//Tell neighbors propensity
		if (left_pid != -1) {
			MPI_Isend(&(propensity), 1, MPI_DOUBLE, left_pid, TAG_PROP_LEFT, MPI_COMM_WORLD, &reqs[0]);
		} else {
			reqs[0] = MPI_REQUEST_NULL;
		}

		if (right_pid != -1) {
			MPI_Isend(&(propensity), 1, MPI_DOUBLE, right_pid, TAG_PROP_RIGHT, MPI_COMM_WORLD, &reqs[1]);
			// recv same info back from them
		} else {
			reqs[1] = MPI_REQUEST_NULL;
		}

		if (left_pid != -1) {
			MPI_Irecv(&(left_pro), 1, MPI_DOUBLE, left_pid, TAG_PROP_RIGHT, MPI_COMM_WORLD, &reqs[2]);
		} else {
			left_pro = 1.0/0.0;
			reqs[2] = MPI_REQUEST_NULL;
		}
		if (right_pid != -1) {
			MPI_Irecv(&(right_pro), 1, MPI_DOUBLE, right_pid, TAG_PROP_LEFT, MPI_COMM_WORLD, &reqs[3]);
		} else {
			right_pro = 1.0/0.0;
			reqs[3] = MPI_REQUEST_NULL;
		}

		MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

		// want to give
		if (propensity > 0){
			//If both neighbors are givers, give out
			surface *leftSfc = NULL, *rightSfc = NULL;
			list_reset_iter(ne_matchings[left_pid]);
			list_reset_iter(ne_matchings[right_pid]);
			if (right_pid != -1) {
				rightSfc = list_get_next(ne_matchings[right_pid]);
			}
			if (left_pid != -1) {
				leftSfc = list_get_next(ne_matchings[left_pid]);
			}
			if (left_pro > 0 && right_pro > 0){
				double center = (pymin + dy*(jmax-jmin)*.5);
				if (center < y_max/2){
					give_take_surface(grid, leftSfc->cells, left_pid, null_list, right_pid);
				}
				else {
					give_take_surface(grid, null_list, left_pid, rightSfc->cells, right_pid);
				}
			}
			// Else give to lowest propensity
			else if (left_pro < right_pro){
				give_take_surface(grid, leftSfc->cells, left_pid, null_list, right_pid);
			}
			else {
				give_take_surface(grid, null_list, left_pid, rightSfc->cells, right_pid);
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
		work = get_work(grid);
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
	
	// do presend to communicate message sizes to expect
	List *move_list;
	int neighbor_id;
	int i,j,k;
	int l,m,n;
	vec3 pos;
	MPI_Request presend_reqs[4*nNeighbors];
	int nTreeSends[nNeighbors], nParticleSends[nNeighbors];
	int nTreeRecvs[nNeighbors], nParticleRecvs[nNeighbors];
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		if (i == 0) { 
			move_list = list_u;
			neighbor_id = id_u; 
		} else if (i == 1) {
			move_list = list_d;
			neighbor_id = id_d; 
		}
		if (neighbor_id == -1){
			presend_reqs[2*i] = MPI_REQUEST_NULL;
			presend_reqs[2*i+1] = MPI_REQUEST_NULL;
			presend_reqs[2*(i+nNeighbors)] = MPI_REQUEST_NULL;
			presend_reqs[2*(i+nNeighbors)+1] = MPI_REQUEST_NULL;
			continue;
		}
		// send move_list, cells within will be converted into these buffers before send
		nTreeSends[i] = list_length(move_list);
		nParticleSends[i] = 0;
		list_reset_iter(move_list);
		while(list_has_next(move_list)) {
			List *parts = ((tree*) list_get_next(move_list))->particles;
			nParticleSends[i] += list_length(parts);
		}
		MPI_Request *tmp = mpi_tree_presend(&nTreeSends[i], &nParticleSends[i], neighbor_id);
		presend_reqs[2*i] = tmp[0];
		presend_reqs[2*i+1] = tmp[1];
		free(tmp);

		tmp = mpi_tree_prerecv(&nTreeRecvs[i], &nParticleRecvs[i], neighbor_id);
		presend_reqs[2*(i+nNeighbors)] = tmp[0];
		presend_reqs[2*(i+nNeighbors)+1] = tmp[1];
	}

	// send list_u up then list_d down (up to 1 can be a real send)
	MPI_Waitall(4*nNeighbors, presend_reqs, MPI_STATUSES_IGNORE);

	MPI_Request all_reqs[6*nNeighbors];
	
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { move_list = list_u; neighbor_id = id_u; send_dir6[i] = 'u'; }
		else if (i == 1) { move_list = list_d; neighbor_id = id_d; send_dir6[i] = 'd'; }
		if (neighbor_id == -1){
			all_reqs[3*i] = MPI_REQUEST_NULL;
			all_reqs[3*i+1] = MPI_REQUEST_NULL;
			all_reqs[3*i+2] = MPI_REQUEST_NULL;
			buff_send_trees[i] = malloc(0);
			buff_send_parts[i] = malloc(0);
			buff_send_part_list_lengths[i] = malloc(0);
			continue;
		}
		
		// fill neighbor_owners
		tree* move_tree;
		list_reset_iter(move_list);
		while (list_has_next(move_list)) {
			move_tree = (tree*) list_get_next(move_list);
			pos = move_tree->loc;
			l = imin + round_i((pos.x - pxmin)/dx);
			m = jmin + round_i((pos.y - pymin)/dy);
			n = kmin + round_i((pos.z - pzmin)/dz);
			if (send_dir6[i] == 'u') {
				for (j = -1; j <= 1; ++j) {
					for (k = -1; k <= 1; ++k) {
						move_tree->neighbor_owners[j+1][k+1] = base_grid[l+k][m+1][n+j]->owner; //yes j and k are used correctly here
					}
				}
			} else if (send_dir6[i] == 'd') {
				for (j = -1; j <= 1; ++j) {
					for (k = -1; k <= 1; ++k) {
						move_tree->neighbor_owners[j+1][k+1] = base_grid[l+k][m-1][n+j]->owner; //yes j and k are used correctly here
					}
				}
			}		
		}
		
		buff_send_trees[i] = malloc(nTreeSends[i] * sizeof(simple_tree));
		buff_send_parts[i] = malloc(nParticleSends[i] * sizeof(particle));
		buff_send_part_list_lengths[i] = malloc(nTreeSends[i] * sizeof(simple_tree));

		// send move_list, cells within will be converted into these buffers before send
		trees_parts_and_lengths = mpi_tree_send(move_list, neighbor_id, buff_send_trees[i], buff_send_parts[i], buff_send_part_list_lengths[i], &send_dir6[i]);
		all_reqs[3*i] = trees_parts_and_lengths[0];
		all_reqs[3*i+1] = trees_parts_and_lengths[1];
		all_reqs[3*i+2] = trees_parts_and_lengths[2];
		free(trees_parts_and_lengths);
		
		// clean out trees being given away by turning them into ghosts, their neighbor ghosts into NULLs or init ghosts, and their neighbor init ghosts into NULLs
		// only works for 1D, whole surfaces at a time
		list_reset_iter(move_list);
		while (list_has_next(move_list)) {
			pos = ((tree*) list_get_next(move_list))->loc;
			l = imin + round_i((pos.x - pxmin)/dx);
			m = jmin + round_i((pos.y - pymin)/dy);
			n = kmin + round_i((pos.z - pzmin)/dz);
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
			++jmin;
			pymin += dy;
		} else if (send_dir6[i] == 'd' && move_list->length != 0) {
			--jmax;
		}
			
	}
	
		
	// buffer receive arrays
	simple_tree* buff_recv_trees[nNeighbors];
	particle* buff_recv_parts[nNeighbors];
	int* buff_recv_part_list_lengths[nNeighbors];
	char recv_dir6[nNeighbors];
	
	// receive from neighbors: either buffers or nothing
	// receive first from up then from down (up to 2 can be real receives)
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { neighbor_id = id_u; }
		else if (i == 1) { neighbor_id = id_d; }
		
		if (neighbor_id == -1){
			all_reqs[3*(i+nNeighbors)] = MPI_REQUEST_NULL;
			all_reqs[3*(i+nNeighbors)+1] = MPI_REQUEST_NULL;
			all_reqs[3*(i+nNeighbors)+2] = MPI_REQUEST_NULL;
			buff_recv_trees[i] = malloc(0);
			buff_recv_parts[i] = malloc(0);
			buff_recv_part_list_lengths[i] = malloc(0);
			continue;
		}

		buff_recv_trees[i] = malloc(nTreeRecvs[i] * sizeof(simple_tree));
		buff_recv_parts[i] = malloc(nParticleRecvs[i] * sizeof(particle));
		buff_recv_part_list_lengths[i] = malloc(nTreeRecvs[i] * sizeof(simple_tree));

		trees_parts_and_lengths = mpi_tree_recv(neighbor_id, buff_recv_trees[i], nTreeRecvs[i], buff_recv_parts[i], nParticleRecvs[i], buff_recv_part_list_lengths[i], &recv_dir6[i]);
		all_reqs[3*(i+nNeighbors)] = trees_parts_and_lengths[0];
		all_reqs[3*(i+nNeighbors)+1] = trees_parts_and_lengths[1];
		all_reqs[3*(i+nNeighbors)+2] = trees_parts_and_lengths[2];
		free(trees_parts_and_lengths);
	}	
	
	// wait for the sends and receives
	MPI_Waitall(6*nNeighbors, all_reqs, MPI_STATUSES_IGNORE);

	// once receives are done, can start unpacking buffers, putting new trees where they belong, while adjusting base_grid as necessary
	// skips empty lists, but always expects a list from each neighbor in the current direction of passing
	List* new_trees;
	tree* new_tree;
	double pymax = pymin+dy*(jmax-jmin); //pxmin is start of ghosts using global x position
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { neighbor_id = id_u; }
		else if (i == 1) { neighbor_id = id_d; }
		
		if (neighbor_id == -1) {
			free(buff_recv_trees[i]);
			free(buff_recv_parts[i]);
			free(buff_recv_part_list_lengths[i]);
			continue;
		}
		//NB: frees the arrays after it unpacks
		new_trees = mpi_tree_unpack(buff_recv_trees[i], buff_recv_parts[i], buff_recv_part_list_lengths[i], nTreeRecvs[i]);
		list_reset_iter(new_trees);
		
		// first new_tree is used to check if base_grid is big enough, we loop through the rest later
		if (list_has_next(new_trees)) {
			new_tree = (tree*) list_get_next(new_trees);
		} else {
			continue; //no new_tree's to see here
		}
		
		// possibly change imin/imax plus possibly resize whole base_grid
		// new_tree only provides spatial to work with, so compare that to pxmin
		if ((new_tree->loc.y-dy/2) < pymin) {	//all new_tree's have same x. using dx/2 to prevent rounding issues
			if (jmin == 0) {
				resize_allocation(base_grid); //resets wi to 2*(imax-imin) and centers around (imin+imax)/2
			} else if (jmin < 0) {
				printf("ERROR: unexpected jmin < 0");
				MPI_Finalize();
				return;
			}
			--jmin;
			pymin -= dy;
		} else if((new_tree->loc.y+3*dy/2) > pymax) {
			if (jmax-1 == wj-1) {
				resize_allocation(base_grid);
			} else if (jmax-1 > wj-1) {
				printf("ERROR: unexpected jmax-1 > wj-1");
				MPI_Finalize();
				return;
			}
			++jmax;
		}
		
		// we're good now: let's insert our new_tree's
		convert_ghost2real_and_reghost(base_grid, new_tree, recv_dir6[i]);
		while (list_has_next(new_trees)) {
			new_tree = (tree*) list_get_next(new_trees);
			convert_ghost2real_and_reghost(base_grid, new_tree, recv_dir6[i]);
		}
		
	}
	
	// free buffers etc.
	for (i = 0; i < nNeighbors; ++i) {
		free(buff_send_trees[i]);
		free(buff_send_parts[i]);
		free(buff_send_part_list_lengths[i]);

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
	int i = imin + round_i((new_tree->loc.x - pxmin)/dx); //round to force correct int value
	int j = jmin + round_i((new_tree->loc.y - pymin)/dy);
	int k = kmin + round_i((new_tree->loc.z - pzmin)/dz);
	
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
				tree_init(get_loc(i+di,j+dj,k+dk), new_tree->neighbor_owners[1+*d1][1+*d2]);
				
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
	// dont want to free the particle list we copied from this
	new_tree->particles = list_init();
	tree_free(new_tree);
} //convert_ghost2real_and_reghost
