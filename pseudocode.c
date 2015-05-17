// All Pseudocode in here

// GRID CELL PASSER
/*
- load balancer just decides 'along the x axis, I'm sending [left or right], or not at all', then decides similarly for y and for z
- then an algorithm figures out which procs this means we're sending to, and constructs lists of the trees to be sent (mpi_tree_send expects lists)
	- for this we look at all non-ghost cells at position imax-2 (1 from edge), since that's what we pass
	- each looks at the owner of the ghost next to it (1 away in chosen x direction) and adds itself to a List for that proc
- communication
	- send each neighbor a List* (mpi_tree_send makes this an array), and make buffers and MPI_request objects
	- receive from each neighbor, and make buffers and MPI_request objects
	- wait all for these MPI_request objects
	- unpack the received buffers (simple_tree*, particle*, and int* become just tree**, an array of tree pointers)
	- assign those tree*s where they belong in their new home based on their loc's
*/

#include "grid.h"


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
	simple_tree* buff_send_trees[nNeighbors];
	particle* buff_send_parts[nNeighbors];
	int* buff_send_part_list_lengths[nNeighbors];
	char send_dir6[nNeighbors]; // the direction of the incoming surface, pointing from giver to taker

	// need to track MPI_request objects for a waitall later, in order to be non-blocking
	MPI_Request req_send_trees[nNeighbors];
	MPI_Request req_send_parts[nNeighbors];
	MPI_Request req_send_lengths[nNeighbors];
	MPI_Request* trees_parts_and_lengths;
	
	// send first up then down (up to 1 can be a real send)
	List* move_list;
	int neighbor_id;
	int i,j,k;
	vec3 pos;
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { move_list = list_u; neighbor_id = id_u; *send_dir6[i] = 'u'; }
		else if (i == 1) { move_list = list_d; neighbor_id = id_d; *send_dir6[i] = 'd' }
		
		if (neighbor_id == -1) continue;

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
			base_grid[pos.x][pox.y][pox.z]->owner = neighbor_id;
			if (send_dir[i] == 'u') {
				for (j = -1; j <= 1; ++j) {
					for (k = -1; k <= 1; ++k) {
						if ( base_grid[pos.x+j][pos.y-1][pos.z+k] != NULL )
							tree_free( base_grid[pos.x+j][pos.y-1][pos.z+k] );
						base_grid[pos.x+j][pos.y-1][pos.z+k] = NULL;
					}
				}
			} else if (send_dir[i] == 'd') {
				for (j = -1; j <= 1; ++j) {
					for (k = -1; k <= 1; ++k) {
						if ( base_grid[pos.x+j][pos.y+2][pos.z+k] != NULL )
							tree_free( base_grid[pos.x+j][pos.y+2][pos.z+k] );
						base_grid[pos.x+j][pos.y+2][pos.z+k] = NULL;
						base_grid[pos.x+j][pos.y+1][pos.z+k]->owner = -2;
					}
				}
			}
		}
		
		// if giving up a surface, change imin and pxmin or imax
		if (send_dir[i] == 'u' && move_list->length != 0) {
			++imin;
			pxmin += dx;
		} else if (send_dir[i] == 'd' && move_list->length != 0) {
			--imax;
		}
			
	}
	
		
	// buffer receive arrays
	simple_tree* buff_recv_trees[nNeighbors];
	particle* buff_recv_parts[nNeighbors];
	int* buff_recv_part_list_lengths[nNeighbors];
	char recv_dir6[nNeighbors];
	int length_of_trees[nNeighbors]; //int for the tree array, pointer so it can be filled in mpi_tree_recv
	
	MPI_Request req_recv_trees[nNeighbors];
	MPI_Request req_recv_parts[nNeighbors];
	MPI_Request req_recv_lengths[nNeighbors];
	
	// receive from neighbors: either buffers or nothing
	// receive first from up then from down (up to 2 can be real receives)
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { neighbor_id = id_u; }
		else if (i == 1) { neighbor_id = id_d; }
		
		if (neighbor_id == -1) continue;

		trees_parts_and_lengths = mpi_tree_recv(neighbor_id, &buff_recv_trees[i], &buff_recv_parts[i], &buff_recv_part_list_lengths[i], &lengths_of_trees[i], &recv_dir6[i]);
		req_recv_trees[i] = trees_parts_and_lengths[0];
		req_recv_parts[i] = trees_parts_and_lengths[1];
		req_recv_lengths[i] = trees_parts_and_lengths[2];
		free(trees_parts_and_lengths);
	}	
	
	// wait for receives
	MPI_Waitall(nNeighbors, req_recv_trees);
	MPI_Waitall(nNeighbors, req_recv_parts);
	MPI_Waitall(nNeighbors, req_recv_lengths);
	
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
		list_deset_iter(new_trees);
		
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
				return -1;
			}
			--imin;
			pxmin -= dx;
		} else if((new_tree->loc.x-dx/2) > pxmax) {
			if (imax == wi-1) {
				resize_allocation(base_grid);
			} else if (imax > wi-1) {
				printf("ERROR: unexpected imax > wi-1");
				MPI_Finalize();
				return -1;
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
	MPI_Waitall(nNeighbors, req_send_trees);
	MPI_Waitall(nNeighbors, req_send_parts);
	MPI_Waitall(nNeighbors, req_send_lengths);
	
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
	int* d1,d2;
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
				base_grid[i+di][j+dj][k+dk] = tree_init(get_loc(i+di,j+dj,k+dk), new_tree->neighbor_owners[1+*d1][1+*d2]);
				
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

------------------------------------------------------------------------
OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD  OLD
------------------------------------------------------------------------

/*
give_cells
	//excerpt
	for (takers)
		taker_ids = {}
		for (taker's ghost cells)
			if (taker_cell->owner == pid && pass_ids.contains(taker_cell->id)) // checking for right owner prob not necessary
				// taker changes
				pass particles from giver_cell to taker_cell // could use morton id to easily find giver cell; should this be per-cell or per-taker-proc?
				taker_ids.add(taker_cell->id)
				
				// giver changes
				giver_cell->owner = taker
				--xmax;
				// how do we deal with now-outdated ghosts? kill them or ignore them? killing means any ghosts at xmax+1 = NULL
			end
		end
		send_message(taker, taker_ids)
	end
*/

update_grid
	foreach (grid_cell)
		update_grid_cell(grid_cell)
	end
	
update_grid_cell(grid_cell)
	if (grid_cell.children != NULL)
		bool want_coarsen = ask_to_coarsen(grid_cell)
		if (want_coarsen)
			coarsen(grid_cell)
		else
			foreach (child)
				update_grid_cell(child)
			end
		end
	else
		bool refined = ask_to_refine(grid_cell)
		(unfinished)