// All Pseudocode in here

// GRID CELL PASSER
/*
some logic
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

/*
// TODO: This is not up to date with how load balancer logic will work
// Correct logic: for each neighbor, say whether you're giving left (into them), right (into them), or neither (possibly receiving), and by how many "faces" if you choose one of the first two
trade_cells(base_grid)
	// dirs (directions) is an array of chars, each is a char dir6 (direction w/ 6 possibilities)
	// each dir6 can be:
	//  'l' = left = x decreasing
	//  'r' = right = x increasing
	//  'u' = up = y decreasing
	//  'd' = down = y increasing
	//  'f' = forward = z decreasing
	//  'b' = back = z increasing
	char[3] dirs = load_balancer();
	for each (dir6 in dirs)
		give_and_take_cells(base_grid, dir6)
	end
end trade_cells
*/

// TODO: written only for x direction, need to generalize to y and z
// determine which trees to send where, and send them. then take trees given to you, adjust your base_grid if necessary, and insert the trees
give_take_surface(tree**** base_grid, List* list_l, int id_l, List* list_r, int id_r) {
	//direction-dependent pointers to variables
	
	int nNeighbors = 2; // for 1D balancing
	//int nNeighbors = sizeof(moving_trees); // should have been statically declared to be size of neighbors so this should work
	
	/* // now getting moving_trees passed in, so this is unnecessary
	
	int neighbors[numNeighbors] = proc neighbors
	
	// create moving_trees: an array of List*s pointing to tree*s. The use of tree*s is not made explicit here, may need to cast later
	// always List*, never just List, for C reasons
	List* moving_trees[nNeighbors]
	for (neighbors)
		*moving_trees[neighbor_id] = list_init();
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
	
	// buffer send arrays. both get malloced in mpi_tree_send
	simple_tree* buff_send_trees[nNeighbors];
	particle* buff_send_parts[nNeighbors];
	int* buff_send_part_list_lengths[nNeighbors];
	int* lengths_of_buffs[nNeighbors][3]; //one length for each of the above 3 buffers

	// TODO: add sending/receiving a char. send a &char, recv a &char. this is the dir6 of the incoming processor
	// need to track MPI_request objects for a waitall later, in order to be non-blocking
	MPI_Request req_send_trees[nNeighbors];
	MPI_Request req_send_parts[nNeighbors];
	MPI_Request req_send_lengths[nNeighbors];
	MPI_Request trees_parts_and_lengths[3];
		
	// send first to left then to right (up to 1 can be a real send)
	List* move_list;
	int neighbor_id;
	int i;
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { move_list = list_l; neighbor_id = id_l; }
		else if (i == 1) { move_list = list_r; neighbor_id = id_r; }
		
		if (neighbor_id == -1) continue;

		// send move_list, cells within will be converted into these buffers before send
		trees_parts_and_lengths = mpi_tree_send(move_list, neighbor_id, &buff_send_trees[i], &buff_send_parts[i], &buff_send_part_list_lengths[i], &lengths_of_buffs[i]);
		req_send_trees[i] = trees_parts_and_lengths[0];
		req_send_parts[i] = trees_parts_and_lengths[1];
		req_send_lengths[i] = trees_parts_and_lengths[2];
	}
		
	// buffer receive arrays
	simple_tree* buff_recv_trees[nNeighbors];
	particle* buff_recv_parts[nNeighbors];
	int* buff_recv_part_list_lengths[nNeighbors];
	
	MPI_Request req_recv_trees[nNeighbors];
	MPI_Request req_recv_parts[nNeighbors];
	MPI_Request req_recv_lengths[nNeighbors];
	
	// receive from neighbors: either buffers or nothing
	// receive first from left then from right (up to 2 can be real receives)
	for (i = 0; i < nNeighbors; ++i) { // nNeighbors = 2 = number of directions to send
		// assign left or right
		if (i == 0) { move_list = list_l; neighbor_id = id_l; }
		else if (i == 1) { move_list = list_r; neighbor_id = id_r; }
		
		if (neighbor_id == -1) continue;

		trees_parts_and_lengths = mpi_tree_recv(neighbor_id, &buff_recv_trees[i], &buff_recv_parts[i], &buff_recv_part_list_lengths[i], &lengths_of_buffs[i]);
		req_recv_trees[i] = trees_parts_and_lengths[0];
		req_recv_parts[i] = trees_parts_and_lengths[1];
		req_recv_lengths[i] = trees_parts_and_lengths[2];
	}
	
	// wait for receives
	MPI_Waitall(nNeighbors, req_recv_trees);
	MPI_Waitall(nNeighbors, req_recv_parts);
	MPI_Waitall(nNeighbors, req_recv_lengths);
	
	// once receives are done, can start unpacking buffers, putting new trees where they belong, while adjusting base_grid as necessary
	List* new_trees;
	tree* new_tree;
	double pxmax = pxmin+dx*(imax-imin); //pxmin is start of ghosts using global x position
	int j,k;
	for (neighbors)
		//new_trees = mpi_tree_unpack(&buff_recv_trees[i], &buff_recv_parts[i], &buff_recv_part_list_lengths[i], &lengths_of_buffs[i]);
		new_trees = mpi_tree_unpack(&buff_recv_trees[i], &buff_recv_parts[i], &buff_recv_part_list_lengths[i], &lengths_of_buffs[i]);
		list_reset_iter(new_trees);
		
		// first new_tree is used to check if base_grid is big enough
		if (list_has_next(*new_trees))
			new_tree = (tree*) list_get_next(new_trees);
		else
			next neighbor // no new_tree's to see here
		end
		
		// possibly change imin/imax plus possibly resize whole base_grid
		// new_tree only provides spatial to work with, so compare that to pxmin
		if ((new_tree->loc.x+dx/2) < pxmin)	//all new_tree's have same x. using dx/2 to prevent rounding issues
			if (imin == 0)
				resize_allocation(base_grid); //resets wi to 2*(imax-imin) and centers around (imin+imax)/2
			else if (imin < 0)
				printf("ERROR: unexpected imin < 0");
				MPI_Finalize();
				return -1;
			end
			--imin;
			pxmin -= dx;
		else if((new_tree->loc.x-dx/2) > pxmax)
			if (imax == wi-1)
				resize_allocation(base_grid);
			else if (imax > wi-1)
				printf("ERROR: unexpected imax > wi-1");
				MPI_Finalize();
				return -1;
			end
			++imax;
		end
		
		// we're good now: let's insert our new_tree's
		convert_ghost2real_and_reghost(base_grid, new_tree, dir6);
		while (list_has_next(*new_trees))
			new_tree = (tree*) list_get_next(new_trees);
			convert_ghost2real_and_reghost(base_grid, new_tree, dir6);
		end
		
	end
	
	// wait for sends
	MPI_Waitall(nNeighbors, req_send_trees);
	MPI_Waitall(nNeighbors, req_send_parts);
	MPI_Waitall(nNeighbors, req_send_lengths);
	
} //end give_and_take_cells


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
	base_grid = new_grid;
} //end resize_allocation


// TODO: write this
// dir6 is direction that tree was passed to get here, can be 'l','r','u','d','f','b'
// init ghosts are always present (for their points), so during convert we should treat them like NULLs and recreate them as necessary
//  - only when moving left, up, or forward are old init ghosts destroyed (not actually positive about this)
//  - when moving left, up, or forward, up to 23 new init ghosts are needed; otherwise up to 14 new init ghosts are needed
void convert_ghost2real_and_reghost(tree**** base_grid, tree* new_tree, char dir6) {
	int i = imin + 1 + (int) round((new_tree->loc.x - pxmin)/dx); //round to force correct int value
	int j = jmin + 1 + (int) round((new_tree->loc.y - pymin)/dy);
	int k = kmin + 1 + (int) round((new_tree->loc.z - pzmin)/dz);
	
	base_grid[i][j][k] = new_tree; //replaces the ghost that was there before
	new_tree->owner = pid;
	
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
	
	/*
	logic
	- points:
		- new_tree also needs to share points
	- 
		
	*/
	
	// create new ghosts
	// when dir6 is left, up, or forward, you're moving into init ghosts and so can use pre-existing ones in your layer. however you need to create a new layer of them.
	// when dir6 is right, down, or back, there are no init ghosts where you're going so you have to make ones in your layer. however there is no new layer to make.
	for (*d1 = -1; *d1 <= 1; ++*d1) {
		for (*d2 = -1; *d2 <= 1; ++*d2) {
			if (base_grid[i+di][j+dj][k+dk] == NULL) {
			
				// malloc
				base_grid[i+di][j+dj][k+dk] = (tree*) malloc(sizeof(tree));
				
				// tree_init, includes setting points[0] and owner
				*(base_grid[i+di][j+dj][k+dk]) = tree_init(get_loc(i+di,j+dj,k+dk), new_tree->neighbor_owners[1+*d1][1+*d2]);
				
			// if init ghost then convert into ghost instead of re-mallocing
			else if (base_grid[i+di][j+dj][k+dk]->owner == -2) {
				base_grid[i+di][j+dj][k+dk]->owner = new_tree->neighbor_owners[1+*d1][1+*d2];
				
			}
			// if neither NULL nor init ghost, this cell is already real or ghost and should be left alone
		}
	}
	// create new init ghosts. then make other 7 points point to neighbor points, but only for recently made ghosts
	int n;
	for (*d1 = -1; *d1 <= 1; ++*d1) {
		for (*d2 = -1; *d2 <= 1; ++*d2) {
			// recently made ghosts can be identified by points[1] == NULL (which is specifically set in tree_init)
			// second condition is to avoid choosing init ghosts as well
			if (base_grid[i+di][j+dj][k+dk]->root->points[1] == NULL && base_grid[i+di][j+dj][k+dk]->owner != -2) {
			
				// create new init ghosts as necessary. n should be thought of as binary for getX etc. (n for "neighbors")
				for (n = 1; n < 8; ++n) {
					if (base_grid[i+di+getX(n)][j+dj+getY(n)][k+dk+getZ(n)] == NULL) {
					
						// malloc
						base_grid[i+di+getX(n)][j+dj+getY(n)][k+dk+getZ(n)] = (tree*) malloc(sizeof(tree));
				
						// tree_init, includes setting points[0] and owner = -2
						*(base_grid[i+di+getX(n)][j+dj+getY(n)][k+dk+getZ(n)]) = tree_init(get_loc(i+di+getX(n),j+dj+getY(n),k+dk+getZ(n)), -2);
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
// All Pseudocode in here

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