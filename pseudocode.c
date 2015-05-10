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

// TODO: needs some work on properly separating the tasks onto their appropriate procs, i.e. giving vs taking

// this one may need more?
trade_cells(base_grid)
	for each (direction) // direction = x, y, z
		give_and_take_cells(direction, base_grid)
	end
end trade_cells

// TODO: written only for x direction, need to generalize to y and z
// determine which trees to send where, and send them. then take trees given to you, adjust your base_grid if necessary, and insert the trees
give_and_take_cells(direction, tree**** base_grid)
	int neighbors[numNeighbors] = proc neighbors
	
	// create moving_trees: an array of List*s pointing to tree*s. The use of tree*s is not made explicit here, may need to cast later
	// always List*, never just List, for C reasons
	List* moving_trees[sizeof(neighbors)]
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
	
	// need a neighbors-sized buffer array filled with simple_tree/particle arrays. both get malloced in mpi_tree_send
	simple_tree* buff_send_trees[sizeof(neighbors)]
	particle* buff_send_parts[sizeof(neighbors)]
	int* buff_send_part_list_lengths[sizeof(neighbors)]
	
	// need to track MPI_request objects for a waitall later, in order to be non-blocking
	MPI_Request req_send_trees[sizeof(neighbors)]
	MPI_Request req_send_parts[sizeof(neighbors)]
	MPI_Request req_send_lengths[sizeof(neighbors)]
	MPI_Request trees_parts_and_lengths[3]

	// send moving_trees, they will be converted into buffers before send
	for (neighbors)
		trees_parts_and_lengths = mpi_tree_send(moving_trees[neighbor_id], neighbor_id, &buff_send_trees[neighbor_id], &buff_send_parts[neighbor_id], &buff_send_part_list_lengths[neighbor_id])
		req_send_trees[neighbor_id] = trees_parts_and_lengths[0]
		req_send_parts[neighbor_id] = trees_parts_and_lengths[1]
		req_send_lengths[neighbor_id] = trees_parts_and_lengths[2]
	end
	
	simple_tree* buff_recv_trees[sizeof(neighbors)]
	particle* buff_recv_parts[sizeof(neighbors)]
	int* buff_recv_part_list_lengths[sizeof(neighbors)]
	int lengths_of_buffs[sizeof(neighbors)][3] //one length for each of the above 3 buffers
	
	MPI_Request req_recv_trees[sizeof(neighbors)]
	MPI_Request req_recv_parts[sizeof(neighbors)]
	MPI_Request req_recv_lengths[sizeof(neighbors)]
	
	// receive from neighbors: either buffers or nothing
	for (neighbors)
		trees_parts_and_lengths = mpi_tree_recv(neighbor_id, &buff_recv_trees[neighbor_id], &buff_recv_parts[neighbor_id], &buff_recv_part_list_lengths[neighbor_id], &lengths_of_buffs[neighbor_id])
		req_recv_trees[neighbor_id] = trees_parts_and_lengths[0]
		req_recv_parts[neighbor_id] = trees_parts_and_lengths[1]
		req_recv_lengths[neighbor_id] = trees_parts_and_lengths[2]
	end
	
	// wait for receives
	MPI_Waitall(sizeof(neighbors), req_recv_trees);
	MPI_Waitall(sizeof(neighbors), req_recv_parts);
	MPI_Waitall(sizeof(neighbors), req_recv_lengths);
	
	// once receives are done, can start unpacking buffers, putting new trees where they belong, while adjusting base_grid as necessary
	List* new_trees;
	tree* new_tree;
	double pxmax = pxmin+dx*(imax-imin);
	int j,k;
	for (neighbors)
		//new_trees = mpi_tree_unpack(&buff_recv_trees[neighbor_id], &buff_recv_parts[neighbor_id], &buff_recv_part_list_lengths[neighbor_id], &lengths_of_buffs[neighbor_id]);
		new_trees = mpi_tree_unpack(&buff_recv_trees[neighbor_id], &buff_recv_parts[neighbor_id], &buff_recv_part_list_lengths[neighbor_id], &lengths_of_buffs[neighbor_id]);
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
		convert_ghost2real_and_reghost(base_grid, new_tree);
		while (list_has_next(*new_trees))
			new_tree = (tree*) list_get_next(new_trees);
			convert_ghost2real_and_reghost(base_grid, new_tree);
		end
		
	end
	
	// wait for sends
	MPI_Waitall(sizeof(neighbors), req_send_trees);
	MPI_Waitall(sizeof(neighbors), req_send_parts);
	MPI_Waitall(sizeof(neighbors), req_send_lengths);
	
end give_and_take_cells


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
				if (i >= imin && i < imax &&
					j >= jmin && j < jmax &&
					k >= kmin && k < kmax) {
					
					new_grid[i][j][k] = base_grid[i-di][j-dj][k-dk]; //translating old tree*s onto new grid, includes some NULLs
					
				} else {
					new_grid[i][j][k] = NULL;
				}
			}
		}
	}
	base_grid = new_grid;
} //end resize_allocation


// TODO: write this. should it even be a separate function or not?
void convert_ghost2real_and_reghost(tree**** base_grid, tree* new_tree) {
	int i = imin + 1 + (int) round((new_tree->loc.x - pxmin)/dx); //round to force correct int value
	int j = jmin + 1 + (int) round((new_tree->loc.y - pymin)/dy);
	int k = kmin + 1 + (int) round((new_tree->loc.z - pzmin)/dz);
	base_grid[i][j][k] = new_tree; //replaces the ghost that was there before
	new_tree->owner = pid;
	
	// create new ghosts around new_tree as needed
	for (neighbor trees of new_tree)
		if (neighbor == NULL)
			init cell	// find owner somehow
			laser(cell)
		end
	end
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