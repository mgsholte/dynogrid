void Balance(tree ****grid){
	tree *curCell = NULL;
	int partwork = 0;
	int cellwork = 0;
	double work, avgwork, mostwork;
	// Calculate amount of work on this processor
	for (i = imin; i < imax; ++i) {
		for (j = jmin; j < jmax; ++j) {
			for (k = kmin; k < kmax; ++k) {
				curCell = grid[i][j][k];
				if (curCell != NULL) {
					if (curCell->owner == pid) {
						partwork += list_length(curCell->part_list);
					}
					cellwork += curCell->numpoints;
				}
			}
		}
	}
	work = .8*partwork + .2*cellwork;

	// Find out how much work there is in total
	MPI_Allredice(&work, &avgwork, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	avgwork = avgwork / numprocs;

	MPI_Allreduce(&work, &mostwork, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	while (mostwork > 1.1*avgwork){
		propensity = work - 1.05*avgwork;
		//Tell neighbors propensity
		//get propensity from neighbors
		if (propensity > 0){
			//calculate propensity left and right
			//Give left or right and give null other way
			//update work based on what was given
		}
		else{
			//Give null left and right
		}
		//Recieve from left and right and update work
		MPI_Allreduce(&work, &mostwork, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
}
