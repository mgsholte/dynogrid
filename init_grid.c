grid_point** init_grid(int nx, int ny) {
	grid_point** grid;
	for (int i = 0; i < ny; i++) {
		grid[i] = (grid_point*) malloc (nx*sizeof(grid_point));
		for (int j = 0; j < nx; j++) {
			grid[i][j].E = {0,0,0};
			grid[i][j].B = {0,0,0};
			// Assumes no more info at each grid point
		}
	}
	return grid;
}