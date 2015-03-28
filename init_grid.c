#include "grid.h"

grid_point** init_grid(int nx, int ny) {
	grid_point** grid = (grid_point **)malloc( ny * sizeof(size_t) ); // allocate an array of pointers to rows
	for (int i = 0; i < ny; ++i) {
		grid[i] = (grid_point *)malloc( nx * sizeof(grid_point) );  // allocate the row
		for (int j = 0; j < nx; j++) {
			// all grid points are initialized with no field
			grid[i][j].E = {0,0,0};
			grid[i][j].B = {0,0,0};
		}
	}
	return grid;
}
