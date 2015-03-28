#include "grid.h"
#include <stdlib.h>

// return a random double value in the range [low,high]
static double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

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

// inits rectangular region of particles, approx evenly distributed
// upper left corner ul
// lower right corner lr
List init_particles(vec2 ul, vec2 lr, int part_per_cell) {
	List particles = list_init();
	for (int row = lr.y; row < ul.y; row++) {
		for (int col = ul.x; col < lr.x; col++) {
			// add protons
			for (int k = 0; k < part_per_cell/2; k++) {
				particle* p = (particle *)malloc(sizeof(particle));
				double x = rand_float(0, dx) + j;
				double y = rand_float(0, dy) + i;
				p = { {x, y},		//position
					  {0, 0, 0},	//momentum
					  BASE_PROTON_MASS * PROTON_WEIGHT		//mass
					  BASE_PROTON_CHARGE * PROTON_WEIGHT	//charge
					  PROTON_WEIGHT		//weight
					};
				list_add(particles, p);
			}
			// add electrons
			for (int k = 0; k < (part_per_cell - part_per_cell/2); k++) {
				particle *p = (particle *)malloc(sizeof(particle));
				double x = rand_float(0, dx) + j;
				double y = rand_float(0, dy) + i;
				p = { {x, y},		//position
					  {0, 0, 0},	//momentum
					  BASE_ELECTRON_MASS * ELECTRON_WEIGHT		//mass
					  BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT	//charge
					  ELECTRON_WEIGHT		//weight
					};
				list_add(particles, p);
			}
		}
	}
	return particles;
}
