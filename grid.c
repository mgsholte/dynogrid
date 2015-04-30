#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "decs.h"
#include "grid.h"
#include "output.h"
#include "dynamics.h"

// return a random double value in the range [low,high]
static inline double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

static inline vec3 get_loc(int ix, int iy, int iz) {
	return (vec3) { ix*dx, iy*dy, iz*dz };
}

// inits all grid points to 0 in E and B
tree*** grid_init() {
	int i, j, k, n;  // loop index vars
	tree ***base_grid = (tree***) malloc( (nx+1) * sizeof(tree**) );  // a 3-dim array holding the roots of the base grid
	// initialize all the trees. note that 1 row in each dimension is a row of 'ghost' trees. these are responsible for allocating and freeing grid_point pointers and nothing else
	for (i = 0; i <= nx; ++i) {
		base_grid[i] = (tree**) malloc( (ny+1) * sizeof(tree*) );
		for (j = 0; j <= ny; ++j) {
			base_grid[i][j] = (tree*) malloc( (nz+1) * sizeof(tree) );
			for (k = 0; k <= nz; ++k) {
				base_grid[i][j][k] = tree_init(get_loc(i,j,k));
			}
			// each tree allocated only 1 grid_point. now set all trees (except ghosts) to point to the points allocated by their neighbors
			for (k = 0; k < nz; ++k) {
				for (n = 1; n < 8; ++n) {
					base_grid[i][j][k].root->points[n] = base_grid[i+(n&1)][j+(n&2)/2][k+(n&4)/4].root->points[0];
				}
			}
		}
	}
	return base_grid;
}

// populates particles randomly within each grid cell inside the bounds of the specified prism. particles are evenly distrubuted across the cells
// if the origin and dims of the prism don't algin exactly to grid points, they will be rounded to the nearest ones
void init_particles(tre base_grid, vec3 origin, vec3 dims, int elec_per_cell) {
	List particles = list_init();
	int iy, ix, iz, k;
	int ix_min = round_i(origin.x/dx), ix_max = round_i(dims.x/dx) + ix_min;
	int iy_min = round_i(origin.y/dy), iy_max = round_i(dims.y/dy) + iy_min;
	int iz_min = round_i(origin.z/dz), iz_max = round_i(dims.z/dz) + iz_min;
	double x,y,z; // the coords of the next particle to add
	particle *p; // pointer to the next particle to add
	for (ix = ix_min; ix < ix_max; ++ix) {
		for (iy = iy_min; iy < iy_max; ++iy) {
			for (iz = iz_min; iz < iz_max; ++iz) {
				// add particles in this cell
				for (k = 0; k < elec_per_cell; ++k) {
					tree cell = base_grid[ix][iy][iz];
					// add a proton
					p = (particle*) malloc(sizeof(particle));

					x = rand_float(0, dx) + cell.loc.x;
					y = rand_float(0, dy) + cell.loc.y;
					z = rand_float(0, dz) + cell.loc.z;
					*p = (particle) { 
						{x, y, z},  //position
						{0, 0, 0},  //momentum
						BASE_PROTON_MASS * PROTON_WEIGHT,   //mass
						BASE_PROTON_CHARGE * PROTON_WEIGHT, //charge
						PROTON_WEIGHT  //weight
						};

					list_add(cell.particles, p);

					// add an electron
					p = (particle*) malloc(sizeof(particle));

					x = rand_float(0, dx) + cell.loc.x;
					y = rand_float(0, dy) + cell.loc.y;
					z = rand_float(0, dz) + cell.loc.z;
					*p = (particle) { 
						{x, y, z},  //position
						{0, 0, 0},  //momentum
						BASE_ELECTRON_MASS * ELECTRON_WEIGHT,   //mass
						BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT, //charge
						ELECTRON_WEIGHT  //weight
						};

					list_add(cell.particles, p);
					list_reset_iter(&cell.particles);
				}
			}
		}
	}
}

void output_grid(int itNum, int numFiles, tree ***base_grid) {
	output_grid_impl(itNum, numFiles, base_grid, "data");
	//TODO: it should be true that itNum == time/dt. maybe we don't need to pass the itNum variable as an argument
	// int itNum = round_i(time/dt);
}
void recursive_execute_coarsen(grid_cell* cell);
void cleanup(grid_cell ***grid_cells) {
	int x,y,z,i;
	grid_cell* cell;
	for (x = 0; x <= nx; x++) {
		for (y = 0; y <= ny; y++) {
			for (z = 0; z <= nz; z++) {
				// coarsen completely (coarsest should have no children after)
				cell = &grid_cells[x][y][z];
				recursive_execute_coarsen(cell);
				
				// only free what was malloc'd per cell in init_grid
				free(cell->points[0]);
			}
			// frees cells and cell.children address
			free(grid_cells[x][y]);
		}
		free(grid_cells[x]);
	}
	free(grid_cells);
}
void recursive_execute_coarsen(grid_cell* cell) {
	// BASE CASE or if at coarsest cell size from the start
	if (cell->children == NULL) {
		return;
	}
	else {
		int child_num;
		for(child_num = 0; child_num < 8; child_num++) {
			recursive_execute_coarsen(cell->children[child_num]);
		}
		execute_coarsen(cell);
	}
}
