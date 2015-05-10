#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "tree.h"
#include "decs.h"
#include "grid.h"
#include "output.h"
#include "dynamics.h"

// return a random double value in the range [low,high]
static inline double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

static inline vec3 get_loc(int ix, int iy, int iz) {
	return (vec3) { (ix-imin)*dx + pxmin, (iy-jmin)*dy + pymin, (iz-kmin)*dz + pzmin };
}

static inline int max(const int x, const int y) {
	return x > y ? x : y;
}

static inline int min(const int x, const int y) {
	return x < y ? x : y;
}

// inits all grid points with E=B=0 for all of them
// isize, etc. are lengths of real cell dimensions. ghost and NULL cells are added in this function
// imin|imax, etc. are the start|end indexes for the real+ghost cell dims
// padding is 50% of isize (etc.) on each side
// ASSUMING GLOBAL numProcs, x_divs, y_divs, z_divs (the user-supplied specs)
tree**** grid_init(int isize, int jsize, int ksize, int x_divs, int y_divs, int z_divs) {
	// global vars (local indices)
	imin = isize/2 - 1; // -1 for ghost on left
	jmin = jsize/2 - 1;
	kmin = ksize/2 - 1;
	imax = imin + isize + 2; // +2 for ghosts on left and right
	jmax = jmin + jsize + 2;
	kmax = kmin + ksize + 2;
	// global vars (global spatial positions)
	pxmin = ((double)(pid % x_divs))/x_divs * x_max - dx;
	pymin = ((double)((pid - pid % x_divs)/x_divs % y_divs))/y_divs * y_max - dy;
	pzmin = ((double)((pid - pid % x_divs - pid % (x_divs * y_divs))/(x_divs * y_divs)))/z_divs * z_max - dz;
	
	// local vars
	int wi = 2*isize; // padding is 50% of isize (etc.) on each side
	int wj = 2*jsize;
	int wk = 2*ksize;
	
	tree ****base_grid = (tree****) malloc( (wi+1) * sizeof(tree***) ); // allocate an array of pointers to rows-depthwise
	int i, j, k, n;
	for (i = 0; i < wi; ++i) {
		base_grid[i] = (tree***) malloc( (wj+1) * sizeof(tree**) );  // allocate the row
		for (j = 0; j < wj; ++j) {
			base_grid[i][j] = (tree**) malloc( (wk+1) * sizeof(tree*) );  // allocate the row
			// each index goes from (ghost cell on left) to (initialization ghost cell on right) which is one past (ghost cell on right) for the
			// purpose of initializing all points. each cell makes just one point, so an extra layer is needed at the end
			for (k = 0; k < wk; ++k) {
				// malloc for real, ghost, and 'init ghost' cells
				if (i >= imin && i <= imax &&
					j >= jmin && j <= jmax &&
					k >= kmin && k <= kmax) {
					
					base_grid[i][j][k] = (tree*) malloc(sizeof(tree));
					
					*(base_grid[i][j][k]) = tree_init(get_loc(i,j,k), pid);

					// 26 possible neighbors could own each ghost cell, but their pids can be constructed from true/false statements. using true->1 and false->0
					int di, dj, dk, owner_id;
					di = (i == imax-1) - (i == imin); // +1 or -1
					dj = (j == jmax-1) - (j == jmin);
					dk = (k == kmax-1) - (k == kmin);
					owner_id = pid;
					owner_id += di;
					owner_id += dj * x_divs;
					owner_id += dk * x_divs * y_divs;
					
					// find bad cases where there is no proper owner_id, i.e. the above algorithm got a bad answer b/c ghost is out of simulation bounds
					if ((di == 1  &&  pxmin + dx * (isize+1.5) >= x_max) || 	// 1.5 is to prevent rounding issues, even though 1 should be enough
						(dj == 1  &&  pymin + dy * (isize+1.5) >= y_max) ||
						(dk == 1  &&  pzmin + dz * (isize+1.5) >= z_max) ||
						(di == -1  &&  pxmin + dx * 0.5 <= 0) ||
						(dj == -1  &&  pymin + dy * 0.5 <= 0) ||
						(dk == -1  &&  pzmin + dz * 0.5 <= 0)) {

						owner_id = -1;
					}
					
				} else {
					base_grid[i][j][k] = NULL;
				}
			}
		}
	}
	// make other 7 grid_point pointers of each grid cell point to grid points allocated by neighbors. except for fake 'init ghost' cells on right/bottom/back boundaries
	for (i = imin; i < imax; ++i) {
		for (j = jmin; j < jmax; ++j) {
			for (k = kmin; k < kmax; ++k) {
				// n should be thought of as binary (n for "neighbors")
				for (n = 1; n < 8; ++n) {
					base_grid[i][j][k]->root->points[n] = base_grid[i+getX(n)][j+getY(n)][k+getZ(n)]->root->points[0];
				}
			}
		}
	}
	return base_grid;
}

// populates particles randomly within each grid cell inside the bounds of the specified prism. particles are evenly distrubuted across the cells
// if the origin and dims of the prism don't algin exactly to grid points, they will be rounded to the nearest ones
void init_particles(tree ****base_grid, vec3 origin, vec3 dims, int elec_per_cell) {
	//TODO: this assumes that ix=0 gives the cell with origin at x=0. this is not necessarily the case looking at grid_init(). each proc needs to know an offset where its grid begins. i.e. ix=0 on my local grid corresponds to ix=offset on the global grid
	int iy, ix, iz, n;
	int ix_min = round_i((origin.x-pxmin)/dx)+imin, ix_max = round_i(dims.x/dx) + ix_min;
	int iy_min = round_i((origin.y-pymin)/dy)+jmin, iy_max = round_i(dims.y/dy) + iy_min;
	int iz_min = round_i((origin.z-pzmin)/dz)+kmin, iz_max = round_i(dims.z/dz) + iz_min;
	// each processor starts no lower than the cells it is responsible for, also no lower than the origin of the prism containing the particles
	// also, exclude ghost cells
	ix_min = max(ix_min, imin+1);
	iy_min = max(iy_min, jmin+1);
	iz_min = max(iz_min, kmin+1);
	// each processor goes no further than the cells it is responsible for, also no furthen than the end of the prism containing the particles
	// also, exclude ghost cells
	ix_max = min(ix_max, imax-1);
	iy_max = min(iy_max, jmax-1);
	iz_max = min(iz_max, kmax-1);
	//NB: it is possible that a processor might not be responsible for adding any particles, in which case the body of the loops below will not execute
	for (ix = ix_min; ix < ix_max; ++ix) {
		for (iy = iy_min; iy < iy_max; ++iy) {
			for (iz = iz_min; iz < iz_max; ++iz) {
				tree cell = *base_grid[ix][iy][iz];
				// add particles in this cell
				for (n = 0; n < elec_per_cell; ++n) {
					double x,y,z; // the coords of the next particle to add
					particle *p; // pointer to the next particle to add

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

					list_add(&cell.particles, p);

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

					list_add(&cell.particles, p);
				}
				// reset the iterator after adding the particles. probably redundant since this should be called anytime before starting an iteration
				list_reset_iter(&cell.particles);
			}
		}
	}
}

void output_grid(int itNum, int numFiles, tree ****base_grid) {
	output_grid_impl(itNum, numFiles, base_grid, "data");
	//TODO: it should be true that itNum == time/dt. maybe we don't need to pass the itNum variable as an argument
	// int itNum = round_i(time/dt);
}

/* haven't yet integrated this function with the tree system
void recursive_execute_coarsen(grid_cell* cell);
//TODO: change these loop bounds to whatever suits the changing grid chunk size. not just imin to imax, right? there's uninitialized cells.
void cleanup(grid_cell ****grid_cells) {
	int x,y,z,i;
	grid_cell* cell;
	for (x = 0; x <= nx; x++) {
		for (y = 0; y <= ny; y++) {
			for (z = 0; z <= nz; z++) {
				// coarsen completely (coarsest should have no children after)
				cell = grid_cells[x][y][z];
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
*/
