#include <stdlib.h>

#include "decs.h"
#include "grid.h"
#include "output.h"

// return a random double value in the range [low,high]
static inline double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

// inits all grid points to 0 in E and B
grid_cell*** init_grid() {
	grid_cell ***grid_cells = (grid_cell***) malloc( (nx+1) * sizeof(grid_cell**) ); // allocate an array of pointers to rows-depthwise
	int i, j, k, n;
	for (i = 0; i <= nx; ++i) {
		grid_cells[i] = (grid_cell**) malloc( (ny+1) * sizeof(grid_cell*) );  // allocate the row
		for (j = 0; j <= ny; ++j) {
			grid_cells[i][j] = (grid_cell*) malloc( (nz+1) * sizeof(grid_cell) );  // allocate the row
			// initialize only the upper-left-forward grid_point (points[0] represents (x,y,z)==000)
			for (k = 0; k <= nz; ++k) {
				grid_cells[i][j][k].points[0] = malloc( sizeof(grid_point) );
				(grid_cells[i][j][k].points[0])->E = (vec3) {0,0,0};
				(grid_cells[i][j][k].points[0])->B = (vec3) {0,0,0};
				grid_cells[i][j][k].children = NULL;
			}
			// make other 7 points point to neighbors
			for (k = 0; k < nz; ++k) {
				if (i<nx && j<ny && k<nz) {
					// n should be thought of as binary (n for "neighbors")
					for (n = 1; n < 8; ++n) {
						grid_cells[i][j][k].points[n] = grid_cells[i+(n&1)][j+(n&2)/2][k+(n&4)/4].points[0];
					}
				}
			}
		}
	}
	return grid_cells;
}

// populates particles randomly within each grid cell inside the bounds of the specified prism. particles are evenly distrubuted across the cells
// if the origin and dims of the prism don't algin exactly to grid points, they will be rounded to the nearest ones
List init_particles(vec3 origin, vec3 dims, int part_per_cell) {
	List particles = list_init();
	int row, col, dep, k;
	int x_i_min = round_i(origin.x/dx), x_i_max = round_i(dims.x/dx) + x_i_min;
	int y_i_min = round_i(origin.y/dy), y_i_max = round_i(dims.y/dy) + y_i_min;
	int z_i_min = round_i(origin.z/dz), z_i_max = round_i(dims.z/dz) + z_i_min;
	double x,y,z; // the coords of the next particle to add
	particle *p; // pointer to the next particle to add
	for (col = x_i_min; col < x_i_max; ++col) {
		for (row = y_i_min; row < y_i_max; ++row) {
			for (dep = z_i_min; dep < z_i_max; ++dep) {
				// add particles in this cell
				for (k = 0; k < part_per_cell/2; ++k) {
					// add a proton
					p = (particle*) malloc(sizeof(particle));

					x = rand_float(0, dx) + col*dx;
					y = rand_float(0, dy) + row*dy;
					z = rand_float(0, dz) + dep*dz;
					*p = (particle) { 
						{x, y, z},  //position
						{0, 0, 0},  //momentum
						BASE_PROTON_MASS * PROTON_WEIGHT,   //mass
						BASE_PROTON_CHARGE * PROTON_WEIGHT, //charge
						PROTON_WEIGHT  //weight
						};

					list_add(&particles, p);

					// add an electron
					p = (particle*) malloc(sizeof(particle));

					x = rand_float(0, dx) + col*dx;
					y = rand_float(0, dy) + row*dy;
					z = rand_float(0, dz) + dep*dz;
					*p = (particle) { 
						{x, y, z},  //position
						{0, 0, 0},  //momentum
						BASE_ELECTRON_MASS * ELECTRON_WEIGHT,   //mass
						BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT, //charge
						ELECTRON_WEIGHT  //weight
						};

					list_add(&particles, p);
				}
			}
		}
	}
	list_reset_iter(&particles);
	return particles;
}

void output_grid(int itNum, int numFiles, grid_point ***grid_points, List particles) {
	output_data3D(itNum, numFiles, grid_points, particles);
	//TODO: it should be true that itNum == time/dt. maybe we don't need to pass the itNum variable as an argument
	// int itNum = round_i(time/dt);
}
          
