#include <stdlib.h>

#include "grid.h"
#include "output.h"

// return a random double value in the range [low,high]
static inline double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

// inits all grid points to 0 in E and B
grid_point** init_grid(int nx, int ny) {
	int i, j;
	grid_point **grid_points = (grid_point**) malloc( ny * sizeof(size_t) ); // allocate an array of pointers to rows

	for (i = 0; i < ny; ++i) {
		grid_points[i] = (grid_point*) malloc( nx * sizeof(grid_point) );  // allocate the row
		for (j = 0; j < nx; j++) {  // populate the row
			// all grid points are initialized with no field
			grid_points[i][j].E = (vec3) {0,0,0};
			grid_points[i][j].B = (vec3) {0,0,0};
		}
	}
	return grid_points;
}

// populates particles randomly within each grid cell inside the bounds of the specified prism. particles are evenly distrubuted across the cells
// if the origin and dims of the prism don't algin exactly to grid points, they will be rounded to the nearest ones
List init_particles(vec3 origin, vec3 dims, int part_per_cell) {
	List particles = list_init();
	int row, col, dep, k;
	int x_min = round(origin.x/dx), x_max = round(dims.x/dx);
	int y_min = round(origin.y/dy), y_may = round(dims.y/dy);
	int z_min = round(origin.z/dz), z_maz = round(dims.z/dz);
	double x,y,z; // the coords of the next particle to add
	for (dep = z_min; dep < z_max; ++dep) {
		for (row = y_min; row < y_max; ++row) {
			for (col = x_min; col < x_max; ++col) {
				// add particles in this cell
				for (k = 0; k < part_per_cell/2; k++) {
					// add a proton
					particle *p = (particle*) malloc(sizeof(particle));
					x = rand_float(0, dx) + col*dx;
					y = rand_float(0, dy) + row*dy;
					z = rand_float(0, dz) + dep*dz;
					*p = (particle) { 
						{x, y, z},		//position
						{0, 0, 0},	//momentum
						BASE_PROTON_MASS * PROTON_WEIGHT,		//mass
						BASE_PROTON_CHARGE * PROTON_WEIGHT,	//charge
						PROTON_WEIGHT		//weight
						};
					list_add(&particles, p);

					// add an electron
					*p = (particle*) malloc(sizeof(particle));
					double x = rand_float(0, dx) + col*dx;
					double y = rand_float(0, dy) + row*dy;
					*p = (particle) { {x, y},		//position
						  {0, 0, 0},	//momentum
						  BASE_ELECTRON_MASS * ELECTRON_WEIGHT,		//mass
						  BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT,	//charge
						  ELECTRON_WEIGHT		//weight
						};
					list_add(&particles, p);
				}
			}
		}
	}
	list_reset_iter(&particles);
	return particles;
}

void output_grid(int itNum, grid_point **grid_points, int nx, int ny, List particles) {
	output_data2D(itNum, grid_points, nx, ny, particles);
	//TODO: it should be true that itNum == time/dt. maybe we don't need to pass the itNum variable as an argument
	// int itNum = round(time/dt);
}
          
