#include <stdlib.h>

#include "grid.h"
#include "output.h"

// return a random double value in the range [low,high]
static double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

// inits all grid points to 0 in E and B
grid_point** init_grid(int nx, int ny) {
	grid_point **grid_points = (grid_point**) malloc( ny * sizeof(size_t) ); // allocate an array of pointers to rows
	int i; int j;
	for (i = 0; i < ny; ++i) {
		grid_points[i] = (grid_point*) malloc( nx * sizeof(grid_point) );  // allocate the row
		for (j = 0; j < nx; j++) {
			// all grid points are initialized with no field
			(grid_points[i][j]).E = (vec3) {0,0,0};
			(grid_points[i][j]).B = (vec3) {0,0,0};
		}
	}
	return grid_points;
}

// inits rectangular region of particles, approx evenly distributed
// upper left corner ul
// lower right corner lr
List init_particles(vec2 ul, vec2 lr, int part_per_cell) {
	List particles = list_init();
	int row; int col; int k;
	for (row = ul.y; row < lr.y; row++) {
		for (col = ul.x; col < lr.x; col++) {
			// add protons
			for (k = 0; k < part_per_cell/2; k++) {
				particle *p = (particle*) malloc(sizeof(particle));
				double x = rand_float(0, dx) + col;
				double y = rand_float(0, dy) + row;
				*p = (particle) { {x, y},		//position
					  {0, 0, 0},	//momentum
					  BASE_PROTON_MASS * PROTON_WEIGHT,		//mass
					  BASE_PROTON_CHARGE * PROTON_WEIGHT,	//charge
					  PROTON_WEIGHT		//weight
					};
				list_add(&particles, p);
			}
			// add electrons
			for (k = 0; k < (part_per_cell - part_per_cell/2); k++) {
				particle *p = (particle*) malloc(sizeof(particle));
				double x = rand_float(0, dx) + col;
				double y = rand_float(0, dy) + row;
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
	list_reset_iter(&particles);
	return particles;
}

void output_grid(int itNum, grid_point **grid_points, int nx, int ny, List particles) {
	output_data2D(itNum, grid_points, nx, ny, particles);
	//TODO: it should be true that itNum == time/dt. maybe we don't need to pass the itNum variable as an argument
	//#define round(x) (int) (x+0.5)
	// int itNum = round(time/dt);
}
          
