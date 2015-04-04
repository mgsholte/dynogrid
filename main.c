#include <stdio.h>

#include "decs.h"
#include "grid.h"
#include "dynamics.h"
#include "list.h"

int main(int argc, char *argv[]) {
	int i;  // loop index variable

	// read as inputs in the future
	int nSteps = 1000;
	//int nx = 100, ny = 100, nz = 100;
	int part_per_cell = 10;
	int output_freq = nSteps/10;
	//vec2 ul = {45, 45}, lr = {55, 55};  // upper left, lower right coordinates defining the rectangle where particles begin in the simulation
	vec3 ulf = {.4*x_max, .4*y_max, .4*z_max};
	vec3 dims = {.2*x_max, .2*y_max, .2*z_max};

	printf("initializing grid and particles\n");

	grid_point ***grid_points = init_grid();
	List particles = init_particles(ulf, dims, part_per_cell);

	printf("finished initializing. beginning simulation\n");

	for(i = 0; i < nSteps; ++i) {
		time = i*dt;
		push_particles(grid_points, particles);
		update_grid(grid_points);  // add the laser, ...
		if (i % output_freq == 0) {
			output_grid((i/output_freq), (nSteps/output_freq), grid_points, particles);
		}
	}
	output_grid((i/output_freq), (nSteps/output_freq), grid_points, particles);

	printf("simulation finished\n");

	return 0;
}
