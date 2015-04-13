#include <stdio.h>

#include "decs.h"
#include "grid.h"
#include "dynamics.h"
#include "list.h"

int main(int argc, char *argv[]) {
	int i;  // loop index variable

	// read as inputs in the future
	int nSteps = ceil(dt/t_end);
	//int nx = 100, ny = 100, nz = 100;
	int part_per_cell = 5;
	// print output every output_freq iterations of the main loop
	int output_freq = nSteps/10;
	// upper left coordinate and dimensions defining the rectangle where particles begin in the simulation
	vec3 ulf = {.46*x_max, .46*y_max, .46*z_max};
	vec3 dims = {.08*x_max, .08*y_max, .08*z_max};

	printf("initializing grid and particles\n");

	grid_point ***grid_points = init_grid();
	List particles = init_particles(ulf, dims, part_per_cell);

	printf("finished initializing. beginning simulation\n");

	for(i = 0; i < nSteps; ++i) {
		time = i*dt;
		push_particles(grid_points, particles);
		update_grid(grid_points);  // add the laser, etc.
		if (i % output_freq == 0) {
			output_grid((i/output_freq), (nSteps/output_freq), grid_points, particles);
		}
	}

	// print final state unless it was already output on the last iteration of the above loop
	if (i % output_freq != 0)
		output_grid((i/output_freq), (nSteps/output_freq), grid_points, particles);
	 
	printf("simulation finished\n");

	return 0;
}
