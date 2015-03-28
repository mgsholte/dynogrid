#include <stdio.h>

#include "decs.h"
#include "grid.h"

int main(int argc, char *argv[]) {
	int i;  // loop index variable

	// read as inputs in the future
	int nSteps = 1000;
	int nx = 100, ny = 100;
	int part_per_cell = 10;
	int output_freq = nSteps/10;

	printf("initializing grid and particles");

	double **grid_points = init_grid(nx, ny);
	particle *particles = init_particles(part_per_cell);

	printf("beginning simulation");

	for(i = 0; i < nSteps; ++i) {
		time = i*dt;
		push_particles(grid_points, particles);
		update_grid(i, grid_points);  // add the laser, ...
		if (i % output_freq == 0) {
			output_grid(i, grid_points, particles);
		}
	}
	output_grid(i, grid_points, particles);

	printf("simulation finished");

	return 0;
}