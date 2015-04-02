#include <stdio.h>

#include "decs.h"
#include "grid.h"
#include "dynamics.h"
#include "list.h"

int main(int argc, char *argv[]) {
	int i;  // loop index variable

	// read as inputs in the future
	int nSteps = 1000;
	int nx = 100, ny = 100;
	int part_per_cell = 10;
	int output_freq = nSteps/10;
	vec2 ul = {45, 45}, lr = {55, 55};  // upper left, lower right coordinates defining the rectangle where particles begin in the simulation

	printf("initializing grid and particles\n");

	grid_point **grid_points = init_grid(nx, ny);
	List particles = init_particles(ul, lr, part_per_cell);

	output_data2D(1, grid_points, nx, ny, particles);

	printf("beginning simulation\n");

	for(i = 0; i < nSteps; ++i) {
		time = i*dt;
		push_particles(grid_points, particles);
		update_grid(grid_points);  // add the laser, ...
		if (i % output_freq == 0) {
			output_data2D(i, grid_points, nx, ny, particles);
		}
	}
	output_grid(i, grid_points, nx, ny, particles);

	printf("simulation finished\n");

	return 0;
}
