#include "decs.h"
#include "grid.h"

int main(int argc, char *argv[]) {
	// read as inputs in the future
	int nSteps = 1000;
	int nx = 100, ny = 100;
	int part_per_cell = 10;
	int output_freq = nSteps/10;

	double **grid_points = init_grid(nx, ny);
	particle *particles = init_particles(part_per_cell);

	int i;
	for(i = 0; i < nSteps; ++i) {
		push_particles(grid_points, particles);
		update_grid(i, grid_points);  // add the laser, ...
		if (i % output_freq == 0) {
			output_grid(grid_points, particles);
		}
	}
	output_grid(grid_points, particles);
}
