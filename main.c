#include <stdio.h>
#include <math.h>

#include "decs.h"
#include "grid.h"
#include "dynamics.h"
#include "list.h"

static double min(const double x, const double y, const double z) {
	return (x < y)
			? (x < z ? x : z)
			: (y < z ? y : z);
}

// define functions declared in decs.h
//void scale_vec(vec3 *v, double factor) {
//	v->x *= factor;	v->y *= factor;	v->z *= factor;
//}

double time = 0.;

int main(int argc, char *argv[]) {
	int i;  // loop index var

	// read as inputs in the future
	const int nSteps = ceil(t_end/dt);
	//int nx = 100, ny = 100, nz = 100;
	const int elec_per_cell = 4;
	// print output every output_freq iterations of the main loop
	const int output_freq = nSteps/10;
	// lower left coordinate and dimensions defining the rectangle where particles begin in the simulation
	const vec3 origin = {.46*x_max, .46*y_max, .46*z_max};
	const vec3 dims = {.08*x_max, .08*y_max, .08*z_max};

	printf("initializing grid and particles\n");

	tree ***base_grid = grid_init();
	// add particles to the specified trees
	init_particles(base_grid, origin, dims, elec_per_cell);

	printf("finished initializing. beginning simulation\n");

	for(i = 0; i < nSteps; ++i) {
		time = i*dt;
		printf("pushing particles\n");
		push_particles(base_grid);
		printf("updating grid\n");
		update_grid(base_grid);  // add the laser, etc.
		if (i % output_freq == 0) {
			printf("outputting grid\n");
			output_grid((i/output_freq), (nSteps/output_freq), base_grid);
		}
	}

	// print final state unless it was already output on the last iteration of the above loop
	if ((i-1) % output_freq != 0)
		output_grid((i/output_freq), (nSteps/output_freq), base_grid);
	 
	printf("simulation finished\n");

	//list_free(particles);
	//cleanup(grid_cells);

	return 0;
}
