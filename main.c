#include <stdio.h>
#include <math.h>

#include <mpi.h>
// #include <omp.h>

#include "decs.h"
#include "grid.h"
#include "dynamics.h"
#include "list.h"
#include "defs.h"

static double min(const double x, const double y, const double z) {
	return (x < y)
			? (x < z ? x : z)
			: (y < z ? y : z);
}

double time = 0.;

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	// Definitions for globals
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	part_per_cell = 5;

	int i;  // loop index var
	int x_divs, y_divs, z_divs, i_size, j_size, k_size, numProcs;
	sscanf (argv[1],"%d",&x_divs);
	sscanf (argv[2],"%d",&y_divs);
	sscanf (argv[3],"%d",&z_divs);
	printf("x_divs is: %d\n", x_divs);
	printf("y_divs is: %d\n", y_divs);
	printf("z_divs is: %d\n", z_divs);
	
	//TODO: calculate i_size, j_size, and k_size based on x_divs, y_divs, and z_divs:
	//....


	int MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
	if(numProcs != x_divs*y_divs*z_divs){
		printf("ERROR! numProcs != x_divs*y_divs*z_divs!\nMoron!!!\n");
		return -1;
	}//end if


	// read as inputs in the future
	const int nSteps = ceil(t_end/dt);
	//int nx = 100, ny = 100, nz = 100;
	part_per_cell = 5; //now a global constant
	// print output every output_freq iterations of the main loop
	const int output_freq = nSteps/10;
	// upper left coordinate and dimensions defining the rectangle where particles begin in the simulation
	const vec3 ulf = {.46*x_max, .46*y_max, .46*z_max};
	const vec3 dims = {.08*x_max, .08*y_max, .08*z_max};

	printf("initializing grid and particles\n");

	grid_cell ****grid_cells = init_grid(i_size, j_size, k_size, x_divs, y_divs, z_divs);
	List particles = init_particles(ulf, dims, part_per_cell);

	printf("finished initializing. beginning simulation\n");

	for(i = 0; i < nSteps; ++i) {
		time = i*dt;
		printf("pushing particles\n");
		push_particles(grid_cells, particles);
		printf("updating grid\n");
		update_grid(grid_cells);  // add the laser, etc.
		if (i % output_freq == 0) {
			printf("outputting grid\n");
			output_grid((i/output_freq), (nSteps/output_freq), grid_cells, particles);
		}
	}

	// print final state unless it was already output on the last iteration of the above loop
	if ((i-1) % output_freq != 0)
		output_grid((i/output_freq), (nSteps/output_freq), grid_cells, particles);
	 
	printf("simulation finished\n");

	list_free(particles);
	cleanup(grid_cells);

	MPI_Finalize();

	return 0;
}
