#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "list.h"
#include "decs.h"
#include "grid.h"


void output_data2D(int itNum, grid_cell ***grid_cells, List particles) {
	printf("ERROR: output_data2D not implemented\n");
	0/0; // method not supposed to be implemented
}

void output_data3D(int itNum, int numFiles, grid_cell ***grid_cells, List particles) {
	char *f_extension = "data";
    int numParticles = list_length(particles);
    // create a file to output the grid data to:
	if (itNum == 0) {
		char params[256];
		sprintf(params, "params.%s", f_extension);

		FILE *params_file = fopen(params, "w"); // write only
		// test to ensure that the file was actually created and exists: 
		if (params_file == NULL){
			printf("1. Error! Could not create file\n"); 
			exit(-1); // must include stdlib.h 
		}//end if

		fprintf(params_file, "PARAMS\n");
		fprintf(params_file, "nx=%d,ny=%d,nz=%d\n", nx, ny, nz);
		fprintf(params_file, "dx=%lg,dy=%lg,dz=%lg\n", dx, dy, dz);
		fprintf(params_file, "numFiles=%d\n", numFiles);
        fprintf(params_file, "numParticles=%d\n", numParticles);
		fclose(params_file);
	}

    // create the file name to output for the grid point data:
    char fname_grid[256];
	sprintf(fname_grid, "%d_grid.%s", itNum, f_extension);

    // create a file to output the grid data to:
    FILE *grid_file = fopen(fname_grid, "w"); // write only
    // test to ensure that the file was actually created and exists: 
    if (grid_file == NULL){
        printf("3. Error! Could not create file\n"); 
        exit(-1); // must include stdlib.h 
    }//end if

    int y;
    int x;
    int z;
    fprintf(grid_file, "GRIDPOINTS\n");
    fprintf(grid_file, "itNum=%d\n", itNum);
    fprintf(grid_file, "Time=%lg\n", time);
    fprintf(grid_file, "TimeStep=%lg\n", dt);
    fprintf(grid_file, "GridSize:nx=%d,ny=%d,nz=%d\n", nx, ny, nz);
	//TODO: update this logic to output grid cells rather than grid points
	fprintf(grid_file, "not yet implemented for adaptive grid");
	/*
	for(x = 0; x < nx; x++){
		for(y = 0; y < ny; y++){
            for(z = 0; z < nz; z++){
                //calculate the L2 norm of the E field vector:
                double E = sqrt(pow((grid_cells[x][y][z].E).x, 2.0) + pow((grid_cells[x][y][z].E).y, 2.0) + pow((grid_cells[x][y][z].E).z, 2.0));
                //calculate the L2 norm of the B field vector:
                double B = sqrt(pow((grid_cells[x][y][z].B).x, 2.0) + pow((grid_cells[x][y][z].B).y, 2.0) + pow((grid_cells[x][y][z].B).z, 2.0));
                
                fprintf(grid_file, "%lg,%lg,%lg,%lg,%lg\n", ((double)x * dx), ((double)y * dy), ((double)z * dz), E, B);
                // fprintf(grid_file, "%d,%d,%d,%lg,%lg\n", x, y, z, E, B);
            }//end innermost for
        }//end middle for
    }//end outer for
	*/
    fclose(grid_file);

    // PARTICLE DATA:
    // create the file name to output for the particle data:
    char fname_particles[256];
	sprintf(fname_particles, "%d_particles.%s", itNum, f_extension);

    // create a file to output the particle data to:
    FILE *particles_file = fopen(fname_particles, "w"); // write only
    // test to ensure that the file was actually created and exists: 
    if (particles_file == NULL){
        printf("3. Error! Could not create file\n"); 
        exit(-1); // must include stdlib.h 
    }//end if

    fprintf(particles_file, "PARTICLES\n");
    fprintf(particles_file, "itNum=%d\n", itNum);
    fprintf(particles_file, "Time=%lg\n", time);
    fprintf(particles_file, "TimeStep=%lg\n", dt);
    fprintf(particles_file, "GridSize:nx=%d,ny=%d,nz=%d\n", nx, ny, nz);
    
    int particle_ct = -1;
    list_reset_iter(&particles);
    while(list_has_next(particles)) {
        particle *ptc = list_get_next(&particles);
        // printf("Anything here");
        //calculate the L2 norm of the p vector:
        double p = sqrt(pow((ptc->p).x, 2.0) + pow((ptc->p).y, 2.0) + pow((ptc->p).z, 2.0));
        /*write one line to file per particle in this format:
            ptcl:#,pos_x,pos_y,pos_z,p
        */
        particle_ct++;
		//HACK: fix me later
        //fprintf(grid_file, "ptcl:%d,%lf,%lf,%lf,%lf\n", particle_ct, (ptc->pos).x, (ptc->pos).y, (ptc->pos).z, p);
        fprintf(particles_file, "ptcl:%d,%lg,%lg,%lg,%lg\n", particle_ct, (ptc->pos).x, (ptc->pos).y, (ptc->pos).z, p);
    }//end while
    fclose(particles_file);
}//end output_data3D function
