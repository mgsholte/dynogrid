#include "decs.h"
#include <stdlib.h>
#include "list.h"

void output_grid(int nx, int ny, int itNum, double **grid_points, List particles){
	// GRIDPOINT DATA:
	// create the file name to output for the grid point data:
	char fname_grid[] = ((char)itNum);
	strcat(fname_grid, "_grid.txt");

	// create a file to output the grid data to:
	FILE *grid_file = fopen(fname_grid, "w"); // write only
    // test to ensure that the file was actually created and exists: 
    if (grid_file == NULL){
    	printf("Error! Could not create file\n"); 
        exit(-1); // must include stdlib.h 
    }//end if

    int i;
    int j;
    fprintf(grid_file, "GRIDPOINTS\n");
    fprintf(grid_file, "itNum=%d\n", itNum);
    fprintf(grid_file, "Time=%lf\n", time);
    fprintf(grid_file, "TimeStep=%lf\n", dt);
    fprintf(grid_file, "GridSize:nx=%d,ny=%d\n", nx, ny);
    for(i = 0; i < ny; i++){
    	for(j = 0; j < nx; j++){
    		/*write one line to file per grid point in this format:
			coordinates:i,j,Ex=lf,Ey=lf,Ez=lf,Bx=lf,By=lf,Bz=lf
			*/
			fprintf(grid_file, "coordinates=%d,%d,Ex=%lf,Ey=%lf,Ez=%lf,Bx=%lf,By=%lf,Bz=%lf\n", i, j, (grid_points[i][j].E).x, (grid_points[i][j].E).y, (grid_points[i][j].E).z, (grid_points[i][j].B).x, (grid_points[i][j].B).y, (grid_points[i][j].B).z);
    	}//end inner for
    }//end outer for
    fclose(grid_file);

    // PARTICLE DATA:
    // create the file name to output for the particle data:
	char fname_particles[] = ((char)itNum);
	strcat(fname_particles, "_particles.txt");

	// create a file to output the particle data to:
	FILE *particles_file = fopen(fname_particles, "w"); // write only
    // test to ensure that the file was actually created and exists: 
    if (particles_file == NULL){
    	printf("Error! Could not create file\n"); 
        exit(-1); // must include stdlib.h 
    }//end if

    fprintf(grid_file, "PARTICLES\n");
    fprintf(grid_file, "itNum=%d\n", itNum);
    fprintf(grid_file, "Time=%lf\n", time);
    fprintf(grid_file, "TimeStep=%lf\n", dt);
    fprintf(grid_file, "GridSize:nx=%d,ny=%d\n", nx, ny);
    
    list_reset_iter(particles);
    while(list_has_next(particles)){
    	/*write one line to file per grid point in this format:
			coordinates:i,j,Ex=lf,Ey=lf,Ez=lf,Bx=lf,By=lf,Bz=lf
			*/
		particle *ptc = (particle*) list_get_next(particles);
		fprintf(grid_file, "pos_x=%lf,pos_y=%lf,p_x=%lf,p_y=%lf,p_z=%lf\n", (ptc->pos).x, (ptc->pos).y, (ptc->pos).x, (ptc->p).x, (ptc->p).y, (ptc->p).z);
	}//end while
    fclose(particles_file);
}//end output_grid function
          