#include "decs.h"
#include <stdlib.h>

void output_grid(int nx, int ny, int itNum, double **grid_points, particle *particles){
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
    fprintf(grid_file, "Time=%d\n", time);
    fprintf(grid_file, "TimeStep=%d\n", dt);
    fprintf(grid_file, "GridSize=%d\n", time);
    for(i = 0; i < ny; i++){
    	for(j = 0; j < nx; j++){
    		/*write one line to file per grid point in this format:
			coordinates:i,j,Ex=lf,Ey=lf,Ez=lf,Bx=lf,By=lf,Bz=lf
			*/
			fprintf(grid_file, "coordinates=%d,%d,Ex=%lf,Ey=%lf,Ez=%lf,Bx=%lf,By=%lf,Bz=%lf\n", i, j, (grid_points[i][j].E).x, (grid_points[i][j].E).y, (grid_points[i][j].E).z, (grid_points[i][j].B).x, (grid_points[i][j].B).y, (grid_points[i][j].B).z); // write to file 



    	}//end inner for

    }//end outer for





}//end output_grid function
          