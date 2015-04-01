#include <stdlib.h>

#include "grid.h"

// return a random double value in the range [low,high]
static double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

grid_point** init_grid(int nx, int ny) {
	grid_point **grid = (grid_point**) malloc( ny * sizeof(size_t) ); // allocate an array of pointers to rows
	for (int i = 0; i < ny; ++i) {
		grid[i] = (grid_point*) malloc( nx * sizeof(grid_point) );  // allocate the row
		for (int j = 0; j < nx; j++) {
			// all grid points are initialized with no field
			grid[i][j].E = {0,0,0};
			grid[i][j].B = {0,0,0};
		}
	}
	return grid;
}

// inits rectangular region of particles, approx evenly distributed
// upper left corner ul
// lower right corner lr
List init_particles(vec2 ul, vec2 lr, int part_per_cell) {
	List particles = list_init();
	for (int row = lr.y; row < ul.y; row++) {
		for (int col = ul.x; col < lr.x; col++) {
			// add protons
			for (int k = 0; k < part_per_cell/2; k++) {
				particle* p = (particle*) malloc(sizeof(particle));
				double x = rand_float(0, dx) + j;
				double y = rand_float(0, dy) + i;
				p = { {x, y},		//position
					  {0, 0, 0},	//momentum
					  BASE_PROTON_MASS * PROTON_WEIGHT		//mass
					  BASE_PROTON_CHARGE * PROTON_WEIGHT	//charge
					  PROTON_WEIGHT		//weight
					};
				list_add(particles, p);
			}
			// add electrons
			for (int k = 0; k < (part_per_cell - part_per_cell/2); k++) {
				particle *p = (particle*) malloc(sizeof(particle));
				double x = rand_float(0, dx) + j;
				double y = rand_float(0, dy) + i;
				p = { {x, y},		//position
					  {0, 0, 0},	//momentum
					  BASE_ELECTRON_MASS * ELECTRON_WEIGHT		//mass
					  BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT	//charge
					  ELECTRON_WEIGHT		//weight
					};
				list_add(particles, p);
			}
		}
	}
	return particles;
}

void output_grid(int itNum, grid_point **grid_points, int nx, int ny, List particles) {
	//TODO: it should be true that itNum == time/dt. maybe we don't need to pass the itNum variable as an argument
	//#define round(x) (int) (x+0.5)
	// int itNum = round(time/dt);

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
    	/*write one line to file per particle in this format:
			pos_x=#,pos_y=#,p_x=#,p_y=#,p_z=#
		*/
		particle *ptc = list_get_next(particles);
		fprintf(grid_file, "pos_x=%lf,pos_y=%lf,p_x=%lf,p_y=%lf,p_z=%lf\n", (ptc->pos).x, (ptc->pos).y, (ptc->pos).x, (ptc->p).x, (ptc->p).y, (ptc->p).z);
	}//end while
    fclose(particles_file);
}//end output_grid function
          
