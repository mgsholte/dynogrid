#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "grid.h"
#include "list.h"

static inline const double norm(const vec3 v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

static inline void testPFile(FILE *pfile, char *fname) {
	if (pfile == NULL) {
		printf("Error! Could not create/open file: %s\n", fname);
		exit(1); // must include stdlib.h 
	}
}

/*	Recursive function to output the grid_point data for each grid_cell:
	Notes:	(1) only data for leaf grid_cell's needs to be output
			(2) each leaf grid_cell only needs to output its (0,0) coordinate grid_point (I think...)
*/
void output_one_cell(grid_cell* cell, double x_spat, double y_spat, double z_spat, int depth, FILE *pfile){
	//BASE CASE:
	if(cell->children == NULL){
		double E = norm(((cell)->points[0])->E);
		double B = norm(((cell)->points[0])->B);
		fprintf(pfile, "%lg,%lg,%lg,%lg,%lg\n", x_spat, y_spat, z_spat, E, B);
	}
	//RECURSIVE STEP:
	else{
		int cn; //child num
		for(cn = 0; cn < 8; cn++){
			output_one_cell(cell->children[cn], x_spat+(cn&1)*dx/pow(2.0,depth+1),
									   y_spat+((cn&2)/2)*dy/pow(2.0,depth+1),
									   z_spat+((cn&4)/4)*dz/pow(2.0,depth+1), depth+1, pfile);
									   // z_spat+((cn&1)/4)*dz/pow(2.0,depth+1), depth+1);
		}//end for
	}//end else
}//end out_one_coarse_cell() function

// grid_points and particles are the stuff to print
// suffix is the suffix used in naming the output files
void output_grid_impl(int itNum, int numFiles, grid_cell ***grid_cells, List particles, const char suffix[]) {
	int suffix_len = strlen(suffix);

	// the # of chars needed to represent the biggest iteration # as a string. all iter #s will be padded to this value
	int width_iter_num = round_i(log10(numFiles)+1);
	// the file name. reused for all files so it needs to be large enough to hold the longest name. +1 at end holds string terminator
    char fname[width_iter_num+11+suffix_len+1];
	FILE *pfile;

	static bool shouldWriteHeader = true;
    // create header file only once
	if (shouldWriteHeader) {
		shouldWriteHeader = false;

		// set name to header and open
		sprintf(fname, "params.%s", suffix);
		pfile = fopen(fname, "w");

		// test to ensure that the file was successfully created
		testPFile(pfile, fname);

		fprintf(pfile, "params.data\n");
		fprintf(pfile, "nx=%d,ny=%d,nz=%d\n", nx+1, ny+1, nz+1);
		fprintf(pfile, "dx=%lg,dy=%lg,dz=%lg\n", dx, dy, dz);
		fprintf(pfile, "numFiles=%d\n", numFiles);
		fclose(pfile);
	}

	// GRID OUTPUT
	// set file name to grid file and open
	sprintf(fname, "%d_grid.%s", itNum, suffix);
	pfile = fopen(fname, "w");

    // test to ensure that the file was actually created and exists: 
	testPFile(pfile, fname);

    int x,y,z;
	// double magE, magB;
	// print |E|, |B| for each grid point
	fprintf(pfile, "|E|, |B|\n");
	// fprintf(pfile, "not yet implemented for adaptive grid\n");
	/* TODO: make work for grid cells
	for(x = 0; x <= nx; x++) {
		for(y = 0; y <= ny; y++) {
            for(z = 0; z <= nz; z++) {
                //calculate the L2 norms of the fields
				magE = norm(grid_points[x][y][z].E);
				magB = norm(grid_points[x][y][z].B);

                fprintf(pfile, "%lg,%lg\n", magE, magB);
            }
        }
    }
	*/

	for(x = 0; x <= nx; x++) {
		for(y = 0; y <= ny; y++) {
            for(z = 0; z <= nz; z++) {
				output_one_cell(&(grid_cells[x][y][z]), x*dx, y*dy, z*dz, 0, pfile);
            }
        }
    }
	
    fclose(pfile);

	// PARTICLE OUTPUT
	// set file name to particle file and open
	sprintf(fname, "%d_particles.%s", itNum, suffix);
	pfile = fopen(fname, "w");

    // test to ensure that the file was actually created and exists: 
	testPFile(pfile, fname);

	fprintf(pfile, "x, y, z, |p|\n");
	// print # of particles as a header
	fprintf(pfile, "%d\n", list_length(particles));
    list_reset_iter(&particles);
	// print x,y,z,|p| for each particle
    while(list_has_next(particles)) {
        particle *ptc = list_get_next(&particles);
        fprintf(pfile, "%lg,%lg,%lg,%lg\n", (ptc->pos).x, (ptc->pos).y, (ptc->pos).z, norm(ptc->p));
    }
    fclose(pfile);
}

