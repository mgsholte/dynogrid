#include <stdlib.h>
#include <math.h>
#include <stdio.h>

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

// grid_points and particles are the stuff to print
// suffix is the suffix used in naming the output files
void output_grid_impl(int itNum, int numFiles, grid_point ***grid_points, List particles, const char suffix[]) {
	int suffix_len = strlen(suffix);

	// the # of chars needed to represent the biggest iteration # as a string. all iter #s will be padded to this value
	int width_iter_num = round_i(log10(numFiles)+1);
	// the file name. reused for all files so it needs to be large enough to hold the longest name. +1 at end holds string terminator
    char fname[width_iter_num+7+suffix_len+1];
	FILE *pfile;

	static bool shouldInitHeader = true;
    // create header file only once
	if (shouldWriteHeader) {
		shouldWriteHeader = false;

		// set name to header and open
		sprintf(fname, "params.%s", suffix);
		pfile = fopen(fname, "w");

		// test to ensure that the file was successfully created
		testPFile(pfile, fname);

		fprintf(pfile, "params.data\n");
		fprintf(pfile, "nx=%d,ny=%d,nz=%d\n", nx, ny, nz);
		fprintf(pfile, "dx=%lg,dy=%lg,dz=%lg\n", dx, dy, dz);
		fprintf(pfile, "numFiles=%d\n", numFiles);
		fprintf(pfile, "numParticles=%d\n", numParticles);
		fclose(pfile);
	}

	// GRID OUTPUT
	// set file name to grid file and open
	sprintf(fname, "%d_grid.%s", itNum, suffix);
	pfile = fopen(fname, "w");

    // test to ensure that the file was actually created and exists: 
	testPFile(pfile, fname);

    int x,y,z;
	double magE, magB;
	// print |E|, |B| for each grid point
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
    fclose(pfile);

	// PARTICLE OUTPUT
	// set file name to particle file and open
	sprintf(fname, "%d_ptcls.%s", itNum, suffix);
	pfile = fopen(fname, "w");

    // test to ensure that the file was actually created and exists: 
	testPFile(pfile, fname);

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
