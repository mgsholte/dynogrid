#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "grid.h"
#include "list.h"
#include "vector.h"
#include "tree.h"

// only way I could see to make output_one_point work with tree_apply_function
static FILE *pfile = NULL;

static inline void testPFile(FILE *pfile, char *fname) {
	if (pfile == NULL) {
		printf("Error! Could not create/open file: %s\n", fname);
		exit(1); // must include stdlib.h 
	}
}

// print x,y,z coord of the point followed by |E|,|B| at the point
// extra_args: pass a FILE* to the file where the point should be written
void output_one_point(grid_point *point, double x, double y, double z) {
	double E = vec3_norm(point->E),
		B = vec3_norm(point->B);
	fprintf(pfile, "%lg,%lg,%lg,%lg,%lg\n", x, y, z, E, B);
}

// the base_grid has the trees which hold the points and particles to print
// suffix is the suffix used in naming the output files
void output_grid_impl(int itNum, int numFiles, tree ****base_grid, const char suffix[]) {
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
		fprintf(pfile, "nx=%d,ny=%d,nz=%d\n", nx, ny, nz);
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

    int ix,iy,iz;
	// double magE, magB;
	// print |E|, |B| for each grid point
	fprintf(pfile, "x, y, z, |E|, |B|\n");
	for(i = imin+1; i < imax-1; ++i) {
		for(j = jmin+1; j < jmax-1; ++j) {
            for(k = kmin+1; k < kmax-1; ++k) {
            	if (base_grid[i][j][k] != NULL) {
					// tell each cell to print all of the points inside it for which it is responsible
            		// TODO: make sure we aren't double outputting (via diff procs) or missing cells at the end of the grid
					tree_apply_fcn(base_grid[i][j][k], &output_one_point);
				}
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
	//fprintf(pfile, "%d\n", list_length(particles));
	for (ix = 0; ix < nx; ++ix) {
		for (iy = 0; iy < ny; ++iy) {
            for (iz = 0; iz < nz; ++iz) {
				// print x,y,z,|p| for each particle in the cell
				List particles = base_grid[ix][iy][iz].particles;
				list_reset_iter(&particles);
				while(list_has_next(particles)) {
					particle *ptc = list_get_next(&particles);
					fprintf(pfile, "%lg,%lg,%lg,%lg\n", (ptc->pos).x, (ptc->pos).y, (ptc->pos).z, vec3_norm(ptc->p));
				}
			}
		}
	}
    fclose(pfile);
}

