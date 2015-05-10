#include <stdlib.h>
#include <math.h>

#include "dynamics.h"

static const double inv_sigma_sq_2 = 1./(2*sigma*sigma);

// changes E and B at the given point and time
bool laser(grid_point *point, double x, double y, double z) {
	double B = B0*cos(freq*time - wavenum*x);
	double E = E0*cos(freq*time - wavenum*x);
	// N is the gaussian distribution factor for 3D
	double pulse_mid = C*time;
	double N = exp(-(pow(x-pulse_mid,2)+pow(y-y_mid,2)+pow(z-z_mid,2))*inv_sigma_sq_2);
	
	point->E = (vec3) { 0, N*E, 0 };
	point->B = (vec3) { 0, 0, N*B };

	return false;
}

//TODO: update for parallel
void grid_update(tree ****base_grid) {
	int ix, iy, iz;
	// set the laser fields at every point
	for (ix = 0; ix < nx; ++ix) {
		for (iy = 0; iy < ny; ++iy) {
			for (iz = 0; iz < nz; ++iz) {
				tree_apply_fcn(base_grid[ix][iy][iz], &laser);
			}
		}
	}
	// now check to see if the grid needs to refine/coarsen
	// each tree is responsible for refining the space between it and the points in the +x/y/z dirs. thus, the last point in each dimension never needs to refine
	for (ix = 0; ix < nx; ++ix) {
		for (iy = 0; iy < ny; ++iy) {
			for (iz = 0; iz < nz; ++iz) {
				tree_update(base_grid[ix][iy][iz]);
			}
		}
	}
}
