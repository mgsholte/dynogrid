#include "dynamics.h"

// Currently the laser pulse is step-function-like in both x and y directions with a wavelength
//  of 4*dx and E,B at magnitudes of 10 (in SI, I don't know if this is valid).
// Improvements will likely be necessary, but this should compile.
void update_grid(grid_point ***grid_points) {
	// could make some of these inputs later
	int xCenter = c * time;
	int yCenter = ny/2;
	int zCenter = nz/2;
	double beam_hw = 2; //half width of laser in y and z
	double pulse_hw = 15; //pulse half width: "amplitude" in x
	
	double EMax = 10; // CAUTION: no idea on the validity of these values
	double BMax = 10;
	
	int x,y,z;
	for (y = yCenter - beam_hw; y < yCenter + beam_hw + 1; y++) {
		for (z = zCenter - beam_hw; z < zCenter + beam_hw + 1; z++) {
			// create a (not very physical) plane wave
			for (x = xCenter - pulse_hw; x < xCenter + pulse_hw + 1; x++) {
				if (x >= 0 && x < x_max) {
					if ((xCenter-x)%4 == 0) { // wavelength is 4*dx
						(grid_points[x][y][z]).E = (vec3) { 0, EMax, 0 };
						(grid_points[x][y][z]).B = (vec3) { 0, 0, BMax };
					} else if ((xCenter-x)%4 == 2) {
						(grid_points[x][y][z]).E = (vec3) { 0, -1*EMax, 0 };
						(grid_points[x][y][z]).B = (vec3) { 0, 0, -1*BMax };
					} else {
						(grid_points[x][y][z]).E = (vec3) { 0, 0, 0 };
						(grid_points[x][y][z]).B = (vec3) { 0, 0, 0 };
					}
				}
			}
			// clean up from earlier time steps
			for (x = xCenter - pulse_hw - 2; x < xCenter - pulse_hw; x++) {
				if (x >= 0 && x < x_max) {
					(grid_points[x][y][z]).E = (vec3) { 0, 0, 0 };
					(grid_points[x][y][z]).B = (vec3) { 0, 0, 0 };
				}

			}
		}
	}
}
