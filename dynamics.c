#include "dynamics.h"

// Currently the laser pulse is step-function-like in both x and y directions with a wavelength
//  of 4*dx and E,B at magnitudes of 10 (in SI, I don't know if this is valid).
// Improvements will likely be necessary, but this should compile.
void update_grid(grid_point **grid_points) {
	// could make some of these inputs later
	// assuming ny = 100 (check with main())
	double yCenter = 50; // ny/2
	double beam_amp = 2; //amplitude in y
	double xCenter = c * time;
	double pulse_hw = 15; //pulse half width: "amplitude" in x
	
	double EMax = 10; // no idea on the validity of these values
	double BMax = 10;
	
	int i; int j;
	for (i = yCenter - beam_amp; i < yCenter + beam_amp; i++) {
		// create a (not very physical) plane wave
		for (j = xCenter - pulse_hw; j < xCenter + pulse_hw; j++) {
			if (j%4 == 0) { // wavelength is 4*dx
				grid_points[i][j].E = { 0, EMax, 0 };
				grid_points[i][j].B = { 0, 0, Bmax };
			} else if (j%4 == 2) {
				grid_points[i][j].E = { 0, -1*EMax, 0 };
				grid_points[i][j].B = { 0, 0, -1*Bmax };
			} else {
				grid_points[i][j].E = { 0, 0, 0 };
				grid_points[i][j].B = { 0, 0, 0 };
			}
		}
		// clean up from earlier time steps
		for (j = xCenter - pulse_hw - 2; j < xCenter - pulse_hw; j++) {
			grid_points[i][j].E = { 0, 0, 0 };
			grid_points[i][j].B = { 0, 0, 0 };
		}
	}
}
