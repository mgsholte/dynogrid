#include <stdlib.h>
#include <math.h>

#include "dynamics.h"

static const double inv_sigma_sq_2 = 1./(2*sigma*sigma);

// changes E and B at the given point and time
void laser(grid_point *point, double x, double y, double z, double t) {
	double B = B0*cos(freq*t - wavenum*x);
	double E = E0*cos(freq*t - wavenum*x);
	// N is the gaussian distribution factor for 3D
	double pulse_mid = C*t;
	double N = exp(-(pow(x-pulse_mid,2)+pow(y-y_mid,2)+pow(z-z_mid,2))*inv_sigma_sq_2);
	
	grid_p->E = (vec3) { 0, N*E, 0 };
	grid_p->B = (vec3) { 0, 0, N*B };
}

void recursive_laser(grid_cell *cell, double x_spat, double y_spat, double z_spat, int depth, double t) {
	//BASE CASE:
	int cn; //cn stands for child number

	if(cell->children == NULL){
		//cell has no children, so apply laser to each of the cell's gridpoints:
		for(cn = 0; cn < 8; cn++){
			laser(cell->points[cn], x_spat+(cn&1)*dx/pow(2.0,depth),
									y_spat+((cn&2)/2)*dy/pow(2.0,depth),
									z_spat+((cn&4)/4)*dz/pow(2.0,depth), time); // ??? (cn&1)/4 or (cn&4)/4
									// z_spat+((cn&1)/4)*dz/pow(2.0,depth), time);
		}//end for
		return;
	}
	//RECURSIVE STEP:
	else{
		//recursively call recursive_laser on each child cell:
		for(cn = 0; cn < 8; cn++){
			recursive_laser(cell->children[cn], x_spat+(cn&1)*dx/pow(2.0,depth+1),
										   		y_spat+((cn&2)/2)*dy/pow(2.0,depth+1),
										   		z_spat+((cn&4)/4)*dz/pow(2.0,depth+1), depth+1, time);
										   		// z_spat+((cn&1)/4)*dz/pow(2.0,depth+1), depth+1, time);
		}//end for
	}//end if else
}//end recursive_laser function

static void update_grid_cell(grid_cell* cell, int x, int y, int z) {
	vec3 h = (vec3) {dx,dy,dz};
	if(!refine(cell, x*dx, y*dy, z*dz, &h))
		coarsen(cell);
}//end update_grid_cell function

void update_grid(grid_cell ***grid_cells) {
	int x,y,z;
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny; y++) {
			for (z = 0; z < nz; z++) {
				recursive_laser(&grid_cells[x][y][z], x*dx, y*dy, z*dz, 0, time);
				update_grid_cell(&grid_cells[x][y][z], x, y, z);
			}
		}
	}
}
