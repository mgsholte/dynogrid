#include <stdlib.h>
#include <math.h>

#include "dynamics.h"

// changes E and B at the given point and time
void laser(grid_point *grid_p, double x, double y, double z, double t) {
	double B = B0*cos(freq*t - wavenum*x);
	double E = E0*cos(freq*t - wavenum*x);
	// N is the gaussian distribution factor for 3D
	double pulse_mid = C*t;
	double N = exp(-(pow(x-pulse_mid,2)+pow(y-y_mid,2)+pow(z-z_mid,2))/(2*pow(sigma,2)));
	
	grid_p->E = (vec3) { 0, N*E, 0 };
	grid_p->B = (vec3) { 0, 0, N*B };
}


void update_grid(grid_cell ***grid_cells) {
	int x,y,z;
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny; y++) {
			for (z = 0; z < nz; z++) {
				update_grid_cell(grid_cells[x][y][z]);
				laser(grid_cells[x][y][z].points[0], x*dx, y*dy, z*dz, time);
			}
		}
	}
}
	
void update_grid_cell(grid_cell* cell){
	coarsen(cell);
	refine(cell);
}
	

/*	A grid_cell should only coarsen if it has exactly 1 level of decendants
	(i.e. has children, but does not have grandchildrend or greatgranchildren, etc.
	AND it meets the coarsening criteria based on its E and B fields */
bool coarsen(grid_cell* cell){
	// BASE CASE:
	if(cell->children == NULL){
		/*	coarsening always happens one level up, so if you don't have any children
			you should never coarsen, thus return false */
		return false;
	}
	// RECURSIVE STEP:
	else if(cell->children != NULL){
		// check to see if there are any decendants beyond immediate children:
		bool have_grandchildren = false;
		//TODO: make a bool8 type that encodes 8 bools in a single char. make functions for manipulating said type
		bool childrens_responses[8];
		bool childs_response;
		int child_num;
		for(child_num = 0; child_num < 8; child_num++){
			childs_response = coarsen(cell->children[child_num]);
			if(childs_response == true){
				have_grandchildren = true;
			}
			childrens_responses[child_num] = childs_response; 
		}
		if(have_grandchildren == true){
			return false;
		} else { // I have children, but no grandchildren
			if(need_to_coarsen(cell) == true){
				execute_coarsen(cell);
				/*	now I'm the smallest, since I just executed the coarsen,
					so I return false to keep the recursion chain going */
				return false;
			} else{
				/* 	I have children but no grandchildren, but no need to coarsen,
					so I return true to break the recursion chain so that my
					parent knows not to coarsen */
				return true;
			}
		}
	} else{
		printf("ERROR! This should never happen\n"); //TODO: remove after debugging
	}
}//end coarsen function

bool need_to_coarsen(grid_cell* cell){
	// TODO: check by points E and B fields to see if I should coarsen...

}

void execute_coarsen(grid_cell* cell){
	// TODO: actually execute the coarsening
}

bool refine(grid_cell *cell){
	
}
	


