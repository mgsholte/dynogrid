#include "dynamics.h"
#include <math.h>

/*
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
*/

// New implementation:
// - iterate over all x? this automatically cleans up
// - will I have explicit limits to y and z in update_grid? or will it be laser()'s
// responsibility? laser could start with "if y or z outside (range) return;". this
// would make laser fully self-contained. there seems no reason to put any laser pulse
// info outside of laser, it would just be messier.

// void update_grid(grid_cell ***grid_cells) {
// 	int x,y,z;
// 	for (x = 0; x < nx; x++) {
// 		for (y = 0; y < ny; y++) {
// 			for (z = 0; z < nz; z++) {
// 				laser(&(grid_cells[x][y][z].points[0]), x*dx, y*dy, z*dz, time);
// 			}
// 		}
// 	}
// }

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

void recursive_laser(grid_cell *cell, double x_spat, double y_spat, double z_spat, int depth, double t){
	//BASE CASE:
	int cn; //cn for child number
	if(cell->children == NULL){
		//cell has no children, so apply laser:
		for(cn = 0; cn < 8; cn++){
			laser(cell->points[cn],)
		}
	}
	//RECURSIVE STEP:
	else{
		//recursively call recursive_laser on each child cell:
		for(cn = 0; cn < 8; cn++){
			recursive_laser(cell->children[cn], x_spat+(cn&1)*dx/pow(2.0,depth+1),
										   		y_spat+((cn&2)/2)*dy/pow(2.0,depth+1),
										   		z_spat+((cn&1)/4)*dz/pow(2.0,depth+1), depth+1, time);
		}//end for
	}//end if else



}//end recursive_laser function


void update_grid(grid_cell ***grid_cells) {
	int x,y,z;
	for (x = 0; x < nx; x++) {
		for (y = 0; y < ny; y++) {
			for (z = 0; z < nz; z++) {
				recursive_laser(&(grid_cells[x][y][z].points[0]), x*dx, y*dy, z*dz, time);
				update_grid_cell(grid_cells[x][y][z], x, y, z);
			}
		}
	}
}
	
void update_grid_cell(grid_cell* cell, int x, int y, int z){
	refine(cell, x*dx, y*dy, z*dz, 0);
	coarsen(cell);
}//end update_grid_cell function
	

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
		// bool chidrens_responses[8];
		bool childs_response;
		int child_num;
		for(child_num = 0; child_num < 8; child_num++){
			childs_response = coarsen(cell->children[child_num]);
			if(childs_response == true){
				have_grandchildren = true;
			}
			// childrens_responses[child_num] = childs_response; 
		}
		if(have_grandchildren == true){
			return false;
		} else { // I have children, but no grandchildren
			if(need_to_refine(cell) != true){
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


void execute_coarsen(grid_cell* cell){
	int i, j;
	for(i = 0; i < 8; i++){
		for(j = 0; j < 8; j++){
			if(i == j){
				//don't free these points bc the parent still needs them...
			}else{
				//free each child's 7 points that are no longer needed:
				free(cell->children[i].points[j]);
			}
		}//end inner for
		//free the child cell:
		free(cell.children[i]);
	}//end outer for
	//free parent's children list:
	free(cell.children);
	cell.children == NULL;
}//end execute_coarsen function

/*	A grid_cell should only refine if it has no children AND it meets
	the refining criteria based on its E and B fields */
void refine(grid_cell* cell, double x_spat, double y_spat, double z_spat, int depth){
	// BASE CASE:
	if(cell->children == NULL){
		/*	check to see if refining is needed...if yes, refine and then 
			recursively call refine on each child cell */
		if(need_to_refine(cell)){
			execute_refine(cell, x_spat, y_spat, z_spat, depth);
			//cell now has children after refining, so refine needs to be called on each child:
			int cn; //child num
			for(cn = 0; cn < 8; cn++){
				refine(cell->children[cn], x_spat+(cn&1)*dx/pow(2.0,depth+1),
										   y_spat+((cn&2)/2)*dy/pow(2.0,depth+1),
										   z_spat+((cn&1)/4)*dz/pow(2.0,depth+1), depth+1);
			}//end for 
		}else{
			return;
		}
	}
	// RECURSIVE STEP:
	else if(cell->children != NULL){
		int cn;
			for(cn = 0; cn < 8; cn++){
				refine(cell->children[cn], x_spat+(cn&1)*dx/pow(2.0,depth+1),
										   y_spat+((cn&2)/2)*dy/pow(2.0,depth+1),
										   z_spat+((cn&1)/4)*dz/pow(2.0,depth+1), depth+1);
			}//end for 
	}else{
		printf("ERROR! This should never happen\n"); //TODO: remove after debugging
	}	
}//end refine function
	
bool need_to_refine(grid_cell* cell){
	int i;
	int j;
	double B_dif_x;
	double B_dif_y;
	double B_dif_z;
	double E_dif_x;
	double E_dif_y;
	double E_dif_z;
	for(i =0; i < 8; i++){
		for(j = 7; j > i; j--){
			B_dif_x = ((cell->points[i])->B).x - ((cell->points[j])->B).x;
			B_dif_y = ((cell->points[i])->B).y - ((cell->points[j])->B).y;
			B_dif_z = ((cell->points[i])->B).z - ((cell->points[j])->B).z;

			E_dif_x = ((cell->points[i])->E).x - ((cell->points[j])->E).x;
			E_dif_y = ((cell->points[i])->E).y - ((cell->points[j])->E).y;
			E_dif_z = ((cell->points[i])->E).z - ((cell->points[j])->E).z;
			if((B_dif_x > THRESHOLD) || (B_dif_y > THRESHOLD) || (B_dif_z > THRESHOLD) || (E_dif_x > THRESHOLD) || (E_dif_y > THRESHOLD) || (E_dif_z > THRESHOLD)){
				return true;
			}
		}//end inner for
	}//end outer for
	return false;
}//end need_to_refine function

void execute_refine(grid_cell* cell, double x_spat, double y_spat, double z_spat, int depth){
	// create 8 new children
	grid_cell *children_cells[8];
	int i, j, k, m;
	for(i  = 0; i < 8; ++i){
		children_cells[i] = (grid_cell*) malloc( sizeof(grid_cell) );
		children_cells[i].children = NULL;
	}
	// creating all needed points (27 total), referencing 8 back to parent's points,
	//  allocating the rest, then having all the child cells reference these 27 "master" points
	grid_point* children_points[3][3][3];
	children_points[0][0][0] = cell.points[0];
	for (j = 0; j < 8; ++j){
		children_points[(j&1)*2][(j&2)][(j&4)/2] = cell.points[j];
	}
	// malloc points not pointing to parent points
	for (j = 0; j < 3; ++j){
		for (k = 0; k < 3; ++k){
			for (m = 0; m < 3; ++m){
				// selects only points not pointing to parent points
				if ( j==1 || k==1 || m==1 ){
					children_points[j][k][m] = (grid_point*) malloc( sizeof(grid_point) );
					// NOTE: if adding more than laser to a grid point, add that here
					laser(children_points[j][k][m], x_spat + j*dx/pow(2,depth+1),
													y_spat + k*dy/pow(2,depth+1),
													z_spat + m*dz/pow(2,depth+1), time);
				}
			}
		}
	}

	for (i = 0; i < 8; ++i){
	    children_cells[i].points[0] = children_points[0+(i&1)][0+(i&2)/2][0+(i&4)/4];
	    children_cells[i].points[1] = children_points[1+(i&1)][0+(i&2)/2][0+(i&4)/4];
	    children_cells[i].points[2] = children_points[0+(i&1)][1+(i&2)/2][0+(i&4)/4];
	    children_cells[i].points[3] = children_points[1+(i&1)][1+(i&2)/2][0+(i&4)/4];
	    children_cells[i].points[4] = children_points[0+(i&1)][0+(i&2)/2][1+(i&4)/4];
	    children_cells[i].points[5] = children_points[1+(i&1)][0+(i&2)/2][1+(i&4)/4];
	    children_cells[i].points[6] = children_points[0+(i&1)][1+(i&2)/2][1+(i&4)/4];
	    children_cells[i].points[7] = children_points[1+(i&1)][1+(i&2)/2][1+(i&4)/4];
	}
}//end execute_refine function

