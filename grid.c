#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "decs.h"
#include "grid.h"
#include "output.h"
#include "dynamics.h"

// return a random double value in the range [low,high]
static inline double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

static void scale_vec(vec3 *v, double factor) {
	v->x *= factor;	v->y *= factor;	v->z *= factor;
}

// inits all grid points to 0 in E and B
// isize, etc. are lengths of real cell dimensions. ghost and NULL cells are added in this function
grid_cell**** init_grid(int isize, int jsize, int ksize, int x_divs, int y_divs, int z_divs) {
	// global vars (local indices)
	imin = isize/2 - 1; // -1 for ghost on left
	jmin = jsize/2 - 1;
	kmin = ksize/2 - 1;
	imax = imin + isize + 2; // +2 for ghosts on left and right
	jmax = jmin + jsize + 2;
	kmax = kmin + ksize + 2;
	// global vars (global spatial positions)
	px_min = ((double)(pid % x_divs))/x_divs * x_max - dx;
	py_min = ((double)((pid - pid % x_divs)/x_divs % y_divs))/y_divs * y_max - dy;
	pz_min = ((double)((pid - pid % x_divs - pid % (x_divs * y_divs))/(x_divs * y_divs)))/z_divs * z_max - dz;
	
	// local vars
	int wi = 2*isize; // padding is 50% of isize (etc.) on each side
	int wj = 2*jsize;
	int wk = 2*ksize;
	
	
	grid_cell ****grid_cells = (grid_cell****) malloc( (wi+1) * sizeof(grid_cell***) ); // allocate an array of pointers to rows-depthwise
	int i, j, k, n;
	int di, dj, dk, owner_id;
	for (i = 0; i < wi; ++i) {
		grid_cells[i] = (grid_cell***) malloc( (wj+1) * sizeof(grid_cell**) );  // allocate the row
		for (j = 0; j < wj; ++j) {
			grid_cells[i][j] = (grid_cell**) malloc( (wk+1) * sizeof(grid_cell*) );  // allocate the row
			// each index goes from (ghost cell on left) to (initialization ghost cell on right) which is one past (ghost cell on right) for the
			// purpose of initializing all points. each cell makes just one point, so an extra layer is needed at the end
			for (k = 0; k < wk; ++k) {
				// malloc for real, ghost, and 'init ghost' cells
				if (i >= imin && i < imax+1 &&
					j >= jmin && j < jmax+1 &&
					k >= kmin && k < kmax+1) {
					
					grid_cells[i][j][k] = (grid_cell*) malloc( sizeof(grid_cell) );
					
					// initialize only the upper-left-forward grid_point (points[0] represents (z,y,x)==000, points[3] repr. (z,y,x)==011)
					grid_cells[i][j][k]->points[0] = malloc( sizeof(grid_point) );
					grid_cells[i][j][k]->points[0]->E = (vec3) {0,0,0};
					grid_cells[i][j][k]->points[0]->B = (vec3) {0,0,0};
					grid_cells[i][j][k]->children = NULL;
					
					// I think this will work, assuming true->1 and false->0
					// 26 possible neighbors could own each ghost cell, but their pids can be constructed from true/false statements
					di = (i == imax-1) - (i == imin); // +1 or -1
					dj = (j == jmax-1) - (j == jmin);
					dk = (k == kmax-1) - (k == kmin);
					owner_id = pid;
					owner_id += di;
					owner_id += dj * x_divs;
					owner_id += dk * x_divs * y_divs;
					
					// find bad cases where there is no proper owner_id, i.e. the above algorithm got a bad answer b/c ghost is out of simulation bounds
					if ((di == 1  &&  px_min + dx * (isize+1.5) >= x_max) || 	// 1.5 is to prevent rounding issues, even though 1 should be enough
						(dj == 1  &&  py_min + dy * (isize+1.5) >= y_max) ||
						(dk == 1  &&  pz_min + dz * (isize+1.5) >= z_max) ||
						(di == -1  &&  px_min + dx * 0.5 <= 0) ||
						(dj == -1  &&  py_min + dy * 0.5 <= 0) ||
						(dk == -1  &&  pz_min + dz * 0.5 <= 0)) {
						
						owner_id = -1;
					}
					
					grid_cells[i][j][k]->owner = owner_id;
				}
			}
		}
	}
	// make other 7 grid_point pointers of each grid cell point to grid points allocated by neighbors. except for fake 'init ghost' cells on right/bottom/back boundaries
	for (i = imin; i < imax; ++i) {
		for (j = jmin; j < jmax; ++j) {
			for (k = kmin; k < kmax; ++k) {
				// n should be thought of as binary (n for "neighbors")
				for (n = 1; n < 8; ++n) {
					grid_cells[i][j][k].points[n] = grid_cells[i+(n&1)][j+(n&2)/2][k+(n&4)/4].points[0];
				}
			}
		}
	}
	return grid_cells;
}

// populates particles randomly within each grid cell inside the bounds of the specified prism. particles are evenly distributed across the cells
// if the origin and dims of the prism don't align exactly to grid points, they will be rounded to the nearest ones
void init_particles(grid_cell**** grid_cells, vec3 origin, vec3 dims, int part_per_cell) {
	//ensure that the origin of the tungsten block aligns evenly with increments of dx, dy, and dz:
	double tungsten_xmin = dx*(round_i(origin.x/dx)), tungsten_xmax = dx*(round_i(dims.x/dx)) + tungsten_xmin;
	double tungsten_ymin = dy*(round_i(origin.y/dy)), tungsten_ymax = dy*(round_i(dims.y/dy)) + tungsten_ymin;
	double tungsten_zmin = dz*(round_i(origin.z/dz)), tungsten_zmax = dz*(round_i(dims.z/dz)) + tungsten_zmin;
	
	int i, j, k, ppc; //i corresonds to x; j to y; k to z; ppc for particles per cell
	double cell_xmin; cell_xmax, cell_ymin, cell_ymax, cell_zmin, cell_zmax;
	particle *p; // pointer to the next particle to add...why declare this outside the loops???
	double x,y,z; // the coords of the next particle to add
	for(i = imin; i < imax; i++){
		cell_xmin = px_min + ((i-imin) * dx);
		cell_xmax = cell_xmin + dx;
		for(j = jmin; j < jmax; j++){
			cell_ymin = py_min + ((j-jmin) * dy);
			cell_ymax = cell_ymin + dy;
			for(k = kmin; k < kmax; k++){
				cell_zmin = pz_min + ((k-kmin) * dz);
				cell_zmax = cell_zmin + dz;
				
				List particles = list_init(); //particles list for the current time step
				List next_particles = list_init(); //particles list for the next time step...will be initiated to empty
				
				if(	(cell_xmax >= tungsten_xmin && cell_xmin <= tungsten_xmax) &&
					(cell_ymax >= tungsten_ymin && cell_ymin <= tungsten_ymax) &&
					(cell_zmax >= tungsten_zmin && cell_zmin <= tungsten_zmax) &&
					(i > imin && i < imax) &&
					(j > jmin && j < jmax) &&
					(k > kmin && k < kmax)
					){	
					
					//then this cell lies within the block of tungsten (and NOT a ghost cell), so it gets particles:
					for (ppc = 0; ppc < part_per_cell/2; ppc++) {
						// add a proton
						p = (particle*) malloc(sizeof(particle));

						x = rand_float(0, dx) + cell_xmin;
						y = rand_float(0, dy) + cell_ymin;
						z = rand_float(0, dz) + cell_zmin;
						*p = (particle) { 
							{x, y, z},  //position
							{0, 0, 0},  //momentum
							BASE_PROTON_MASS * PROTON_WEIGHT,   //mass
							BASE_PROTON_CHARGE * PROTON_WEIGHT, //charge
							PROTON_WEIGHT  //weight
							};

						list_add(particles, p);

						// add an electron
						p = (particle*) malloc(sizeof(particle));

						x = rand_float(0, dx) + cell_xmin;
						y = rand_float(0, dy) + cell_ymin;
						z = rand_float(0, dz) + cell_zmin;
						*p = (particle) { 
							{x, y, z},  //position
							{0, 0, 0},  //momentum
							BASE_ELECTRON_MASS * ELECTRON_WEIGHT,   //mass
							BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT, //charge
							ELECTRON_WEIGHT  //weight
							};

						list_add(particles, p);
						list_reset_iter(&particles); //do we need to do this or no??? (Joel asking, not sure if this is needed)
					}//end ppc for loop
				}//end if
				//adding initial list of particles to cell:
				(grid_cells[i][j][k])->part_list = particles;
				//adding blank, but initialized list of particles for next time step:
				(grid_cells[i][j][k])->next_list = next_particles;
			}//end k for loop
		}//end j for loop
	}// end i for loop

	// //old, sequential version...keeping around for now just in case we need to compare/revert anything from the above:
	// int iy, ix, iz, k;
	// int ix_min = round_i(origin.x/dx), ix_max = round_i(dims.x/dx) + ix_min;
	// int iy_min = round_i(origin.y/dy), iy_max = round_i(dims.y/dy) + iy_min;
	// int iz_min = round_i(origin.z/dz), iz_max = round_i(dims.z/dz) + iz_min;
	// double x,y,z; // the coords of the next particle to add
	// particle *p; // pointer to the next particle to add
	// for (ix = ix_min; ix < ix_max; ++ix) {
	// 	for (iy = iy_min; iy < iy_max; ++iy) {
	// 		for (iz = iz_min; iz < iz_max; ++iz) {
	// 			// add particles in this cell
	// 			for (k = 0; k < part_per_cell/2; ++k) {
	// 				// add a proton
	// 				p = (particle*) malloc(sizeof(particle));

	// 				x = rand_float(0, dx) + ix*dx;
	// 				y = rand_float(0, dy) + iy*dy;
	// 				z = rand_float(0, dz) + iz*dz;
	// 				*p = (particle) { 
	// 					{x, y, z},  //position
	// 					{0, 0, 0},  //momentum
	// 					BASE_PROTON_MASS * PROTON_WEIGHT,   //mass
	// 					BASE_PROTON_CHARGE * PROTON_WEIGHT, //charge
	// 					PROTON_WEIGHT  //weight
	// 					};

	// 				list_add(particles, p);

	// 				// add an electron
	// 				p = (particle*) malloc(sizeof(particle));

	// 				x = rand_float(0, dx) + ix*dx;
	// 				y = rand_float(0, dy) + iy*dy;
	// 				z = rand_float(0, dz) + iz*dz;
	// 				*p = (particle) { 
	// 					{x, y, z},  //position
	// 					{0, 0, 0},  //momentum
	// 					BASE_ELECTRON_MASS * ELECTRON_WEIGHT,   //mass
	// 					BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT, //charge
	// 					ELECTRON_WEIGHT  //weight
	// 					};

	// 				list_add(particles, p);
	// 			}
	// 		}
	// 	}
	// }
	// list_reset_iter(&particles);
	// return particles;
}//end init_particles() function

static bool need_to_coarsen(grid_cell* cell) {
	int i, j;
	grid_point a, b;
	for(i =0; i < 8; i++) {
		for(j = 7; j > i; j--) {
			a = *(cell->points[i]);
			b = *(cell->points[j]);
			if (fabs(a.B.x - b.B.x) < 0.20*THRESHOLD_B
				&& fabs(a.B.y - b.B.y) < 0.20*THRESHOLD_B
				&& fabs(a.B.z - b.B.z) < 0.20*THRESHOLD_B
				&& fabs(a.E.x - b.E.x) < 0.20*THRESHOLD_E
				&& fabs(a.E.y - b.E.y) < 0.20*THRESHOLD_E
				&& fabs(a.E.z - b.E.z) < 0.20*THRESHOLD_E) {
					return true;
				}
		}//end inner for
	}//end outer for
	return false;
}//end need_to_refine function

static void execute_coarsen(grid_cell* cell) {
	int i, j;
	for(i = 0; i < 8; i++) {
		for(j = 0; j < 8; j++) {
			if(i != j && ((j&i) == i)) { //don't free points i==j bc the parent still needs them. also, other children may have already freed one of your points so don't double free
				//free each child's 7 points that are no longer needed:
				free(cell->children[i]->points[j]);
			}
			free(cell->children[i]->children);
		}//end inner for
		//free the child cell:
		free(cell->children[i]);
	}//end outer for
	//free parent's children list:
	free(cell->children);
	cell->children = NULL;
}//end execute_coarsen function

bool coarsen(grid_cell* cell) {
	// called on a leaf. leaves can't be coarsened
	if(cell->children == NULL) {
		/*	coarsening always happens one level up, so if you don't have any children
			you should never coarsen, thus return false */
		return false;
	}
	// called on an internal node. check for grand-children before coarsening
	else {
		// check to see if there are any decendants beyond immediate children:
		bool have_grandchildren = false;
		int child_num;
		for(child_num = 0; child_num < 8; child_num++){
			have_grandchildren = coarsen(cell->children[child_num])
				? true
				: have_grandchildren;
		}

		if(have_grandchildren) {
			// don't coarsen a node that has grandchildren
			return true;
		} else if(need_to_coarsen(cell)) {
			execute_coarsen(cell);
			/*	now I'm the smallest, since I just executed the coarsen,
				so I return false to keep the recursion chain going */
			return false;
		} else {
			/* 	I have children but no grandchildren, but no need to coarsen,
				so I return true to break the recursion chain so that my
				parent knows not to coarsen */
			return true;
		}
	}
}//end coarsen function

static bool need_to_refine(grid_cell* cell) {
	int i, j;
	grid_point a, b;
	for(i =0; i < 8; i++) {
		for(j = 7; j > i; j--) {
			a = *(cell->points[i]);
			b = *(cell->points[j]);
			if (fabs(a.B.x - b.B.x) > THRESHOLD_B
				|| fabs(a.B.y - b.B.y) > THRESHOLD_B
				|| fabs(a.B.z - b.B.z) > THRESHOLD_B
				|| fabs(a.E.x - b.E.x) > THRESHOLD_E
				|| fabs(a.E.y - b.E.y) > THRESHOLD_E
				|| fabs(a.E.z - b.E.z) > THRESHOLD_E) {
					return true;
				}
		}//end inner for
	}//end outer for
	return false;
}//end need_to_refine function

void execute_refine(grid_cell* cell, double x_spat, double y_spat, double z_spat, vec3 *h) {
	// create 8 new children
	grid_cell **children_cells;
	children_cells = (grid_cell**) malloc( 8*sizeof(grid_cell*) );
	int i, j, k, m;
	for(i  = 0; i < 8; ++i){
		children_cells[i] = (grid_cell*) malloc( sizeof(grid_cell) );
		children_cells[i]->children = NULL;
	}

	// creating all needed points (27 total), referencing 8 back to parent's points,
	//  allocating the rest, then having all the child cells reference these 27 "master" points
	grid_point* children_points[3][3][3];
	for (j = 0; j < 8; ++j){
		// consider 2x2x2 grid super cell. there will be 3 grid_points in each direction.
		// the elem children_points[i][j][k] holds the point at z=i, y=j, x=k
		children_points[(j&4)/2][(j&2)][(j&1)*2] = cell->points[j];
	}
	// malloc points not pointing to parent points
	for (j = 0; j < 3; ++j) {
		for (k = 0; k < 3; ++k) {
			for (m = 0; m < 3; ++m) {
				// selects only points not pointing to parent points
				if ( j==1 || k==1 || m==1 ) {
					children_points[j][k][m] = (grid_point*) malloc( sizeof(grid_point) );
					// NOTE: if adding more than laser to a grid point, add that here
					laser(children_points[j][k][m], x_spat + j*h->x,
													y_spat + k*h->y,
													z_spat + m*h->z, time);
				}
			}
		}
	}

	for (i = 0; i < 8; ++i){
	    children_cells[i]->points[0] = children_points[0+(i&4)/4][0+(i&2)/2][0+(i&1)];
	    children_cells[i]->points[1] = children_points[0+(i&4)/4][0+(i&2)/2][1+(i&1)];
	    children_cells[i]->points[2] = children_points[0+(i&4)/4][1+(i&2)/2][0+(i&1)];
	    children_cells[i]->points[3] = children_points[0+(i&4)/4][1+(i&2)/2][1+(i&1)];
	    children_cells[i]->points[4] = children_points[1+(i&4)/4][0+(i&2)/2][0+(i&1)];
	    children_cells[i]->points[5] = children_points[1+(i&4)/4][0+(i&2)/2][1+(i&1)];
	    children_cells[i]->points[6] = children_points[1+(i&4)/4][1+(i&2)/2][0+(i&1)];
	    children_cells[i]->points[7] = children_points[1+(i&4)/4][1+(i&2)/2][1+(i&1)];
	}

	// finally, make the cell the parent of the newly created children
	cell->children = children_cells;
}//end execute_refine function

/*	A grid_cell should only coarsen if it has exactly 1 level of decendants
	(i.e. has children, but does not have grandchildrend or greatgranchildren, etc.
	AND it meets the coarsening criteria based on its E and B fields */
/*	A grid_cell should only refine if it has no children AND it meets
	the refining criteria based on its E and B fields */
bool refine(grid_cell* cell, double x_spat, double y_spat, double z_spat, vec3 *h) {
	// called refine on a leaf. refine if refinement condition is met
	if(cell->children == NULL) {
		if(need_to_refine(cell)) {
			execute_refine(cell, x_spat, y_spat, z_spat, h);
			return true;
		} else {
			return false;
		}
	}
	// cell either started as an internal node or has become one via refinement. pass the refine call to the children to see if they need to split
	int cn;
	scale_vec(h,0.5);
	for(cn = 0; cn < 8; cn++) {
		refine(cell->children[cn], x_spat+(cn&1)*h->x,
								   y_spat+((cn&2)/2)*h->y,
								   z_spat+((cn&4)/4)*h->z, h);
								   // z_spat+((cn&1)/4)*dz/pow(2.0,depth+1), depth+1);
	}//end for 
	scale_vec(h,2.);
	return false;
}
	
void output_grid(int itNum, int numFiles, grid_cell ****grid_cells, List particles) {
	output_grid_impl(itNum, numFiles, grid_cells, particles, "data");
	//TODO: it should be true that itNum == time/dt. maybe we don't need to pass the itNum variable as an argument
	// int itNum = round_i(time/dt);
}

//TODO: change these loop bounds to whatever suits the changing grid chunk size. not just imin to imax, right? there's uninitialized cells.
void cleanup(grid_cell ****grid_cells) {
	int x,y,z,i;
	grid_cell* cell;
	for (x = 0; x <= nx; x++) {
		for (y = 0; y <= ny; y++) {
			for (z = 0; z <= nz; z++) {
				// coarsen completely (coarsest should have no children after)
				cell = grid_cells[x][y][z];
				recursive_execute_coarsen(cell);
				
				// only free what was malloc'd per cell in init_grid
				free(cell->points[0]);
			}
			// frees cells and cell.children address
			free(grid_cells[x][y]);
		}
		free(grid_cells[x]);
	}
	free(grid_cells);
}
void recursive_execute_coarsen(grid_cell* cell) {
	// BASE CASE or if at coarsest cell size from the start
	if (cell->children == NULL) {
		return;
	}
	else {
		int child_num;
		for(child_num = 0; child_num < 8; child_num++) {
			recursive_execute_coarsen(cell->children[child_num]);
		}
		execute_coarsen(cell);
	}
}
	
	
	
	
