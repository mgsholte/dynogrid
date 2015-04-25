#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "decs.h"
#include "grid.h"
#include "output.h"

// return a random double value in the range [low,high]
static inline double rand_float( double low, double high ) {
	return ( (double)rand() * (high - low) ) / (double)	RAND_MAX + low;
}

static void scale_vec(vec3 *v, double factor) {
	v->x *= factor;
	v->y *= factor;
	v->z *= factor;
}

// inits all grid points to 0 in E and B
grid_cell*** init_grid() {
	grid_cell ***grid_cells = (grid_cell***) malloc( (nx+1) * sizeof(grid_cell**) ); // allocate an array of pointers to rows-depthwise
	int i, j, k, n;
	for (i = 0; i <= nx; ++i) {
		grid_cells[i] = (grid_cell**) malloc( (ny+1) * sizeof(grid_cell*) );  // allocate the row
		for (j = 0; j <= ny; ++j) {
			grid_cells[i][j] = (grid_cell*) malloc( (nz+1) * sizeof(grid_cell) );  // allocate the row
			// initialize only the upper-left-forward grid_point (points[0] represents (z,y,x)==000, points[3] repr. (z,y,x)==011)
			for (k = 0; k <= nz; ++k) {
				grid_cells[i][j][k].points[0] = malloc( sizeof(grid_point) );
				grid_cells[i][j][k].points[0]->E = (vec3) {0,0,0};
				grid_cells[i][j][k].points[0]->B = (vec3) {0,0,0};
				grid_cells[i][j][k].children = NULL;
			}
		}
	}
	// make other 7 grid_point pointers of each grid cell point to grid points allocated by neighbors. except for right/upper/back boundaries
	for (i = 0; i < nx; ++i) {
		for (j = 0; j < ny; ++j) {
			for (k = 0; k < nz; ++k) {
				// n should be thought of as binary (n for "neighbors")
				for (n = 1; n < 8; ++n) {
					grid_cells[i][j][k].points[n] = grid_cells[i+(n&1)][j+(n&2)/2][k+(n&4)/4].points[0];
				}
			}
		}
	}
	return grid_cells;
}

// populates particles randomly within each grid cell inside the bounds of the specified prism. particles are evenly distrubuted across the cells
// if the origin and dims of the prism don't algin exactly to grid points, they will be rounded to the nearest ones
List init_particles(vec3 origin, vec3 dims, int part_per_cell) {
	List particles = list_init();
	int row, col, dep, k;
	int x_i_min = round_i(origin.x/dx), x_i_max = round_i(dims.x/dx) + x_i_min;
	int y_i_min = round_i(origin.y/dy), y_i_max = round_i(dims.y/dy) + y_i_min;
	int z_i_min = round_i(origin.z/dz), z_i_max = round_i(dims.z/dz) + z_i_min;
	double x,y,z; // the coords of the next particle to add
	particle *p; // pointer to the next particle to add
	for (col = x_i_min; col < x_i_max; ++col) {
		for (row = y_i_min; row < y_i_max; ++row) {
			for (dep = z_i_min; dep < z_i_max; ++dep) {
				// add particles in this cell
				for (k = 0; k < part_per_cell/2; ++k) {
					// add a proton
					p = (particle*) malloc(sizeof(particle));

					x = rand_float(0, dx) + col*dx;
					y = rand_float(0, dy) + row*dy;
					z = rand_float(0, dz) + dep*dz;
					*p = (particle) { 
						{x, y, z},  //position
						{0, 0, 0},  //momentum
						BASE_PROTON_MASS * PROTON_WEIGHT,   //mass
						BASE_PROTON_CHARGE * PROTON_WEIGHT, //charge
						PROTON_WEIGHT  //weight
						};

					list_add(&particles, p);

					// add an electron
					p = (particle*) malloc(sizeof(particle));

					x = rand_float(0, dx) + col*dx;
					y = rand_float(0, dy) + row*dy;
					z = rand_float(0, dz) + dep*dz;
					*p = (particle) { 
						{x, y, z},  //position
						{0, 0, 0},  //momentum
						BASE_ELECTRON_MASS * ELECTRON_WEIGHT,   //mass
						BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT, //charge
						ELECTRON_WEIGHT  //weight
						};

					list_add(&particles, p);
				}
			}
		}
	}
	list_reset_iter(&particles);
	return particles;
}

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
			if(i != j && (j&i == i)) { //don't free points i==j bc the parent still needs them. also, other children may have already freed one of your points so don't double free
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
		bool childs_response;
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
				//printf("children_points[%d][%d][%d] = %p\n", j,k,m, children_points[j][k][m]);
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
		for (j = 0; j < 8; ++j) {
			if (children_cells[i]->points[j] == NULL) {
				printf("cell %d, point %d is null\n", i, j);
				printf("goodbye\n");
			}
		}
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
	
void output_grid(int itNum, int numFiles, grid_cell ***grid_cells, List particles) {
	output_grid_impl(itNum, numFiles, grid_cells, particles, "data");
	//TODO: it should be true that itNum == time/dt. maybe we don't need to pass the itNum variable as an argument
	// int itNum = round_i(time/dt);
}
