#include "tree.h"

tree tree_init(vec3 loc) {
	TreeNode *root = (TreeNode*) calloc(sizeof(TreeNode));
	root->points[0] = (grid_point*) calloc(sizeof(grid_point)); // the point starts with E = B = 0
	
	return (tree) { root, loc, list_init() };
}

// cn is the index of the node in its parent's children array
static void tree_node_apply_fcn(TreeNode *node, int cn, void (*f)(grid_point *), double x, double y, double z, vec3 *h) {
	int i;
	// i=0 point covered by parent already
	for (i = 1; i < 8; ++i) {
		if ((cn&i) == cn) {  // avoids calling f twice on 1 point shared by 2 children
			f(node->points[i], x, y, z);  // apply f to current points
		}
	}
	// now apply to all children if they exist
	if (node->children[0]) {  // if 1 child exitsts, then all 8 exist
		scale_vec(h, 0.5); // spacing is halved every level down
		for (i = 0; i < 8; ++i) {
			//NB: h is half what it should be so need to multiply by 2*h. this
			// is accounted for already in the call below
			tree_node_apply_fcn(node->children[i], i, f, x+((i&4)/2*h.x), y+((i&2)*h.y), z+((i&1)*2*h.z), h);
		}
		scale_vec(h, 2.); // spacing was halved going down the tree, now it needs to be doubled going back up
	}
}

// apply the function f to every point in the tree
void tree_apply_fcn(tree t, void (*f)(grid_point *)) {
	vec3 h = (vec3) { dx, dy, dz };
	tree_node_apply_fcn(t.root, 0, f, loc.x, loc.y, loc.z, &h);
}

static bool need_to_coarsen(TreeNode *cell) {
	int i, j;
	grid_point a, b;
	// compare all pairs of points in the cell
	for(i =0; i < 8; ++i) {
		for(j = 7; j > i; --j) {
			a = *(cell->points[i]);
			b = *(cell->points[j]);
			// if any pair of points in the cell are more than
			// 0.2*threhold different in value, then don't coarsen. 
			// 1*threshold is used to decide when to refine
			if (!(fabs(a.B.x - b.B.x) < 0.20*THRESHOLD_B
				&& fabs(a.B.y - b.B.y) < 0.20*THRESHOLD_B
				&& fabs(a.B.z - b.B.z) < 0.20*THRESHOLD_B
				&& fabs(a.E.x - b.E.x) < 0.20*THRESHOLD_E
				&& fabs(a.E.y - b.E.y) < 0.20*THRESHOLD_E
				&& fabs(a.E.z - b.E.z) < 0.20*THRESHOLD_E)) {
					return false;
			}
		}
	}
	return true;
}

static bool need_to_refine(TreeNode *cell) {
	int i, j;
	grid_point a, b;
	for(i =0; i < 8; i++) {
		for(j = 7; j > i; j--) {
			a = *(cell->points[i]);
			b = *(cell->points[j]);
			// if any pair of points in the cell are too different
			// in value, then split the cell
			if (fabs(a.B.x - b.B.x) > THRESHOLD_B
				|| fabs(a.B.y - b.B.y) > THRESHOLD_B
				|| fabs(a.B.z - b.B.z) > THRESHOLD_B
				|| fabs(a.E.x - b.E.x) > THRESHOLD_E
				|| fabs(a.E.y - b.E.y) > THRESHOLD_E
				|| fabs(a.E.z - b.E.z) > THRESHOLD_E) {
					return true;
			}
		}
	}
	return false;
}

static void coarsen(TreeNode *cell) {
	int i, j;
	if (!cell->children[0]) {
		return;
	}
	// 1st, coarsen your children so you have no grandchildren
	for(i = 0; i < 8; i++) {
		coarsen(cell->children[i]);
	}
	// remove your children
	for(i = 0; i < 8; i++) {
		for(j = 0; j < 8; j++) {
			if(i != j && ((j&i) == i)) { //don't free points i==j bc the parent still needs them. also, other children may have already freed one of your points so don't double free
				//free each child's 7 points that are no longer needed:
				free(cell->children[i]->points[j]);
			}
			// free childrens' child pointers
			free(cell->children[i]->children[j]);
		}
		// free your child pointers
		free(cell->children[i]);
		cell->children[i] = NULL;
	}
}

static void tree_node_update(TreeNode *node) {
	// points have little E,B difference at this level so it should be even less on the finer grid. remove all children and coarsen to this level
	if (need_to_coarsen(node)) {
		coarsen(node);
	} else if (node->children[0]) { // find the leaf nodes to see if we need to refine
		int i;
		for (i = 0; i < 8; ++i) {
			tree_node_update(node->children[i]);
		}
	} else if (need_to_refine(node)) {  // found a leaf node
		refine(node);
	}
}

void tree_update(tree t) {
	tree_node_update(t.root);
}

void refine(TreeNode* cell, double x, double y, double z, vec3 *h) {
	int i, j, k, m;
	// allocate the 8 new children cells
	tree *children_cells[8];
	for (i = 0; i < 8; ++i) {
		children_cells[i] = (TreeNode*) calloc( sizeof(TreeNode) );
	}

	// allocate all 19 new points. then have the new child cells reference these 19 new points and the original 8 points from their parent
	grid_point* children_points[3][3][3];
	for (i = 0; i < 8; ++i) {
		// consider 2x2x2 grid super cell. there will be 3 grid_points in each direction b/c 2 adjacent cells share the points on their boundary
		// the elem children_points[i][j][k] holds the point at x=i, y=j, z=k
		children_points[(j&1)*2][(j&2)][(j&4)/2] = cell->points[j];
	}
	// malloc points not pointing to parent points
	for (j = 0; j < 3; ++j) {
		for (k = 0; k < 3; ++k) {
			for (m = 0; m < 3; ++m) {
				// selects only points not pointing to parent points
				if ( j==1 || k==1 || m==1 ) {
					children_points[j][k][m] = (grid_point*) calloc( sizeof(grid_point) );
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
