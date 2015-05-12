#include <stdlib.h>

#include "tree.h"
#include "math.h"

// loc = location of the point in the cell with the smallest coordinates
// owner = the rank of the mpi process which owns this cell
tree* tree_init(vec3 loc, int owner) {
	TreeNode *root = (TreeNode*) malloc(sizeof(TreeNode));
	root->points[0] = (grid_point*) malloc(sizeof(grid_point)); // the point starts with E = B = 0
	root->points[0]->E = (vec3) { 0., 0., 0. };
	root->points[0]->B = (vec3) { 0., 0., 0. };
	root->points[1] = NULL; // helpful for identifying recently made tree
	root->children[0] = NULL; // don't need to set others b/c they'll be set later

	tree *ans = (tree*) malloc(sizeof(tree));
	
	*ans = (tree) { root, loc, list_init(), list_init(), owner, 
			{{-1, -1, -1}, //TODO: real init of this var
			 {-1, -1, -1},
			 {-1, -1, -1}}
		 };
	return ans;
}

// fwd decl for tree_free
static void coarsen(TreeNode *node);

void tree_free(tree *t) {
	TreeNode *root = t->root;
	// remove all children from the root
	coarsen(root);

	// remove only the 0-th grid_point since that is what was initialized by tree_init
	free(root->points[0]);
	root->points[0] = NULL;

	// free the root itself
	free(root);
	t->root = NULL;

	// free the 2 particle lists
	list_free(t->particles);
	list_free(t->new_particles);

	// finally, free the tree itself
	free(t);
}

// see notes for tree_apply_fcn
// cn is the index of the node in its parent's children array
static void tree_node_apply_fcn(TreeNode *node, int cn, bool (*f)(grid_point *,double,double,double), double x, double y, double z, vec3 *h) {
	int i;
	//TODO: does this call f on any point twice?
	// i=0 point covered by parent already
	for (i = 1; i < 8; ++i) {
		if ((cn&i) == cn) {  // avoids calling f twice on 1 point shared by 2 children
			if(f(node->points[i], x, y, z))  // apply f to current points
				return; // the function requested to abort iteration
		}
	}
	// now apply to all children if they exist
	if (node->children[0]) {  // if 1 child exitsts, then all 8 exist
		vec3_scale(h, 0.5); // spacing is halved every level down
		for (i = 0; i < 8; ++i) {
			//TODO: optimization test. can this be replaced with 2*getX(i) and produce the same code?
			//NB: h is half what it should be so need to multiply by 2*h. this
			// is accounted for already in the call below
			tree_node_apply_fcn(node->children[i], i, f, x+((i&4)/2*h->x), y+((i&2)*h->y), z+((i&1)*2*h->z), h);
		}
		vec3_scale(h, 2.); // spacing was halved going down the tree, now it needs to be doubled going back up
	}
}

// apply the function f to every point in the tree
void tree_apply_fcn(tree *t, bool (*f)(grid_point *,double,double,double)) {
	vec3 h = (vec3) { dx, dy, dz };
	if(!f(t->root->points[0], t->loc.x, t->loc.y, t->loc.z))
		tree_node_apply_fcn(t->root, 0, f, t->loc.x, t->loc.y, t->loc.z, &h);
}

static bool need_to_coarsen(TreeNode *cell) {
	int i, j;
	grid_point a, b;
	// compare all distinct pairs of points in the cell
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
	// compare all distinct pairs of points in the cell
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
	// can't coarsen a cell with no children
	if (cell->children[0] == NULL) {
		return;
	}
	// 1st, coarsen your children so you have no grandchildren
	for(i = 0; i < 8; i++) {
		coarsen(cell->children[i]);
	}

	// now, remove your children
	for(i = 0; i < 8; i++) {
		TreeNode *pChild = cell->children[i];
		for(j = 0; j < 8; j++) {
			// don't free points i==j bc the parent still needs them. also, other children may have already freed one of your points so don't double free
			if(i != j && ((j&i) == i)) { 
				free(pChild->points[j]);
			}
			pChild->points[j] = NULL;
			//NB: grandchildren already freed by the call to coarsen(cell->children[i])
		}
		// free your child pointers
		free(pChild);
		cell->children[i] = NULL;
	}
}

static void refine(TreeNode* cell, double x, double y, double z, vec3 *h) {
	int i, j, k;
	// allocate the 8 new children cells
	TreeNode *children_cells[8];
	for (i = 0; i < 8; ++i) {
		children_cells[i] = (TreeNode*) malloc(sizeof(TreeNode));
		children_cells[i]->children[0] = NULL;
	}

	// allocate all 19 new points. then have the new child cells reference these 19 new points and the original 8 points from their parent for a total of 27 points in the new 2x2x2 supercell
	// consider 2x2x2 grid super cell. there will be 3 grid_points in each direction b/c 2 adjacent cells share the points on their boundary
	// the elem children_points[i][j][k] holds the point given by the vector x=i, y=j, z=k where i|j|k=2 means move 2 points in that direction, i.e. the verticies of the original, unsplit cell
	grid_point* children_points[3][3][3];
	for (i = 0; i < 3; ++i) {
		bool isI = (i == 1);
		for (j = 0; j < 3; ++j) {
			bool isJ = (j == 1);
			for (k = 0; k < 3; ++k) {
				bool isK = (k == 1);
				// gives the index of the closest parent cell with indicies rounded down
				//TODO: these divide by 2s are really just an extra shift
				#define closest_parent(i,j,k) (((i/2)<<2) + ((j/2)<<1) + (k/2))
				if (isI || isJ || isK) {  // these are the child points
					// TODO: init points
					children_points[i][j][k] = (grid_point*) malloc(sizeof(grid_point));
					grid_point *pChildPoint = children_points[i][j][k];
					pChildPoint->E = (vec3) { 0., 0., 0. };
					pChildPoint->B = (vec3) { 0., 0., 0. };
					// set value at grid point by averaging appropriate parent points (the nearest neighbors)
					unsigned int ip, jp, kp;  // the parent point direction
					int basep = closest_parent(i,j,k);
					for (ip = 0; ip <= isI; ++ip) {
						for (jp = 0; jp <= isJ; ++jp) {
							for (kp = 0; kp <= isK; ++kp) {
								grid_point *pParentPoint = cell->points[basep+getIdx(ip,jp,kp)];
								vec3_add(&(pChildPoint->E), pParentPoint->E);
								vec3_add(&(pChildPoint->B), pParentPoint->B);
							}
						}
					}
					// the number of NN is 2^(# of 1s) where # of 1s counts how many of i,j,k are equal to 1
					int nAveraged = ((1<<isI)<<isJ)<<isK;
					vec3_scale(&(pChildPoint->E), 1./nAveraged);
					vec3_scale(&(pChildPoint->B), 1./nAveraged);
				} else {  // these are the parent points
					children_points[i][j][k] = cell->points[closest_parent(i,j,k)];
				}
				#undef closest_parent
			}
		}
	}

	// tell each new child which points it needs
	for (i = 0; i < 8; ++i) {
		int ix = getX(i),
			iy = getY(i),
			iz = getZ(i);
		for (j = 0; j < 8; ++j) {  // 8 points of the i-th child
			children_cells[i]->points[j] = children_points[ix+getX(j)][iy+getY(j)][iz+getZ(j)];
		}
	}

	// finally, tell the parent the new cells are its children
	for (i = 0; i < 8; ++i) {
		cell->children[i] = children_cells[i];
	}
}

static void tree_node_update(TreeNode *node, double x, double y, double z, vec3 *h) {
	// points have little E,B difference at this level so it should be even less on the finer grid. remove all children and coarsen to this level
	if (need_to_coarsen(node)) {
		coarsen(node);
	} else if (node->children[0]) { // find the leaf nodes to see if we need to refine/coarsen there
		vec3_scale(h, 0.5); // spacing is halved every level down
		int i;
		for (i = 0; i < 8; ++i) {
			tree_node_update(node->children[i], x+getX(i), y+getY(i), z+getZ(i), h);
		vec3_scale(h, 2.); // spacing was halved going down the tree, now it needs to be doubled going back up
		}
	} else if (need_to_refine(node)) {  // found a leaf node, might need to refine
		vec3_scale(h, 0.5); // spacing is halved every level down
		refine(node, x, y, z, h);
		vec3_scale(h, 2.); // spacing was halved going down the tree, now it needs to be doubled going back up
	}
}

void tree_update(tree *t) {
	vec3 h = (vec3) { dx, dy, dz };
	tree_node_update(t->root, t->loc.x, t->loc.y, t->loc.z, &h);
}
