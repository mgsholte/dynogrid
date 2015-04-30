#ifndef TREE_H
#define TREE_H

#include "decs.h"
#include "list.h"

typedef struct TreeNode {
	grid_point *points[8];  // the 8 vertices of the cell the node represents
	struct TreeNode *children[8];  // the 8 sub-cells created by halving the cell in each dimension
} TreeNode;

// an octree data structure which represents a single cell at the coarsest level of the simulation. each node holds the verticies of the cell at its level of refinement
typedef struct {
	TreeNode *root;
	vec3 loc; // the (x,y,z) location of the lower-left-front point of the root node
	List particles; // the particles within this cell
} tree;

// create a tree only initializing 1 (0-th) of the 8 grid_points in the cell
tree tree_init();

// apply the function f to every point in the tree. f is a fcn that
// accepts a pointer to the grid point and 3 doubles (x,y,z) giving its location
void tree_apply_fcn(tree t, void (*f)(grid_point *, vec3 loc));

// refine/coarsen the tree if necessary
void tree_update(tree t, tree neighbors[3]);

#endif //TREE_H
