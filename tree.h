#ifndef TREE_H
#define TREE_H

#include "decs.h"
#include "list.h"

// extract direction encoded in index
#define getX(i) ((i&4)/4)
#define getY(i) ((i&2)/2)
#define getZ(i) (i&1)
#define getIdx(x,y,z) ((x<<2)+(y<<1)+z)

typedef struct TreeNode {
	grid_point *points[8];  // the 8 vertices of the cell the node represents
	struct TreeNode *children[8];  // the 8 sub-cells created by halving the cell in each dimension
} TreeNode;

// an octree data structure which represents a single cell at the coarsest level of the simulation. each node holds the verticies of the cell at its level of refinement
typedef struct {
	TreeNode *root;
	vec3 loc; // the (x,y,z) location of the lower-left-front point of the root node in absolute coords
	List *particles; // the particles within this cell
	List *new_particles; // the particles added to this cell for the next time step
	int owner; // the rank of the processor which owns this cell
	int neighbor_owners[3][3]; // exclusively for tree passing, used to construct new ghost cells (1 new tree = up to 9 new ghosts)
	int nPoints; // the # of points which this cell is responsible for
} tree;

// Simplified tree strcture for use in mpi passing routine...an octree data structure which represents a single cell at the coarsest level of the simulation. each node holds the verticies of the cell at its level of refinement
typedef struct {
	vec3 loc; // the (x,y,z) location of the lower-left-front point of the root node
	int owner; // the rank of the processor which owns this cell
	int nPoints;
	int neighbor_owners[9];
} simple_tree;

// allocate a tree only initializing 1 (0-th) of the 8 grid_points in the cell
// you must link the other 7 grid points yourself.
// all calls to tree_init should be matched by a call to tree_free
tree* tree_init(vec3 loc, int owner);

// every call to tree_init should be matched by a call to this when the tree is no longer needed
void tree_free(tree *t);

// apply the function f to every point in the tree. f is a fcn that
// accepts a pointer to the grid point and 3 doubles (x,y,z) giving its location
// it should return true if you want to continue iterating. return false to abort the iteration
void tree_apply_fcn(tree *t, bool (*f)(grid_point *,double,double,double));

// refine/coarsen the tree if necessary
void tree_update(tree *t);

#endif //TREE_H
