#ifndef GRID_H
#define GRID_H

#include "decs.h"
#include "list.h"
#include "tree.h"

//typedef struct {
//	tree ***base;
//} grid;

// return the initialized grid. make a grid with 'nx' points along x, 'ny' points along y, and 'nz' points along z. these 3 parameters are global variables in decs.h
tree**** grid_init(int isize, int jsize, int ksize, int x_divs, int y_divs, int z_divs);

void grid_free(tree ****t);

// Populate each cell with its list of particles. Put 'elec_per_cell' electrons and as many protons in each grid cell within the rectangular prism with the specified origin and dimensions. the origin is the point of the prism with the smallest coords
void init_particles(tree ****base_grid, vec3 origin, vec3 dims, int elec_per_cell);

// print fields and particles at iteration i
void output_grid(int itNum, int numFiles, tree ****base_grid);
//void cleanup(grid_cell ***grid_cells);

vec3 get_loc(int ix, int iy, int iz);

#endif //GRID_H
