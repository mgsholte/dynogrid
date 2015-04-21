#ifndef GRID_H
#define GRID_H

#include "decs.h"
#include "list.h"

// return the initialized grid. make a grid with 'nx' points along x, 'ny' points along y, and 'nz' points along z. these 3 parameters are global variables in decs.h
grid_cell*** init_grid();

// return the list of particles. put 'part_per_cell' particles in each grid cell within the rectangular prism with the specified origin and dimensions. the origin is the upper-left-front-most point of the prism
List init_particles(vec3 origin, vec3 dims, int part_per_cell);

// print fields and particles at iteration i
void output_grid(int itNum, int numFiles, grid_cell ***grid_cells, List particles);

#endif //GRID_H
