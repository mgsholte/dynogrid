#ifndef GRID_H
#define GRID_H

#include "decs.h"
#include "list.h"

// return the initialized grid. make a grid with 'nx' points along x, 'ny' points along y, and 'nz' points along z. these 3 parameters are global variables in decs.h
grid_cell*** init_grid();


bool coarsen(grid_cell *cell);
bool refine(grid_cell* cell, double x_spat, double y_spat, double z_spat, vec3 *h);

// print fields and particles at iteration i
void output_grid(int itNum, int numFiles, grid_cell ***grid_cells, List particles);
void cleanup(grid_cell ***grid_cells);
#endif //GRID_H
