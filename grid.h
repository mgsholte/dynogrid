#ifndef GRID_H
#define GRID_H

#include "decs.h"

// return the initialized grid. make a grid with 'nx' points along x and 'ny' points along y
grid_point** init_grid(int nx, int ny);

// return the list of particles. put 'part_per_cell' particles in each grid
particle* init_particles(vec2 ul, vec2 lr, int part_per_cell);

void update_grid(grid_point **grid_points);

// add the laser contribution to the fields at the current time step
void add_laser(grid_point **grid_points);

void push_particles(grid_point **grid_points, particle *particles);

// print the fields and particles at iteration i
void output_grid(int i, grid_point **grid_points, particle *particles);

#endif //GRID_H
