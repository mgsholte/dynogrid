#ifndef GRID_H
#define GRID_H

#include "decs.h"
#include "list.h"

// return the initialized grid. make a grid with 'nx' points along x and 'ny' points along y
grid_point** init_grid(int nx, int ny);

// return the list of particles. put 'part_per_cell' particles in each grid within the rectangle with upper left corner 'ul' and lower right corner 'lr'
List init_particles(vec2 ul, vec2 lr, int part_per_cell);

void update_grid(grid_point **grid_points);

// add the laser contribution to the fields at the current time step
void add_laser(grid_point **grid_points);

void push_particles(grid_point **grid_points, List particles);

// print the fields and particles at iteration i
void output_grid(int itNum, grid_point **grid_points, int nx, int ny, List particles);

#endif //GRID_H
