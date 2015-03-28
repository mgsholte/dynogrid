#ifndef GRID_H
#define GRID_H

#include "decs.h"

// return the initialized grid. make a grid with 'nx' points along x and 'ny' points along y
double** init_grid(int nx, int ny);

// return the list of particles. put 'part_per_cell' particles in each grid
particle* init_particles(int part_per_cell);

void update_grid(double **grid_points);

// add the laser contribution to the fields at the current time step
void add_laser(double **grid_points);

void push_particles(double **grid_points, particle *particles);

// print the fields and particles at iteration i
void output_grid(int i, double **grid_points, particle *particles);

#endif //GRID_H
