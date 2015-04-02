#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "grid.h"

void update_grid(grid_point **grid_points);

void push_particles(grid_point **grid_points, List particles);

#endif //DYNAMICS_H
