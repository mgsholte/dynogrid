#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "grid.h"

void update_grid(grid_point **grid_points);

void push_particles(grid_point ***grid, List part_list);

#endif //DYNAMICS_H
