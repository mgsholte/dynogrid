#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "grid.h"

void update_grid(grid_point ***grid_points);

void laser(grid_point *grid_p, double x, double y, double z, double t);

void push_particles(grid_point ***grid, List part_list);

#endif //DYNAMICS_H
