#ifndef DYNAMICS_H
#define DYNAMICS_H

void update_grid(grid_point **grid_points);

// add the laser contribution to the fields at the current time step
void add_laser(grid_point **grid_points);

void push_particles(grid_point **grid_points, List particles);

#endif //DYNAMICS_H
