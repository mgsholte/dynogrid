#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "grid.h"

//laser pulse parameters
#define sigma (15e-6)
#define lambda (5e-6) // 1 micron wavelength
#define pi (3.14159265358979)
#define wavenum (2*pi/lambda) // k, wavenumber
#define freq (wavenum*C)
#define y_mid (y_max/2.0)
#define z_mid (z_max/2.0)

void update_grid(grid_cell ***grid_cells);

void laser(grid_point *grid_p, vec3 loc, double t);
void recursive_laser(grid_cell *cell, double x_spat, double y_spat, double z_spat, int depth, double t);

void push_particles(grid_cell ***grid, List part_list);

#endif //DYNAMICS_H
