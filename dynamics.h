#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "grid.h"

//laser pulse parameters
#define sigma (15e-6)
#define lambda (5e-6) // 1 micron wavelength
#define pi (3.14159265358979)
#define wavenum (2*pi/lambda) // k, wavenumber
#define freq (wavenum*c)
#define E0 (8.68e13) // SI units, used Intensity=10^21 W/(cm)^2
#define B0 (2.895e5) // used B0=E0/c
#define y_mid (y_max/2.0)
#define z_mid (z_max/2.0)

void update_grid(grid_point ***grid_points);

void laser(grid_point *grid_p, double x, double y, double z, double t);

void push_particles(grid_point ***grid, List part_list);

#endif //DYNAMICS_H
