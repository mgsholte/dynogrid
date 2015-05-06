#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "grid.h"
#include "tree.h"

//laser pulse parameters
#define sigma (15e-6)
#define lambda (5e-6) // 1 micron wavelength
#define pi (3.14159265358979)
#define wavenum (2*pi/lambda) // k, wavenumber
#define freq (wavenum*C)
#define y_mid (y_max/2.0)
#define z_mid (z_max/2.0)

void update_grid(tree ****base_grid);

void laser(grid_point *grid_p, double x, double y, double z);

void push_particles(tree ****base_grid);

#endif //DYNAMICS_H
