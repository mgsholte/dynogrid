#ifndef DECS_H
#define DECS_H

#include <mpi.h>

#include "vector.h"

//Definitions go here to avoid clutter elsewhere

typedef enum { false, true } bool;

// Define the particle structure
typedef struct {
	vec3 pos, p;  // position, momentum
    double mass, charge, weight;  // weight is the number of actual particles this instance represents
} particle;

typedef struct {
	vec3 E, B; // electric, magnetic fields
	// vec3 J; // current density
	// double rho; // average charge density at grid point
} grid_point;

// Declare some globals
extern int imin, imax, jmin, jmax, kmin, kmax;	//Processor minimum indicies
extern int pid;	//Processor ID
extern double pxmin, pymin, pzmin;	//Processor minimum x, y, and z

extern MPI_Datatype mpi_vec3, mpi_particle, mpi_tree, mpi_tree_node, mpi_grid_point;

#define C (3e8)

double time; // changes every iteration

#define t_end (3.3333e-13)

#define x_max (1e-4)
#define y_max (1e-4)
#define z_max (1e-4)

// number of cells in each direction; nx+1 is number of grid points
#define nx (30) 
#define ny (30)
#define nz (30)

#define dx (x_max/nx)
#define dy (y_max/ny)
#define dz (z_max/nz)

// Use the smallest of dx, dy, or dz!!!
#define dt (dz/C)

#define round_i(x) ((int) (x+0.5))
//inline int round_i(double x) {
//	return (int) (x+0.5);
//}

#define PROTON_WEIGHT (5.0)
#define ELECTRON_WEIGHT (1.0)
#define BASE_PROTON_MASS (1.672e-27)
#define BASE_ELECTRON_MASS (9.109e-31)
#define BASE_PROTON_CHARGE (1.602e-19)
#define BASE_ELECTRON_CHARGE (-1.0*BASE_PROTON_CHARGE)

#define E0 (8.68e13) // SI units, used Intensity=10^21 W/(cm)^2
#define B0 (2.895e5) // used B0=E0/c
#define THRESHOLD_E (0.1*E0) // threshold for refining E
#define THRESHOLD_B (0.1*B0) // threshold for refining B

#endif //DECS_H
