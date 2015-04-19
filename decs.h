#ifndef DECS_H
#define DECS_H

//Definitions go here to avoid clutter elsewhere

typedef enum { false, true } bool;

// a 2D vector
typedef struct {
	double x,y;
} vec2;

// a 3D vector
typedef struct {
	double x,y,z;
} vec3;

// Define the particle structure
typedef struct {
	vec3 pos, p;  // position, momentum
    double mass, charge, weight;  // weight is the number of actual particles this instance represents
	//particle *next;  // treat the particles as a list. this points to the next particle in the list
} particle;

typedef struct {
	vec3 E, B; // electric, magnetic fields
	// vec3 J; // current density
	// double rho; // average charge density at grid point
} grid_point;

typedef struct {
	grid_point *points[8];
	struct grid_cell **children;
} grid_cell;


#define C (3e8)

double time; // changes every iteration

#define t_end (3.3333e-13)

#define x_max (1e-4)
#define y_max (1e-4)
#define z_max (1e-4)

// number of cells in each direction; nx+1 is number of grid points
#define nx (20) 
#define ny (20)
#define nz (20)

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

#endif //DECS_H
