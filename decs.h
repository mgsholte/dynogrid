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
	vec3 p, pos;  // momentum, position
    double mass, charge, weight;  // weight is the number of actual particles this instance represents
	//particle *next;  // treat the particles as a list. this points to the next particle in the list
} particle;

typedef struct {
	vec3 E, B; // electric, magnetic fields
	// vec3 J; // current density
	// double rho; // average charge density at grid point
} grid_point;

#define c (3e8)

// the length of a single iteration in seconds
#define dt (1e-18)

double time; // changes every iteration

#define dx (c*dt)
#define dy (c*dt)

#define round(x) ((int) (x+0.5))

#define PROTON_WEIGHT (5.0)
#define ELECTRON_WEIGHT (1.0)
#define BASE_PROTON_MASS (1.672e-27)
#define BASE_ELECTRON_MASS (9.109e-31)
#define BASE_PROTON_CHARGE (1.602e-19)
#define BASE_ELECTRON_CHARGE (-1.0*BASE_PROTON_CHARGE)
        
#endif //DECS_H
