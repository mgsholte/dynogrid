#ifndef DECS_H
#define DECS_H

//Definitions go here to avoid clutter elsewhere

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
	vec2 pos;  // position
	vec3 p;  // momentum
    const double mass, charge, weight;  // weight is the number of actual particles this instance represents
	//particle *next;  // treat the particles as a list. this points to the next particle in the list
} particle;

typedef struct {
	vec3 E, B; // electric, magnetic fields
	// vec3 J; // current density
	// double rho; // average charge density at grid point
} grid_point;

const double dt = 1e-18; // the length of a single iteration in seconds
double time; // changes every iteration

const double c = 3e8;  // speed of light in m/s

// width and height of a cell. aka distance between adjacent grid points
//const double dx = 3e-10; // defined as c*dt so laser moves 1 dx per 1 dt
//const double dy = 3e-10;
const double dx = c*dt; // defined as c*dt so laser moves 1 dx per 1 dt
const double dy = c*dt;

const double PROTON_WEIGHT = 5.0;
const double ELECTRON_WEIGHT = 1.0;
const int BASE_PROTON_MASS = 1.672e-27; //kg
const int BASE_ELECTRON_MASS = 9.109e-31; //kg
const int BASE_PROTON_CHARGE = 1.602e-19; //C
const int BASE_ELECTRON_CHARGE = -1*PROTON_CHARGE; //C

#endif //DECS_H
