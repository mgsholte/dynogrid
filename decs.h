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
	vec2 pos;  // position
	vec3 p;  // momentum
    double mass, charge, weight;  // weight is the number of actual particles this instance represents
	//particle *next;  // treat the particles as a list. this points to the next particle in the list
} particle;

typedef struct {
	vec3 E, B; // electric, magnetic fields
	// vec3 J; // current density
	// double rho; // average charge density at grid point
} grid_point;

#define c 3e8

// the length of a single iteration in seconds
#define dt 1e-18

double time; // changes every iteration

const double dx, dy;

const double PROTON_WEIGHT;
const double ELECTRON_WEIGHT;
const int BASE_PROTON_MASS;
const int BASE_ELECTRON_MASS;
const int BASE_PROTON_CHARGE;
const int BASE_ELECTRON_CHARGE;

#endif //DECS_H
