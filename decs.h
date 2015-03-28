#ifndef DECS_H
#define DECS_H

//Definitions go here to avoid clutter elsewhere

typedef struct {
	double x,y;
} vec2;

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

double dt = 1e-18; // the time step
double time; // changes every iteration

#endif //DECS_H
