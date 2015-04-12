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

#define c (3e8)

double time; // changes every iteration

#define x_max (1e-4)
#define y_max (1e-4)
#define z_max (1e-4)

#define nx (100)
#define ny (100)
#define nz (100)

#define dx (x_max/nx)
#define dy (y_max/ny)
#define dz (z_max/nz)

// Use the smallest of dx, dy, or dz!!!
#define dt (dz/c)

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

//laser pulse parameters
#define sigma (5)
#define lambda (10e-6) // 10 micron wavelength
#define pi (3.14159265358979)
#define wavenum (2*pi/lambda) // k, wavenumber
#define freq (wavenum*c)
#define E0 (8.68e13) // SI units, used Intensity=10^21 W/(cm)^2
#define B0 (2.895e5) // used B0=E0/c
#define y_mid (y_max/2)
#define z_mid (z_max/2)

        
#endif //DECS_H
