#include "decs.h"

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
const int BASE_ELECTRON_CHARGE = -1*BASE_PROTON_CHARGE; //C

