#include "decs.h"
#include "grid.h"
#include "list.h"

void output_data2D(int itNum, grid_point ***grid_points, List particles);
void output_grid_impl(int itNum, int numFiles, grid_point ***grid_points, List particles, const char suffix[]);
