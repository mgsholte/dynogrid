#define "grid.h"

// inits rectangular region of particles, approx evenly distributed
// upper left corner ul
// lower right corner lr
List init_particles(vec2 ul, vec2 lr, int part_per_cell) {
	List particles = list_init();
	for (int i = ul.y; i < lr.y; i++) {
		for (int j = lr.x; j < ul.x; j++) {
			// add protons
			for (int k = 0; k < part_per_cell/2; k++) {
				particle* p = (particle *)malloc(sizeof(particle));
				double x = rand_float(0, dx) + j;
				double y = rand_float(0, dy) + i;
				p = { {x, y},		//position
					  {0, 0, 0},	//momentum
					  BASE_PROTON_MASS * PROTON_WEIGHT		//mass
					  BASE_PROTON_CHARGE * PROTON_WEIGHT	//charge
					  PROTON_WEIGHT		//weight
					};
				list_add(particles, p);
			}
			// add electrons
			for (int k = 0; k < (part_per_cell - part_per_cell/2); k++) {
				particle *p = (particle *)malloc(sizeof(particle));
				double x = rand_float(0, dx) + j;
				double y = rand_float(0, dy) + i;
				p = { {x, y},		//position
					  {0, 0, 0},	//momentum
					  BASE_ELECTRON_MASS * ELECTRON_WEIGHT		//mass
					  BASE_ELECTRON_CHARGE * ELECTRON_WEIGHT	//charge
					  ELECTRON_WEIGHT		//weight
					};
				list_add(particles, p);
			}
		}
	}
	return particles;
}
