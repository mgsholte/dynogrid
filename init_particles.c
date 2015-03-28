// inits rectangular region of particles, approx evenly distributed
// upper left corner ul
// lower right corner lr
particle* init_particles(vec2 ul, vec2 lr, int part_per_cell) {
	int nCells = (ul.y - lr.y) * (lr.x - ul.x);
	List particles = list_init();
	for (i = 0; i < nCells; i++) {
		// add protons
		for (j = 0; j < part_per_cell/2; j++) {
			particle p = (particle) malloc (sizeof(particle));
			
		 	// check against struct of particle
		}
		// add electrons
		for (j = 0; j < (part_per_cell - part_per_cell/2); j++) {
			