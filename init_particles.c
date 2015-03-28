// inits rectangular region of particles, approx evenly distributed
// upper left corner ul
// lower right corner lr
particle* init_particles(vec2 ul, vec2 lr, int part_per_cell) {
	int nCells = (ul.y - lr.y) * (lr.x - ul.x);
	particle* particles;
	for (i = 0; i < nParticles; i++) {
		 // check against struct of particle