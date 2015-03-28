// inits rectangular region of particles, approx evenly distributed
// top left corner (topx, topy)
// bottom right corner (botx, boty)
particle* init_particles(int topx, int topy, int botx, int boty, int part_per_cell) {
	int nParticles = (topy - boty) * (botx - topx) * part_per_cell;
	particle* particles;
	particles = (particle*) malloc (nParticles(sizeof(particle)));
	for (i = 0; i < nParticles; i++) {
		particles[i] = // check against struct of particle