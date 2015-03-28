#include "decs.h"
#include "math.h"

// Some globals that I am assuming exist
// dx, dy -- cell size
// dt -- time step

// Particle pusher!!!!
void push_particles(double **grid_points, particle *particles);{
	particle *curr = particles;
	//  Just copying multiplication factors from EPOCH
	idx = 1./dx;
	idy = 1./dy;
	idt = 1./dt;
	dto2 = dt/2.;
	dtco2 = c * dto2;
	dtfac = .5*dt ;//times some weighting factor, that I think is 1 anyways
	third = 1./3.;
	idty = idt*idy;
	idtx = idt*idy;
	idxy = idx*idy;

	//loop over all the particles
	while (curr != NULL){
		ux = curr->p->x * ipart_mc;
		uy = curr->p->y * ipart_mc;
		uz = curr->p->z * ipart_mc;

		//Calculate velocity
		root = dtco2 / sqrt(ux*ux + uy*uy + uz*uz + 1.);

		//Move half timestep
		curr->pos->x += ux * root;
		curr->pos->y += uy * root;

		//Do interpolation to find e and b here.


		// Do the Borris rotation
		

		//Move full timestep
		
		

		//move on to the next one
		curr = curr->next
	}
