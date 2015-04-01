#include <math.h>

#include "decs.h"
#include "dynamics.h"

// Some globals that I am assuming exist
// dx, dy -- cell size
// dt -- time step

// Interpolation helper
vec3 interp(vec3 *field00, vec3 *field01, vec3 *field10, vec *field11, double xlf, double yuf){
    vec3 interped;
    int i;
    double norm = 1./(dx*dy);//Make this a global for speed
    double xrf = 1. - xlf;
    double ydf = 1.- yuf;
    // Not sure if these are proper vec3 accesses.  May need to unroll.
    for (i=0;i<3;i++){
        // Bilinear interpolation from Wikipedia.  There are certainly
        // better interpolation methods to use for PIC codes.
        interped[i] = norm * ( field00[i]*xrf*yuf + field10[i]*xlf*yuf + field01*xrf*ydf + field11[i]*xlf*ydf);
    }
}

// Particle pusher!!!!
void push_particles(grid_point **grid_points, particle *particles);{

    // Declare variables!
    particle *curr = particles;
    double ux, uy, uz;
    double root;
    double xl, yu, xlf, yuf;
    vec3 E, B;
    double ux, uy, uz;
    double uxm, uym, uzm;
    double uxp, uyp, uzp;
    double part_mc, ipart_mc;
    //  Just copying multiplication factors from EPOCH
    double idx = 1./dx;
    double idy = 1./dy;
    double idt = 1./dt;
    double dto2 = dt/2.;
    double dtco2 = c * dto2;
    double dtfac = .5*dt ;//times some weighting factor, that I think is 1 anyways
    double third = 1./3.;
    double idty = idt*idy;
    double idtx = idt*idy;
    double idxy = idx*idy;
    
    //loop over all the particles
    while (curr != NULL){
        part_mc = c*curr->mass;
        ipart_mc = 1./part_mc;

        // u is gamma*v, see Birdsall and Langdon sectoin 15-4
	ux = curr->px * ipart_mc;
	uy = curr->py * ipart_mc;
	uz = curr->pz * ipart_mc;
	//Calculate velocity
	root = dtco2 / sqrt(ux*ux + uy*uy + uz*uz + 1.);

	//Move half timestep
	curr->x += ux * root;
	curr->y += uy * root;

	//Do interpolation to find e and b here.
        // x-left and y-up indices
        xl = floor(curr->x * idx);
        yu = floor(curr->y * idy);
        xlf = curr->x % dx;
        yuf = curr->y % dy;
        E = interp(grid[xl][yu].E, grid[xl][yu+1].E, grid[xl+1][yu], grid[xl+1][yu+1].E, xlf, yuf);
        B = interp(grid[xl][yu].B, grid[xl][yu+1].B, grid[xl+1][yu], grid[xl+1][yu+1].B, xlf, yuf);
        
        // Update momenta to u_-, from Birdsall and Langdor
        uxm = ux + cmratio * E.x;
        uym = uy + cmratio * E.y;
        uzm = uz + cmratio * E.z;

        // Do the Borris rotation and update momenta to u_+, from Birdsall and Langdon
        root = cmratio = sqrt(uxm*uxm + uym*uym + uzm*uzm + 1.0);

        taux = B.x*root;
        taux2 *= taux;
        tauy = B.y*root;
        tauy2 *= tauy;
        tauz = B.z*root;
        touz2 *= tauz;

        tau = 1. / (1. + taux2 + tauy2 + tauz2);

        uxp = ((1. + taux2 - tauy2 - tauz2)*uxm + 2.*((taux*tauy + tauz)*uym + (taux*tauz - tauy)*uzm)) * tau;
        uyp = ((1. - taux2 + tauy2 - tauz2)*uym + 2.*((tauy*tauz + taux)*uzm + (tauy*taux - tauz)*uxm)) * tau;
        uzp = ((1. - taux2 - tauy2 + tauz2)*uzm + 2.*((tauz*taux + tauy)*uxm + (tauz*tauy - taux)*uym)) * tau;

        // Full momentum push
        ux = uxp + cmratio*E.x;
        uy = uyp + cmratio*E.y;
        uz = uzp + cmratio*E.z;

        // Full push
	root = dtco2 / sqrt(ux*ux + uy*uy + uz*uz + 1.);

        curr->x += ux*root;
        curr->y += uy*root;

        // Store
        part->px = uxp*part_mc;
        part->py = uyp*part_mc;
        part->pz = uzp*part_mc;

		
        //This is where the current and charge density would be calculatted.

	//move on to the next one
	curr = curr->next
    }
}
