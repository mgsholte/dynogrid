#include "decs.h"
#include "math.h"
#include "list.h"

// Some globals that I am assuming exist
// dx, dy -- cell size
// dt -- time step

// Interpolation helper
vec3 interp(vec3 field00, vec3 field01, vec3 field10, vec3 field11, double xlf, double yuf) {
    vec3 interped;
    double norm = 1./(dx*dy);//Make this a global for speed
    double xrf = 1. - xlf;
    double ydf = 1.- yuf;
    // Not sure if these are proper vec3 accesses.  May need to unroll.
        // Bilinear interpolation from Wikipedia.  There are certainly
        // better interpolation methods to use for PIC codes.
    interped.x = norm * ( field00.x*xrf*yuf + field10.x*xlf*yuf + field01.x*xrf*ydf + field11.x*xlf*ydf);
    interped.y = norm * ( field00.y*xrf*yuf + field10.y*xlf*yuf + field01.y*xrf*ydf + field11.y*xlf*ydf);
    interped.z = norm * ( field00.z*xrf*yuf + field10.z*xlf*yuf + field01.z*xrf*ydf + field11.z*xlf*ydf);
	return interped;
}

// 3D Interpalation helper
vec3 interp3(vec3 field000, vec3 field001, vec3 field010, vec3 field100, vec3 field011, vec3 field101, vec3 field110, vec3 field111, double xlf, double yuf, double znf){
	vec3 interped;
	double norm = 1./(dx*dy*dz);
	double xrf = 1.-xlf;
	double ydf = 1. - yuf;
	double zff = 1. - znf;
	interped.x = norm * (field000.x*xlf*yuf*znf + field001.x*xlf*yuf*zff + field010.x*xlf*ydf*znf + field100.x*xrf*yuf*znf + field011.x*xlf*ydf*zff + field101.x*xrf*yuf*zff + field110.x*xrf*ydf*znf + field111.x*xrf*ydf*zff);
	interped.y = norm * (field000.y*xlf*yuf*znf + field001.y*xlf*yuf*zff + field010.y*xlf*ydf*znf + field100.y*xrf*yuf*znf + field011.y*xlf*ydf*zff + field101.y*xrf*yuf*zff + field110.y*xrf*ydf*znf + field111.y*xrf*ydf*zff);
	interped.z = norm * (field000.z*xlf*yuf*znf + field001.z*xlf*yuf*zff + field010.z*xlf*ydf*znf + field100.z*xrf*yuf*znf + field011.z*xlf*ydf*zff + field101.z*xrf*yuf*zff + field110.z*xrf*ydf*znf + field111.z*xrf*ydf*zff);
	return interped;
}

// Particle pusher!!!!
void push_particles(grid_point ***grid, List part_list) {

	list_reset_iter(&part_list);

    // Declare variables!
	if (!list_has_next(part_list))
			return;
    particle *curr = list_get_next(&part_list);
	Node *prev = part_list.sentinel;
    double ux, uy, uz;
    double root;
    int xl, yu, zn;
	double xlf, yuf, znf;
    vec3 E, B;
    double uxm, uym, uzm;
    double uxp, uyp, uzp;
    double part_mc, ipart_mc;
    //  Just copying multiplication factors from EPOCH
    double idx = 1./dx;
    double idy = 1./dy;
	double idz = 1./dz;
    double idt = 1./dt;
    double dto2 = dt/2.;
    double dtco2 = c * dto2;
    double dtfac = .5*dt ;//times some weighting factor, that I think is 1 anyways
    double third = 1./3.;
    double idty = idt*idy;
    double idtx = idt*idy;
    double idxy = idx*idy;
	double tau, taux, taux2, tauy, tauy2, tauz, tauz2;
    
    //loop over all the particles
    while (true){
		double cmratio = curr->charge/curr->mass;
        part_mc = c*curr->mass;
        ipart_mc = 1./part_mc;

        // u is gamma*v, see Birdsall and Langdon sectoin 15-4
		ux = (curr->p).x * ipart_mc;
		uy = (curr->p).y * ipart_mc;
		uz = (curr->p).z * ipart_mc;
		//Calculate velocity
		root = dtco2 / sqrt(ux*ux + uy*uy + uz*uz + 1.);

		//Move half timestep
		(curr->pos).x += ux * root;
		(curr->pos).y += uy * root;
		(curr->pos).z += uz * root;

		// Check if out of bounds
		if ((((curr->pos).x <= 0 || (curr->pos).y <= 0) || (curr->pos).z <= 0) || (((curr->pos).x >= x_max || (curr->pos).y >= y_max) || (curr->pos).z >= z_max)){
			list_pop(&part_list, prev);
			//move on to the next one
			if (list_has_next(part_list)){
				curr = list_get_next(&part_list);
				continue;
			}
			else{
				list_reset_iter(&part_list);
				return;
			}
		}

		//Do interpolation to find e and b here.
        // x-left and y-up indices
        xl = floor((curr->pos).x * idx);
        yu = floor((curr->pos).y * idy);
		zn = floor((curr->pos).z * idz);
        xlf = ((curr->pos).x - xl*dx) / dx;
        yuf = ((curr->pos).y - yu*dy) / dy;
		znf = ((curr->pos).z - zn*dz) / dz;
        E = interp3(grid[xl][yu][zn].E, grid[xl][yu][zn+1].E, grid[xl][yu+1][zn].E, grid[xl+1][yu][zn].E, grid[xl][yu+1][zn+1].E, grid[xl+1][yu][zn+1].E, grid[xl][yu+1][zn+1].E, grid[xl+1][yu+1][zn+1].E, xlf, yuf, znf);
        B = interp3(grid[xl][yu][zn].B, grid[xl][yu][zn+1].B, grid[xl][yu+1][zn].B, grid[xl+1][yu][zn].B, grid[xl][yu+1][zn+1].B, grid[xl+1][yu][zn+1].B, grid[xl][yu+1][zn+1].B, grid[xl+1][yu+1][zn+1].B, xlf, yuf, znf);
        
        // Update momenta to u_-, from Birdsall and Langdon
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
        tauz2 *= tauz;

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

        (curr->pos).x += ux*root;
        (curr->pos).y += uy*root;
		(curr->pos).z += uz*root;

		// Check if out of bounds
		if ((((curr->pos).x <= 0 || (curr->pos).y <= 0) || (curr->pos).z <= 0) || (((curr->pos).x >= x_max || (curr->pos).y >= y_max) || (curr->pos).z >= z_max)){
			list_pop(&part_list, prev);
			//move on to the next one
			if (list_has_next(part_list)){
				curr = list_get_next(&part_list);
				continue;
			}
			else{
				list_reset_iter(&part_list);
				return;
			}
		}

        // Store
        (curr->p).x = ux*part_mc;
        (curr->p).y = uy*part_mc;
        (curr->p).z = uz*part_mc;

		
        //This is where the current and charge density would be calculatted.

		//move on to the next one
		if (list_has_next(part_list)){
			prev = part_list.iter;
			curr = list_get_next(&part_list);
		}
		else{
			list_reset_iter(&part_list);
			return;
		}
    }
}
