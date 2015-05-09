#include <stdlib.h>
#include <math.h>

#include "decs.h"
#include "list.h"
#include "dynamics.h"
#include "mpicomm.h"

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
    // double norm = 1./(dx*dy*dz);
    double norm = 1.;
    double xrf = 1.-xlf;
    double ydf = 1. - yuf;
    double zff = 1. - znf;
    interped.x = norm * (field000.x*xlf*yuf*znf + field001.x*xlf*yuf*zff + field010.x*xlf*ydf*znf + field100.x*xrf*yuf*znf + field011.x*xlf*ydf*zff + field101.x*xrf*yuf*zff + field110.x*xrf*ydf*znf + field111.x*xrf*ydf*zff);
    interped.y = norm * (field000.y*xlf*yuf*znf + field001.y*xlf*yuf*zff + field010.y*xlf*ydf*znf + field100.y*xrf*yuf*znf + field011.y*xlf*ydf*zff + field101.y*xrf*yuf*zff + field110.y*xrf*ydf*znf + field111.y*xrf*ydf*zff);
    interped.z = norm * (field000.z*xlf*yuf*znf + field001.z*xlf*yuf*zff + field010.z*xlf*ydf*znf + field100.z*xrf*yuf*znf + field011.z*xlf*ydf*zff + field101.z*xrf*yuf*zff + field110.z*xrf*ydf*znf + field111.z*xrf*ydf*zff);
    return interped;
}

// Particle pusher!!!!
static void push_one_cell(tree cell) {
	List part_list = cell.particles;

	list_reset_iter(&part_list);
	if (!list_has_next(part_list))
			return;

    // Declare variables!
    double ux, uy, uz;
    double root;
    int xl, yu, zn;
    double xrf, ydf, zff;
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
    double dtco2 = C * dto2;
    double dtfac = .5*dt ;//times some weighting factor, that I think is 1 anyways
    double third = 1./3.;
    double idty = idt*idy;
    double idtx = idt*idy;
    double idxy = idx*idy;
    double tau, taux, taux2, tauy, tauy2, tauz, tauz2;
    
	particle *curr;
    //loop over all the particles
    while (list_has_next(part_list)) {
		curr = (particle*) list_get_next(&part_list);

		double cmratio = curr->charge/curr->mass;
        part_mc = C*curr->mass;
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
		// This check stays the same when parallel, since each processor will have the appropriate ghost cells
		if ((((curr->pos).x <= 0 || (curr->pos).y <= 0) || (curr->pos).z <= 0) || (((curr->pos).x >= x_max || (curr->pos).y >= y_max) || (curr->pos).z >= z_max)){
			list_pop(&part_list);
			continue;
		}

		//Do interpolation to find e and b here.
        // x-left, y-up, and z-near indices
		// Subtract out the local min to get the correct indicies
        xl = floor(((curr->pos).x - pxmin) * idx);
        yu = floor(((curr->pos).y - pymin) * idy);
		zn = floor(((curr->pos).z - pzmin) * idz);
		// x-right fraction, ...
		// This stays the same for parallel, I think
        xrf = ((curr->pos).x - xl*dx) / dx;
        ydf = ((curr->pos).y - yu*dy) / dy;
		zff = ((curr->pos).z - zn*dz) / dz;
        /*E = interp3(grid[xl][yu][zn].E, grid[xl][yu][zn+1].E, grid[xl][yu+1][zn].E, grid[xl+1][yu][zn].E, grid[xl][yu+1][zn+1].E, grid[xl+1][yu][zn+1].E, grid[xl][yu+1][zn+1].E, grid[xl+1][yu+1][zn+1].E, xrf, ydf, zff);
        B = interp3(grid[xl][yu][zn].B, grid[xl][yu][zn+1].B, grid[xl][yu+1][zn].B, grid[xl+1][yu][zn].B, grid[xl][yu+1][zn+1].B, grid[xl+1][yu][zn+1].B, grid[xl][yu+1][zn+1].B, grid[xl+1][yu+1][zn+1].B, xrf, ydf, zff);*/

		TreeNode *cellIter = cell.root;
		//Find the finest cell that contains the particle
		while (cellIter->children != NULL){
			if (xrf < .5){
				if (ydf < .5){
					if (zff < .5){
						cellIter = cellIter->children[0];
						xrf*=2;
						ydf*=2;
						zff*=2;
					}
					else{
						cellIter = cellIter->children[4];
						xrf*=2;
						ydf*=2;
						zff=(zff-.5)*2;
					}
				}
				else{
					if (zff < .5){
						cellIter = cellIter->children[2];
						xrf*=2;
						ydf=(ydf-.5)*2;
						zff*=2;
					}
					else{
						cellIter = cellIter->children[6];
						xrf*=2;
						ydf=(ydf-.5)*2;
						zff=(zff-.5)*2;
					}
				}
			}
			else{
				if (ydf < .5){
					if (zff < .5){
						cellIter = cellIter->children[1];
						xrf=(xrf-.5)*2;
						ydf*=2;
						zff*=2;
					}
					else{
						cellIter = cellIter->children[5];
						xrf=(xrf-.5)*2;
						ydf*=2;
						zff=(zff-.5)*2;
					}
				}
				else{
					if (zff < .5){
						cellIter = cellIter->children[3];
						xrf=(xrf-.5)*2;
						ydf=(ydf-.5)*2;
						zff*=2;
					}
					else{
						cellIter = cellIter->children[7];
						xrf=(xrf-.5)*2;
						ydf=(ydf-.5)*2;
						zff=(zff-.5)*2;
					}
				}
			}
		}

        //Do interpolation with the new grid_cell
        E = interp3(cellIter->points[0]->E, cellIter->points[1]->E, cellIter->points[2]->E, cellIter->points[4]->E, cellIter->points[3]->E, cellIter->points[5]->E, cellIter->points[6]->E, cellIter->points[7]->E, 1.-xrf, 1.-ydf, 1.-zff);
        B = interp3(cellIter->points[0]->B, cellIter->points[1]->B, cellIter->points[2]->B, cellIter->points[4]->B, cellIter->points[3]->B, cellIter->points[5]->B, cellIter->points[6]->B, cellIter->points[7]->B, 1.-xrf, 1.-ydf, 1.-zff);
        
        // Update momenta to u_-, from Birdsall and Langdon
        uxm = ux + cmratio * E.x;
        uym = uy + cmratio * E.y;
        uzm = uz + cmratio * E.z;

        // Do the Borris rotation and update momenta to u_+, from Birdsall and Langdon
        root = cmratio / sqrt(uxm*uxm + uym*uym + uzm*uzm + 1.0);

        taux = B.x*root;
        taux2 = taux*taux;
        tauy = B.y*root;
        tauy2 = tauy*tauy;
        tauz = B.z*root;
        tauz2 = tauz*tauz;

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
			list_pop(&part_list);
			continue;
		}

        // Store
        (curr->p).x = ux*part_mc;
        (curr->p).y = uy*part_mc;
        (curr->p).z = uz*part_mc;

		
        //This is where the current and charge density would be calculatted.


		// Pass the particles to neighbor cells if necessary
        // Ending x-left, y-up, and z-near indices
		// Subtract out the local min to get the correct indicies
		int xle, yue, zne;
        xle = floor(((curr->pos).x - pxmin) * idx);
        yue = floor(((curr->pos).y - pymin) * idy);
		zne = floor(((curr->pos).z - pzmin) * idz);

		// Check if cell has changed
		// Guarenteed to still be in a cell or ghost cell controled by proc
		if (xle != xl || yue != yu || zne != zn){
			// add curr to the next_list of grid[xle][yue][zne]
			particle_pass(&(grid[xle][yue][zne]->next_list), &(grid[xl][yu][zn]->part_list), curr);
		}

    } 
}

// Calls the pusher and cleans up afterward
void push_particles(tree ****grid) {
	tree *curCell = NULL;
	int i,j,k;
	for (i = imin; i < imax; ++i) {
		for (j = jmin; j < jmax; ++j) {
			for (k = kmin; k < kmax; ++k) {
				// Check if valid cell
				curCell = grid[i][j][k];
				if (curCell != NULL) {
					if (curCell->owner == pid) {
						push_one_cell(*curCell);
					}
				}
			}
		}
	}
	
	//Allocate buffers and do the sends and recvs
	neighbor *neighbors; // array of lists of particle lists to send to each neighbor proc
	neighbors = (neighbor*) calloc( nProcs*sizeof(neighbor) );
	// loop over each ghost cell.
	// loop over entire grid to find the ghost cells, this is the easiest way to find them
	// this can't be integrated with above identical loop since the push must be completed before checking to see which pushed things need to be sent to neighbors
	for (i = imin; i < imax; ++i) {
		for (j = jmin; j < jmax; ++j) {
			for (k = kmin; k < kmax; ++k) {
				curCell = grid[i][j][k];
				if (curCell != NULL) {
					int owner = curCell->owner;
					if (owner != pid) {
						if (neighbors[owner] == NULL) {
							neighbors[owner] = neighbor_init(i);
						}
						neighbor_add_cell(&(neighbors[owner]), curCell);
					}
				}
			}
		}
	}

	// array to the requests so that we can wait for all receives to finish
	MPI_Request *cell_count_requests = (MPI_Request*) malloc( nProcs*sizeof(MPI_Request) );

	// send # of cell we will send
	for (i = 0; i < nProcs; ++i) {
		if (neighbors[i] != NULL) {
			cell_count_requests[i] = neighbor_send_cell_count(neighbors[i]);
		} else {
			cell_count_requests[i] = MPI_REQUEST_NULL;
		}
	}

	//TODO: can we ignore status for waitall?
	MPI_Waitall(cell_count_requests, MPI_STATUSES_IGNORE);

	MPI_Request **cell_requests = (MPI_Request**) malloc( nProcs*sizeof(MPI_Request) );
	// send the cells themselves
	for (i = 0; i < nProcs; ++i) {
		if (neighbors[i] != NULL) {
			cell_requests[i] = neighbor_send_cells(neighbors[i]);
		}
	}

	// wait to finish recving data from all your neighbors
	for (i = 0; i < nProcs; ++i) {
		if (neighbors[i] != NULL) {
			//TODO: can we ignore status for waitall?
			MPI_Waitall(cell_requests[i], MPI_STATUSES_IGNORE);
		}
	}

	// for buffs that hane recieved
	//		do the unpacking
	//
	// Also do some frees
	// loop over your i-th neighbor
	for (i = 0; i < nProcs; ++i) {
		neighbor n = neighbors[i];
		if (n == NULL) { // only receive from actual neighbors
			continue;
		}
		// loop over the j-th cell you received from your i-th neighbor
		int iCell;
		for (iCell = 0; iCell < n.ncellrecvs; ++iCell) {
			// figure out which cell to put the particles in based on the coords of the 1st particle sent
			tree destination = determine_cell_to_put_in(n.recvbufs[iCell][0].pos);
			// loop over the k-th particle received from the j-th cell
			int iPart;
			for (iPart = 0; iPart < n.recvlen[iCell]; ++iPart) {
				// add the receivde particle to the appropriate cell's particle list
				list_add(destination, n.recbufs[iCell][iPart])
			}
			free(n.recvbufs[iCell]);
		}
		// now free all of the remaining buffers, recvbufs already freed when saving the particles
		//TODO: we can't save them because we don't know how long they will need to be next time step. maybe we can use realloc?
		for (iCell = 0; i < n.ncellsends; ++iCell) {
			free(n.sendbufs[i]);
		}
		free(n.sendbufs);
		free(n.recvbufs);
	}

	
	//Combine all the lists
	//This is the last step
	for (i=imin; i<imax; i++){
		for (j=jmin; j<jmax; j++){
			for (k=kmin; k<kmax; k++){
				curCell = grid[i][j][k];
				if (curCell != NULL){
					// Add the next_list to the current list
					list_combine(&(grid[i][j][k]->part_list), &(grid[i][j][k]->next_list));
				}
			}
		}
	}

}


