#include <stdlib.h>
#include <math.h>
#include "list.h"
#include "decs.h"
#include "grid.h"

void output_data2D(int itNum, grid_point **grid_points, int nx, int ny, List particles) {
	// create the file name to output for the grid point data:
	char fname_grid[] = ((char)itNum);
	strcat(fname_grid, "_grid.txt");

	// create a file to output the grid data to:
	FILE *grid_file = fopen(fname_grid, "w"); // write only
    // test to ensure that the file was actually created and exists: 
    if (grid_file == NULL){
    	printf("Error! Could not create file\n"); 
        exit(-1); // must include stdlib.h 
    }//end if

    int i;
    int j;
    fprintf(grid_file, "GRIDPOINTS\n");
    fprintf(grid_file, "itNum=%d\n", itNum);
    fprintf(grid_file, "Time=%lf\n", time);
    fprintf(grid_file, "TimeStep=%lf\n", dt);
    fprintf(grid_file, "GridSize:nx=%d,ny=%d\n", nx, ny);
    for(i = 0; i < ny; i++){
    	for(j = 0; j < nx; j++){
    		//calculate the L2 norm of the E field vector:
            double E = pow(pow((grid_points[i][j].E).x, 2.0) + pow((grid_points[i][j].E).y, 2.0) + pow((grid_points[i][j].E).z, 2.0), 0.5);
            //calculate the L2 norm of the B field vector:
            double B = pow(pow((grid_points[i][j].B).x, 2.0) + pow((grid_points[i][j].B).y, 2.0) + pow((grid_points[i][j].B).z, 2.0), 0.5);
            
            /*write one line to file per grid point in this format:
			xcoord,ycoord,E,B
			*/
			fprintf(grid_file, "%d,%d,%lf,%lf", j, i, E, B);
    	}//end inner for
    }//end outer for
    fclose(grid_file);

    // PARTICLE DATA:
    // create the file name to output for the particle data:
	char fname_particles[] = ((char)itNum);
	strcat(fname_particles, "_particles.txt");

	// create a file to output the particle data to:
	FILE *particles_file = fopen(fname_particles, "w"); // write only
    // test to ensure that the file was actually created and exists: 
    if (particles_file == NULL){
    	printf("Error! Could not create file\n"); 
        exit(-1); // must include stdlib.h 
    }//end if

    fprintf(grid_file, "PARTICLES\n");
    fprintf(grid_file, "itNum=%d\n", itNum);
    fprintf(grid_file, "Time=%lf\n", time);
    fprintf(grid_file, "TimeStep=%lf\n", dt);
    fprintf(grid_file, "GridSize:nx=%d,ny=%d\n", nx, ny);
    
    int particle_ct = -1;
    list_reset_iter(particles);
    while(list_has_next(particles)){
        //calculate the L2 norm of the p vector:
        double p = pow(pow((ptc->p).x, 2.0) + pow((ptc->p).y, 2.0) + pow((ptc->p).z, 2.0), 0.5);
    	/*write one line to file per particle in this format:
			ptcl:#,pos_x,pos_y,p
		*/
        particle_ct++;
		particle *ptc = (particle*) list_get_next(particles);
		fprintf(grid_file, "ptcl:%d,%lf,%lf,%lf\n", particle_ct, (ptc->pos).x, (ptc->pos).y, p);
	}//end while
    fclose(particles_file);
}//end output_grid function

void output_data3D(int itNum, grid_point **grid_points, int nx, int ny, int nz List particles) {
    // create the file name to output for the grid point data:
    char fname_grid[] = ((char)itNum);
    strcat(fname_grid, "_grid.txt");

    // create a file to output the grid data to:
    FILE *grid_file = fopen(fname_grid, "w"); // write only
    // test to ensure that the file was actually created and exists: 
    if (grid_file == NULL){
        printf("Error! Could not create file\n"); 
        exit(-1); // must include stdlib.h 
    }//end if

    int y;
    int x;
    int z;
    fprintf(grid_file, "GRIDPOINTS\n");
    fprintf(grid_file, "itNum=%d\n", itNum);
    fprintf(grid_file, "Time=%lf\n", time);
    fprintf(grid_file, "TimeStep=%lf\n", dt);
    fprintf(grid_file, "GridSize:nx=%d,ny=%d,nz=%d\n", nx, ny, nz);
    for(y = 0; y < ny; y++){
        for(x = 0; x < nx; x++){
            for(z = 0; z < nz, z++){
                //calculate the L2 norm of the E field vector:
                double E = pow(pow((grid_points[i][j].E).x, 2.0) + pow((grid_points[i][j].E).y, 2.0) + pow((grid_points[i][j].E).z, 2.0), 0.5);
                //calculate the L2 norm of the B field vector:
                double B = pow(pow((grid_points[i][j].B).x, 2.0) + pow((grid_points[i][j].B).y, 2.0) + pow((grid_points[i][j].B).z, 2.0), 0.5);
                
                /*write one line to file per grid point in this format:
                xcoord,ycoord,zcoord,E,B
                */
                fprintf(grid_file, "%d,%d,%d,%lf,%lf", x, y, z, E, B);
            }//end innermost for
        }//end middle for
    }//end outer for
    fclose(grid_file);

    // PARTICLE DATA:
    // create the file name to output for the particle data:
    char fname_particles[] = ((char)itNum);
    strcat(fname_particles, "_particles.txt");

    // create a file to output the particle data to:
    FILE *particles_file = fopen(fname_particles, "w"); // write only
    // test to ensure that the file was actually created and exists: 
    if (particles_file == NULL){
        printf("Error! Could not create file\n"); 
        exit(-1); // must include stdlib.h 
    }//end if

    fprintf(grid_file, "PARTICLES\n");
    fprintf(grid_file, "itNum=%d\n", itNum);
    fprintf(grid_file, "Time=%lf\n", time);
    fprintf(grid_file, "TimeStep=%lf\n", dt);
    fprintf(grid_file, "GridSize:nx=%d,ny=%d,nz=%d\n", nx, ny, nz);
    
    int particle_ct = -1;
    list_reset_iter(particles);
    while(list_has_next(particles)){
        //calculate the L2 norm of the p vector:
        double p = pow(pow((ptc->p).x, 2.0) + pow((ptc->p).y, 2.0) + pow((ptc->p).z, 2.0), 0.5);
        /*write one line to file per particle in this format:
            ptcl:#,pos_x,pos_y,pos_z,p
        */
        particle_ct++;
        particle *ptc = (particle*) list_get_next(particles);
        fprintf(grid_file, "ptcl:%d,%lf,%lf,%lf,%lf\n", particle_ct, (ptc->pos).x, (ptc->pos).y, (ptc->pos).z, p);
    }//end while
    fclose(particles_file);
}//end output_data3D function

    //output is too verbose in this version of the function...
// void output_grid(int itNum, grid_point **grid_points, int nx, int ny, List particles) {
//     // create the file name to output for the grid point data:
//     char fname_grid[] = ((char)itNum);
//     strcat(fname_grid, "_grid.txt");

//     // create a file to output the grid data to:
//     FILE *grid_file = fopen(fname_grid, "w"); // write only
//     // test to ensure that the file was actually created and exists: 
//     if (grid_file == NULL){
//         printf("Error! Could not create file\n"); 
//         exit(-1); // must include stdlib.h 
//     }//end if

//     int i;
//     int j;
//     fprintf(grid_file, "GRIDPOINTS\n");
//     fprintf(grid_file, "itNum=%d\n", itNum);
//     fprintf(grid_file, "Time=%lf\n", time);
//     fprintf(grid_file, "TimeStep=%lf\n", dt);
//     fprintf(grid_file, "GridSize:nx=%d,ny=%d\n", nx, ny);
//     for(i = 0; i < ny; i++){
//         for(j = 0; j < nx; j++){
//             /*write one line to file per grid point in this format:
//             coordinates:i,j,Ex=lf,Ey=lf,Ez=lf,Bx=lf,By=lf,Bz=lf
//             */
//             fprintf(grid_file, "coordinates=%d,%d,Ex=%lf,Ey=%lf,Ez=%lf,Bx=%lf,By=%lf,Bz=%lf\n", i, j, (grid_points[i][j].E).x, (grid_points[i][j].E).y, (grid_points[i][j].E).z, (grid_points[i][j].B).x, (grid_points[i][j].B).y, (grid_points[i][j].B).z);
//         }//end inner for
//     }//end outer for
//     fclose(grid_file);

//     // PARTICLE DATA:
//     // create the file name to output for the particle data:
//     char fname_particles[] = ((char)itNum);
//     strcat(fname_particles, "_particles.txt");

//     // create a file to output the particle data to:
//     FILE *particles_file = fopen(fname_particles, "w"); // write only
//     // test to ensure that the file was actually created and exists: 
//     if (particles_file == NULL){
//         printf("Error! Could not create file\n"); 
//         exit(-1); // must include stdlib.h 
//     }//end if

//     fprintf(grid_file, "PARTICLES\n");
//     fprintf(grid_file, "itNum=%d\n", itNum);
//     fprintf(grid_file, "Time=%lf\n", time);
//     fprintf(grid_file, "TimeStep=%lf\n", dt);
//     fprintf(grid_file, "GridSize:nx=%d,ny=%d\n", nx, ny);
    
//     int particle_ct = -1;
//     list_reset_iter(particles);
//     while(list_has_next(particles)){
//         /*write one line to file per particle in this format:
//             pos_x=#,pos_y=#,p_x=#,p_y=#,p_z=#
//         */
//         particle_ct++;
//         particle *ptc = (particle*) list_get_next(particles);
//         fprintf(grid_file, "particle:%d,pos_x=%lf,pos_y=%lf,p_x=%lf,p_y=%lf,p_z=%lf\n", particle_ct, (ptc->pos).x, (ptc->pos).y, (ptc->pos).x, (ptc->p).x, (ptc->p).y, (ptc->p).z);
//     }//end while
//     fclose(particles_file);
// }//end output_grid function
          
