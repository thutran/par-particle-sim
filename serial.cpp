#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <iostream>

// default max number of particles that can co-exist in 1 cell at a time
#if !defined (CELL_MAX_PARTICLES)
#define CELL_MAX_PARTICLES 4
#endif

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;
    double size, cell_size, cutoff = get_cutoff();
    int grid_width;
    std::vector<std::vector<particle_t*>> cells;

    // BEGIN arguments--------------------
    if( find_option( argc, argv, "-h" ) >= 0 ){
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;
    // END arguments--------------------


    // BEGIN initialize--------------------
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    size = get_size();
    // cell_size = (2*sqrt(CELL_MAX_PARTICLES) + 0.01)*cutoff ;
    cell_size = (2*sqrt(CELL_MAX_PARTICLES) )*cutoff ;
    grid_width = (int)ceil(size/cell_size);

    // cells.reserve(grid_width*grid_width);
    init_grid(n, particles, cell_size, grid_width, cells);

    // for (int i=n-1; i<n; ++i){
    //     printf("x %f\t y %f\t cell x %d\t cell y %d\n", particles[0].x, particles[0].y, particles[0].cell_x, particles[0].cell_y);
    //     printf("x %f\t y %f\t cell x %d\t cell y %d\n", particles[i].x, particles[i].y, particles[i].cell_x, particles[i].cell_y);
    // }
    // printf("size vector particle_t: %d\n", sizeof(std::vector<particle_t*>));

    // std::vector<std::vector<particle_t*> > cells;
    // END initialize--------------------
    

    // BEGIN simulate--------------------
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ ){
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        //
        //  compute forces
        //
        for( int i = 0; i < n; i++ ){
            particles[i].ax = particles[i].ay = 0;

            //  _00_|_01_|_02_
            //  _10_|__i_|_12_
            //   20 | 21 | 22

            const double par_x = particles[i].x;
            const double par_y = particles[i].y;
            const int cell_x = particles[i].cell_x;
            const int cell_y = particles[i].cell_y;
            const int cell_own = cell_y*grid_width + cell_x;
            const double x_mid_cell = (cell_x + 0.5)*cell_size;
            const double y_mid_cell = (cell_y + 0.5)*cell_size;
            int cell_i;

            // own cell
            cell_i = cell_own;
            if (cells[cell_i].size() > 1){
                // for (auto it=cells[cell_i].begin(); it!=cells[cell_i].end(); ++it){
                //     // std::cout << "iterator: "<< *it << "\t particle: " << &particles[i] << "\n";
                //     // std::cout << "iterator-p: "<< (*it)->x << "\t particle-p: " << (&particles[i])->x << "\n";
                //     // std::cout << "iterator-p: "<< (*(*it)).x << "\t particle-p: " << particles[i].x << "\n";
                //     // test if iterator point to the particle itself
                //     if (*it != &particles[i]){
                //         apply_force(particles[i], (*(*it)), &dmin, &davg, &navg);
                //     }
                // }

                visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
            }

            // left
            if (fequal(par_y, y_mid_cell) && fequal(par_x, cell_x*cell_size)) {
                // printf("%s\n", "particle at y_mid_cell left edge");
                // if having left neighbor
                if (cell_x > 0){
                    cell_i = cell_own - 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
            }
            
            // // right
            if (fequal(par_y, y_mid_cell) && fequal(par_x, cell_x*cell_size + cell_size)) {
                // printf("%s\n", "particle at y_mid_cell right edge");
                // if having right neighbor
                if (cell_x < (grid_width-1)) {
                    cell_i = cell_own + 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
            }

            // bottom
            if (fequal(par_x, x_mid_cell) && fequal(par_y, cell_y*cell_size)) {
                // printf("%s\n", "particle at x_mid_cell bottom edge");
                // if having bottom neighbor
                if (cell_y > 0) {
                    cell_i = cell_own - grid_width;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
            }

            // top
            if (fequal(par_x, x_mid_cell) && fequal(par_y, cell_y*cell_size + cell_size)) {
                // printf("%s\n", "particle at x_mid_cell top edge");
                // if having top neighbor
                if (cell_y < (grid_width-1)) {
                    cell_i = cell_own + grid_width;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
            }

            // bottom-left
            if (fless(par_x, x_mid_cell) && fless(par_y, y_mid_cell)) {
                // left neighbor
                if (cell_x > 0){
                    cell_i = cell_own - 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
                // bottom neighbor
                if (cell_y > 0) {
                    cell_i = cell_own - grid_width;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
                // bottom-left neighbor
                if (cell_x > 0 && cell_y > 0){
                    cell_i = cell_own - grid_width - 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
            }

            // bottom-right
            if (fgreater(par_x, x_mid_cell) && fless(par_y, y_mid_cell)) {
                // right neighbor
                if (cell_x < (grid_width-1)) {
                    cell_i = cell_own + 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
                // bottom neighbor
                if (cell_y > 0) {
                    cell_i = cell_own - grid_width;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
                // bottom-right neighbor
                if (cell_x < (grid_width-1) && cell_y > 0){
                    cell_i = cell_own - grid_width + 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
            }

            // top-right
            if (fgreater(par_x, x_mid_cell) && fgreater(par_y, y_mid_cell)) {
                // right neighbor
                if (cell_x < (grid_width-1)) {
                    cell_i = cell_own + 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
                // top neighbor
                if (cell_y < (grid_width-1)) {
                    cell_i = cell_own + grid_width;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
                // top-right neighbor
                if (cell_x < (grid_width-1) && cell_y < (grid_width-1)) {
                    cell_i = cell_own + grid_width + 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
            }

            // top-left
            if (fless(par_x, x_mid_cell) && fgreater(par_y, y_mid_cell)) {
                // left neighbor
                if (cell_x > 0){
                    cell_i = cell_own - 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
                // top neighbor
                if (cell_y < (grid_width-1)) {
                    cell_i = cell_own + grid_width;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
                // top-left neighbor
                if (cell_x > 0 && cell_y < (grid_width-1)) {
                    cell_i = cell_own + grid_width - 1;
                    if (!cells[cell_i].empty())
                        visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);
                }
            }

        
        }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move_to_cell(particles[i], cell_size, grid_width, cells);

        // for( int i = 0; i < n; i++ ) 
        //     move( particles[i] );		



        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    // END simulate--------------------

    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);


    if( find_option( argc, argv, "-no" ) == -1 ){
        if (nabsavg) absavg /= nabsavg;
        // 
        //  -the minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);

    // printf("x_%f\t y_%f\t cell_x_%d\t cell_y_%d\n", particles[n-1].x, particles[n-1].y, particles[n-1].cell_x, particles[n-1].cell_y);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
