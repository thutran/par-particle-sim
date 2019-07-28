#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
#include <limits>

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
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


    double size, cell_size, cutoff = get_cutoff();
    int grid_width;
    std::vector<std::vector<particle_t*>> cells;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    // list of old/new cell locations corresponding to the list of particles (to be used in MOVE)
    int *cells_old = (int*) malloc( n * sizeof(int) );
    int *cells_new = (int*) malloc( n * sizeof(int) );

    size = get_size();
    cell_size = (2*sqrt(CELL_MAX_PARTICLES) )*cutoff + std::numeric_limits<double>::epsilon();
    grid_width = (int)ceil(size/cell_size);

    // #pragma omp single
    for (int i = 0; i < grid_width*grid_width; ++i){
        cells.emplace_back(std::vector<particle_t*>());
    }
    #pragma omp parallel for
    for (int i=0; i<n; ++i){
        int cell_x = (int)floor(particles[i].x/cell_size);
        int cell_y = (int)floor(particles[i].y/cell_size);

        particles[i].cell_x = cell_x;
        particles[i].cell_y = cell_y;

        cells_old[i] = cell_y*grid_width + cell_x;
        cells_new[i] = cell_y*grid_width + cell_x;

        #pragma omp critical
        cells[cell_y*grid_width + cell_x].push_back(&particles[i]);
    }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    #pragma omp parallel private(dmin) 
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        //
        //  compute all forces
        //
        #pragma omp for reduction (+:navg) reduction(+:davg)
        for( int i = 0; i < n; i++ ){
            particles[i].ax = particles[i].ay = 0;

            //  _20_|_21_|_22_
            //  _10_|__i_|_12_
            //   00 | 01 | 02

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
            if (cells[cell_i].size() > 1)
                visit_cell_and_apply_force(particles[i], cells, cell_i, &dmin, &davg, &navg);

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
        #pragma omp for
        for (int i=0; i<n; ++i){
            cells_old[i] = particles[i].cell_y*grid_width + particles[i].cell_x;
            
            move( particles[i] );
            
            particles[i].cell_x = (int)floor(particles[i].x/cell_size);
            particles[i].cell_y = (int)floor(particles[i].y/cell_size);
            cells_new[i] = particles[i].cell_y*grid_width + particles[i].cell_x;
        }

        //
        // check if changing cells
        //
        #pragma omp single
        for (int i=0; i<n; ++i){
            int cell_old=cells_old[i];
            int cell_new=cells_new[i];

            if (cell_new != cell_old){
                for (int j=0; j<cells[cell_old].size(); ++j){

                    if (cells[cell_old].at(j) == &particles[i]){
                        cells[cell_old].at(j) = cells[cell_old].back();
                        cells[cell_old].back() = nullptr;
                        cells[cell_old].pop_back();
                        j = n+1; // break the for loop
                    }
                }
                cells[cell_new].push_back(&particles[i]);
            }
        }        

        // #pragma omp barrier

        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
          if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    }

    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
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
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

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
