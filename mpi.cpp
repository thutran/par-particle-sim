#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <cstring>
// #include <memory>

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
 
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
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

    // paritcles array
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    // particle_t *particles_prev = (particle_t*) malloc( n * sizeof(particle_t) );
    std::vector<std::vector<particle_t>> rows;

    // grid sizes
    double size, cell_size, cutoff = get_cutoff();
    int grid_width;
    
    set_size( n );
    size = get_size();
    cell_size = (2*sqrt(CELL_MAX_PARTICLES) )*cutoff + std::numeric_limits<double>::epsilon();
    grid_width = (int)ceil(size/cell_size);
    

    //////////////////////////////////////////////////////////////
    //  set up MPI
    //
    MPI_Request request;
    MPI_Status status;

    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( sizeof(particle_t), MPI_BYTE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

    
    // local grid structure
    std::vector<std::vector<particle_t* > > local_cells;
    std::vector<std::vector<particle_t > > local_rows;
    std::vector<particle_t > above_row;
    std::vector<particle_t > below_row;
    
    //
    // prepare resource partitioning
    //
    // proc 
    //   2 |    row 4   |
    //   1 |    row 3   |
    //   1 |    row 2   |
    //   0 |    row 1   |   
    //   0 |    row 0   |
    // total number of particles in a row of the grid
    int *particles_in_row = (int*) malloc(grid_width *sizeof(int));
    // capacity: each proc is responsible for >=0 rows of the grid
    int row_per_proc = grid_width < n_proc ? 1 : (int)ceil((double)grid_width/n_proc);
    // row offset to mark starting row for each proc
    int *row_offset_arr = (int*) malloc((n_proc+1) * sizeof(int));
    for (int i = 0; i < n_proc+1; ++i){
            row_offset_arr[i] = min(i*row_per_proc, grid_width);
    }
    // total number of rows to be assigned to each proc
    int *row_total_arr = (int*) malloc(n_proc * sizeof(int));
    for (int i = 0; i < n_proc; ++i){
        row_total_arr[i] = row_offset_arr[i+1] - row_offset_arr[i];
    }
    int local_row_offset = row_offset_arr[rank];
    // printf("rank %d: grid_width %d, row_per_proc %d, row_offset %d, row_total %d\n", rank, grid_width, row_per_proc, local_row_offset, row_total_arr[rank]);
        


    //
    // all ranks shall know the grid structure
    //
    for (int i = 0; i < grid_width*grid_width; ++i){
        local_cells.emplace_back(std::vector<particle_t*>());
    }
    for (int i = 0; i < row_total_arr[rank]; ++i){
        local_rows.emplace_back(std::vector<particle_t>());
        // printf("rank %d local_rows %d size %d | particles in row %d\n", rank, i, local_rows[i].size(), particles_in_row[local_row_offset + i]);
    }
    


    //
    // initialize particles at rank 0 only
    //
    if( rank == 0 ){
        init_particles( n, particles );
        // memcpy (&particles_prev[0], &particles[0], n * sizeof(particle_t));

        for (int i = 0; i < grid_width; ++i){
            rows.emplace_back(std::vector<particle_t>());
        }
        for (int i=0; i<n; ++i){
            int cell_x = (int)floor(particles[i].x/cell_size);
            int cell_y = (int)floor(particles[i].y/cell_size);

            particles[i].cell_x = cell_x;
            particles[i].cell_y = cell_y;
        }
    }

    //
    //  MAIN sim
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        // prepare/reset the main rows vector at rank 0
        if (rank == 0){
            //  save current step if necessary (slightly different semantics than in other codes)
            if( find_option( argc, argv, "-no" ) == -1 )
              if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );

            // fill in the main rows vector from the particles array
            for (int i = 0; i < grid_width; ++i){
                if (!rows[i].empty())
                    std::vector<particle_t>().swap(rows[i]);
            }
            for (int i=0; i<n; ++i){
                rows[particles[i].cell_y].push_back(particles[i]);
            }
            // count number of particles in rows
            for (int i = 0; i < grid_width; ++i){
                particles_in_row[i] = rows[i].size();
            }
        }

        // broadcast updated the total number of particles in each row to all proc
        MPI_Bcast(particles_in_row, grid_width, MPI_INT, 0, MPI_COMM_WORLD);

        // reset local_cells, local_rows, above_row, below_row
        for (int i = 0; i < grid_width*grid_width; ++i){
            if (!local_cells[i].empty())
                std::vector<particle_t*>().swap(local_cells[i]);
        }
        for (int i = 0; i < row_total_arr[rank]; ++i){
            if (!local_rows[i].empty())
                std::vector<particle_t>().swap(local_rows[i]);
        }
        if (!above_row.empty())
            std::vector<particle_t>().swap(above_row);
        if (!below_row.empty())
            std::vector<particle_t>().swap(below_row);

        //
        // prepare row vectors to the correct capacity provided by particles_in_row
        //
        // resize local_rows to the correct capacity size
        for (int i = 0; i < row_total_arr[rank]; ++i){
            local_rows[i].resize(particles_in_row[local_row_offset + i]);
            // printf("rank %d local_rows %d: %d\n", rank, i, local_rows[i].size());
        }
        // resize below_row with the number of particles in the top row of the (rank-1) processor
        if (rank > 0){
            if (row_total_arr[rank] >0)
                below_row.resize(particles_in_row[local_row_offset - 1]);
        }
        // resize above_row with the number of particles in the bottom row of the (rank+1) processor
        if (rank < (n_proc - 1)){
            if (row_total_arr[rank+1] >0)
                above_row.resize(particles_in_row[row_offset_arr[rank+1] ]);
        }

        //
        // allocate rows to proc
        //
        // send rows from rank 0 to local_rows in all proc including itself
        if (rank==0){
            // printf("rank 0 rows.back().size() size before send: %d\n", rows.back().size());
            for (int i = 0; i < n_proc; ++i){
                int row_offset = row_offset_arr[i];
                    for (int j = 0; j < row_total_arr[i]; ++j)
                        if (particles_in_row[row_offset + j] >0)
                            MPI_Isend(&rows[row_offset + j][0], particles_in_row[row_offset + j], PARTICLE, i, 100, MPI_COMM_WORLD, &request);
                            // MPI_Send(&rows[row_offset + j][0], particles_in_row[row_offset + j], PARTICLE, i, 100, MPI_COMM_WORLD );
            }
            // printf("rank 0 rows.back().size() size after send: %d\n", rows.back().size());
        } 
        //
        // fill local_rows at each proc
        for (int j = 0; j < row_total_arr[rank]; ++j)
            if (particles_in_row[local_row_offset + j] >0)
                // MPI_Irecv(&local_rows[j][0], particles_in_row[row_offset + j], PARTICLE, 0, 100, MPI_COMM_WORLD, &request);
                MPI_Recv(&local_rows[j][0], particles_in_row[local_row_offset + j], PARTICLE, 0, 100, MPI_COMM_WORLD, &status);
        //
        // rows from neighbors
        //
        // send particles to the processor right below itself
        if (rank > 0)
            if (row_total_arr[rank] >0)
                MPI_Isend(&local_rows[0][0], particles_in_row[local_row_offset], PARTICLE, rank-1, 101, MPI_COMM_WORLD, &request);
        // 
        // receive particles from the processor right above itself
        if (rank < (n_proc - 1))
            if (row_total_arr[rank+1] >0)
                MPI_Recv(&above_row[0], particles_in_row[row_offset_arr[rank+1]], PARTICLE, rank+1, 101, MPI_COMM_WORLD, &status);

        //
        // send particles to the processor right above itself
        if (rank < (n_proc - 1))
            if (row_total_arr[rank+1] >0)
                MPI_Isend(&local_rows.back().at(0), particles_in_row[row_offset_arr[rank+1] - 1], PARTICLE, rank+1, 102, MPI_COMM_WORLD, &request);
        //
        // receive particles from the processor right below itself
        if (rank > 0)
            if (row_total_arr[rank] >0)
                MPI_Recv(&below_row[0], particles_in_row[local_row_offset - 1], PARTICLE, rank-1, 102, MPI_COMM_WORLD, &status);
        
        // MPI_Barrier(MPI_COMM_WORLD);
        //
        // fill local_cells with local_rows, above_row, below_row
        //
        for (int i = 0; i < local_rows.size(); ++i)
            fill_grid(local_rows[i].size(), &local_rows[i][0], grid_width, local_cells);
        if (rank >0)
            fill_grid(below_row.size(), &below_row[0], grid_width, local_cells);
        if (rank < (n_proc-1))
            fill_grid(above_row.size(), &above_row[0], grid_width, local_cells);


        // MPI_Barrier(MPI_COMM_WORLD);
        //------------------------------------------------//
        //
        // apply force
        //
        for (int i = 0; i < local_rows.size(); ++i){
            for (int j = 0; j < local_rows[i].size(); ++j){
                local_rows[i][j].ax = local_rows[i][j].ay = 0;

                //  _20_|_21_|_22_
                //  _10_|__i_|_12_
                //   00 | 01 | 02

                const double par_x = local_rows[i][j].x;
                const double par_y = local_rows[i][j].y;
                const int cell_x = local_rows[i][j].cell_x;
                const int cell_y = local_rows[i][j].cell_y;
                const int cell_own = cell_y*grid_width + cell_x;
                const double x_mid_cell = (cell_x + 0.5)*cell_size;
                const double y_mid_cell = (cell_y + 0.5)*cell_size;
                int cell_i;

                // own cell
                cell_i = cell_own;
                if (local_cells[cell_i].size() > 1)
                    visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);

                // left
                if (fequal(par_y, y_mid_cell) && fequal(par_x, cell_x*cell_size)) {
                    // printf("%s\n", "particle at y_mid_cell left edge");
                    // if having left neighbor
                    if (cell_x > 0){
                        cell_i = cell_own - 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                }
                
                // // right
                if (fequal(par_y, y_mid_cell) && fequal(par_x, cell_x*cell_size + cell_size)) {
                    // printf("%s\n", "particle at y_mid_cell right edge");
                    // if having right neighbor
                    if (cell_x < (grid_width-1)) {
                        cell_i = cell_own + 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                }

                // bottom
                if (fequal(par_x, x_mid_cell) && fequal(par_y, cell_y*cell_size)) {
                    // printf("%s\n", "particle at x_mid_cell bottom edge");
                    // if having bottom neighbor
                    if (cell_y > 0) {
                        cell_i = cell_own - grid_width;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                }

                // top
                if (fequal(par_x, x_mid_cell) && fequal(par_y, cell_y*cell_size + cell_size)) {
                    // printf("%s\n", "particle at x_mid_cell top edge");
                    // if having top neighbor
                    if (cell_y < (grid_width-1)) {
                        cell_i = cell_own + grid_width;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                }

                // bottom-left
                if (fless(par_x, x_mid_cell) && fless(par_y, y_mid_cell)) {
                    // left neighbor
                    if (cell_x > 0){
                        cell_i = cell_own - 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                    // bottom neighbor
                    if (cell_y > 0) {
                        cell_i = cell_own - grid_width;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                    // bottom-left neighbor
                    if (cell_x > 0 && cell_y > 0){
                        cell_i = cell_own - grid_width - 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                }

                // bottom-right
                if (fgreater(par_x, x_mid_cell) && fless(par_y, y_mid_cell)) {
                    // right neighbor
                    if (cell_x < (grid_width-1)) {
                        cell_i = cell_own + 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                    // bottom neighbor
                    if (cell_y > 0) {
                        cell_i = cell_own - grid_width;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                    // bottom-right neighbor
                    if (cell_x < (grid_width-1) && cell_y > 0){
                        cell_i = cell_own - grid_width + 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                }

                // top-right
                if (fgreater(par_x, x_mid_cell) && fgreater(par_y, y_mid_cell)) {
                    // right neighbor
                    if (cell_x < (grid_width-1)) {
                        cell_i = cell_own + 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                    // top neighbor
                    if (cell_y < (grid_width-1)) {
                        cell_i = cell_own + grid_width;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                    // top-right neighbor
                    if (cell_x < (grid_width-1) && cell_y < (grid_width-1)) {
                        cell_i = cell_own + grid_width + 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                }

                // top-left
                if (fless(par_x, x_mid_cell) && fgreater(par_y, y_mid_cell)) {
                    // left neighbor
                    if (cell_x > 0){
                        cell_i = cell_own - 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                    // top neighbor
                    if (cell_y < (grid_width-1)) {
                        cell_i = cell_own + grid_width;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                    // top-left neighbor
                    if (cell_x > 0 && cell_y < (grid_width-1)) {
                        cell_i = cell_own + grid_width - 1;
                        if (!local_cells[cell_i].empty())
                            visit_cell_and_apply_force(local_rows[i][j], local_cells, cell_i, &dmin, &davg, &navg);
                    }
                }
            }
        }
        //--------------end apply force------------------//

        //
        // move
        //
        for (int i = 0; i < local_rows.size(); ++i){
            for (int j = 0; j < local_rows[i].size(); ++j){
                move( local_rows[i][j] );
                
                local_rows[i][j].cell_x = (int)floor(local_rows[i][j].x/cell_size);
                local_rows[i][j].cell_y = (int)floor(local_rows[i][j].y/cell_size);
            }
        }


        // gather to rank 0
        //
        // send from local_rows
        for (int j = 0; j < row_total_arr[rank]; ++j)
            if (particles_in_row[local_row_offset + j] >0)
                // MPI_Send(&local_rows[j][0], particles_in_row[local_row_offset + j], PARTICLE, 0, 103, MPI_COMM_WORLD);
                MPI_Isend(&local_rows[j][0], particles_in_row[local_row_offset + j], PARTICLE, 0, 103, MPI_COMM_WORLD, &request);

        // assemble at rank 0
        if (rank==0){
            int n_copied_particles = 0;
            // printf("rank 0 rows.back().size() size before send: %d\n", rows.back().size());
            for (int i = 0; i < n_proc; ++i){
                int row_offset = row_offset_arr[i];
                    for (int j = 0; j < row_total_arr[i]; ++j){
                        if (particles_in_row[row_offset + j] > 0){
                            // MPI_Recv(&rows[row_offset + j][0], particles_in_row[row_offset + j], PARTICLE, i, 103, MPI_COMM_WORLD, &status);
                            MPI_Recv(&particles[n_copied_particles], particles_in_row[row_offset + j], PARTICLE, i, 103, MPI_COMM_WORLD, &status);
                            // MPI_Recv(&particles_prev[n_copied_particles], particles_in_row[row_offset + j], PARTICLE, i, 103, MPI_COMM_WORLD, &status);
                            // replace data in particles with those in rows
                            // memmove(particles + (row_offset | j ? 1 : 0)*particles_in_row[row_offset + j-1], &rows[row_offset + j][0], particles_in_row[row_offset + j]*sizeof(particle_t));
                            // memcpy(&particles_prev[n_copied_particles], &rows[row_offset + j][0], particles_in_row[row_offset + j] * sizeof(particle_t));
                            n_copied_particles += particles_in_row[row_offset + j];
                        }
                    }
            }
            // memcpy (&particles[0], &particles_prev[0], n * sizeof(particle_t));
            // reset the main rows vector
            for (int i = 0; i < grid_width; ++i)
                std::vector<particle_t>().swap(rows[i]);
            // for (int i=0; i<n; ++i)
            //     rows[particles[i].cell_y].push_back(particles[i]);
        } 








        if( find_option( argc, argv, "-no" ) == -1 )
        {
          
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

 
          if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }

        // MPI_Barrier(MPI_COMM_WORLD); 
    } //--------------end MAIN sim------------------//
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, number of nodes %d, simulation time = %g seconds", n, n_proc, simulation_time);

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
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    free(particles_in_row);
    free(row_total_arr);
    free(row_offset_arr);
    // free(cell_total_arr);
    // free(cell_offset_arr);
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    ////////////////////////////////////////////////////////////////////

    free( particles );
    // free( particles_prev );
    
    return 0;
}
