#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include <iostream>

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005
// #define dt      0.005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
    // printf("sqrt(n) %f\t sqrt(density) %f\t size%f\n", sqrt(n), sqrt(density), size);
    // printf("size: %f\n", size);
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        // p[i].x = size*(0.0+(k%sx))/(0+sx);
        // p[i].y = size*(0.0+(k/sx))/(0+sy);
        // printf("x: %f\t y: %f\n", p[i].x, p[i].y);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg){
    // printf("%s\n", "in apply_force");

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // case neighbor is too far away
    if( r2 > cutoff*cutoff ){
        // printf("p_x%f\t p_y%f\t nb_x%f\t nb_y%f\t\n", particle.x, particle.y, neighbor.x, neighbor.y);
        return;
    }
    // case neighbor != particle
	if (r2 != 0){
        // prepare to update min distance between 2 particles
        if (r2/(cutoff*cutoff) < *dmin * (*dmin))
            *dmin = sqrt(r2)/cutoff;
        // prepare to update average distance between 2 particles
        (*davg) += sqrt(r2)/cutoff;
        (*navg) ++;

        // if (*dmin < 0.4)
            // printf("p_x%f\t p_y%f\t nb_x%f\t nb_y%f\t\n", particle.x, particle.y, neighbor.x, neighbor.y);
    }
		
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
	
    //
    //  very simple short-range repulsive force
    //
    // F=ma <--> a=F/m
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;

    // printf("dmin_%f\n", *dmin);
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }

    // p.cell_x = (int)(p.x/cutoff);
    // p.cell_y = (int)(p.y/cutoff);
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}



//----------------------------------------------
//
// get bounds of x and y
double get_size(){
    return size;
}

// get cutoff value
double get_cutoff(){
    return cutoff;
}

// initialize grid
// each particle will be assigned to a cell in the grid
// each cell is a square, size cutoff-cutoff
// grid is a square, size size-size
void init_grid(const int& n, particle_t* p, const double& cell_size, const int& grid_width, std::vector<std::vector<particle_t*>> &cells){
    for (int i = 0; i < grid_width*grid_width; ++i){
        cells.emplace_back(std::vector<particle_t*>());
    }
    for (int i=0; i<n; ++i){
        int cell_x = (int)floor(p[i].x/cell_size);
        int cell_y = (int)floor(p[i].y/cell_size);

        p[i].cell_x = cell_x;
        p[i].cell_y = cell_y;

        cells[cell_y*grid_width + cell_x].push_back(&p[i]);
    }

    // for (int i=0; i<grid_width*grid_width; ++i)
    //     printf("cell_%d size %d\n", i, cells[i].size() );
}

// openmp version
// apply force to all particles in a cell
void visit_cell_and_apply_force(particle_t &particle, std::vector<std::vector<particle_t*>> &cells, const int& cell_i, double *dmin, double *davg, int *navg){
    for (auto it=cells[cell_i].begin(); it!=cells[cell_i].end(); ++it){
        // test if iterator point to the particle itself
        if (*it != &particle){
            // printf("%s\n", "in visit_cell_and_apply_force, call apply_force");
            apply_force(particle, (*(*it)), dmin, davg, navg);

        }
    }
}

// serial version
void move_to_cell(particle_t &p, const double& cell_size, const int& grid_width, std::vector<std::vector<particle_t*>> &cells){
    const int cell_old = p.cell_y*grid_width + p.cell_x;

    move(p);
    p.cell_x = (int)floor(p.x/cell_size);
    p.cell_y = (int)floor(p.y/cell_size);
    
    const int cell_new = p.cell_y*grid_width + p.cell_x;
    
    // change cell
    if (cell_new != cell_old){
        for (int i=0; i<cells[cell_old].size(); ++i){
            if (cells[cell_old].at(i) == &p){
                cells[cell_old].at(i) = cells[cell_old].back();
                cells[cell_old].back() = nullptr;
                cells[cell_old].pop_back();
                break;
            }
        }
        // std::cout << &p  << "\n";
        cells[cell_new].push_back(&p);
    }
}

// put particles into corresponding cells
void fill_grid(const int& p_array_size, particle_t* p, const int& grid_width, std::vector<std::vector<particle_t*>> &cells){
    for (int i=0; i<p_array_size; ++i){
        cells[p[i].cell_y*grid_width + p[i].cell_x].push_back(&p[i]);
    }
    // for (int i=0; i<grid_width*grid_width; ++i)
    //     printf("cell_%d size %d\n", i, cells[i].size() );
}