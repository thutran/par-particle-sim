#include <vector>
#include <cmath>
#include <limits>

#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }
inline bool fequal (const double& a, const double& b){ return fabs(a-b) <= std::numeric_limits<double>::epsilon(); }
inline bool fgreater (const double& a, const double& b){ return (a-b) > std::numeric_limits<double>::epsilon(); }
inline bool fless (const double& a, const double& b){ return (b-a) > std::numeric_limits<double>::epsilon(); }

//
//  saving parameters
//
const int NSTEPS = 200;
// const int NSTEPS = 100;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;

  int cell_x;
  int cell_y;
} particle_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );



//----------------------------------------------
//
// get constants
//
// get bounds of x and y
double get_size();
// cutoff value
double get_cutoff(); 

//
// grid
//
// initialize grid
// void init_grid(const int& n, particle_t *p, std::vector<particle_t*> *cells);
void init_grid(const int& n, particle_t* p, const double& cell_size, const int& grid_width, std::vector<std::vector<particle_t*>> &cells);
// void refesh_grid(const int& n, particle_t *p, std::vector<particle_t*> *cells);

void visit_cell_and_apply_force(particle_t &particle, std::vector<std::vector<particle_t*>> &cells, const int& cell_i, double *dmin, double *davg, int *navg);

void move_to_cell(particle_t &p, const double& cell_size, const int& grid_width, std::vector<std::vector<particle_t*>> &cells);

// put particles into corresponding cells (MPI version)
void fill_grid(const int& p_array_size, particle_t* p, const int& grid_width, std::vector<std::vector<particle_t*>> &cells);

#endif


