#include <p4est.h>

/*

x=-l     x=-d-lc  x=-d     x=d      x=d+a    x=l

o--------o--------o--------o--------o--------o  y=domain_height
|        |        |        |        |        | 
|   10   |   11   |   12   |   13   |    14  |
o--------o--------o--------o--------o--------o  y=h2
|   5    |   6    |   7    |   8    |    9   |
o--------o--------o--------o--------o--------o  y=h1
|   0    |   1    |   2    |   3    |    4   |
o--------o--------o--------o--------o--------o  y=0

*/

constexpr double domain_width  = 4.0;
constexpr double l             = domain_width;
constexpr double d             = .2;
constexpr double lc            = .1;
constexpr double a             = .5;
constexpr double domain_height = 4.0;
constexpr double h3            = domain_height;
constexpr double h1            = 1.0e-3;
constexpr double h2            = 5.0e-2;
constexpr p4est_topidx_t simple_conn_num_vertices = 24;
constexpr p4est_topidx_t simple_conn_num_trees    = 15;
constexpr double simple_conn_p[simple_conn_num_vertices*2] =
  {
    -l, 0.0,      -d-lc, 0.0,     -d, 0.0,      d,0.0,      d+a, 0.0,      l, 0.0,
    -l, h1,       -d-lc, h1,      -d, h1,       d,h1,       d+a, h1,       l, h1,
    -l, h2,       -d-lc, h2,      -d, h2,       d,h2,       d+a, h2,       l, h2,
    -l, h3,       -d-lc, h3,      -d, h3,       d,h3,       d+a, h3,       l, h3,
  };

constexpr p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {
    
    1, 2,  8,  7, 1,
    2, 3,  9,  8, 1,
    3, 4, 10,  9, 1,
    4, 5, 11, 10, 1,
    5, 6, 12, 11, 1,

    6+1, 6+2, 6+ 8, 6+ 7, 6+1,
    6+2, 6+3, 6+ 9, 6+ 8, 6+1,
    6+3, 6+4, 6+10, 6+ 9, 6+1,
    6+4, 6+5, 6+11, 6+10, 6+1,
    6+5, 6+6, 6+12, 6+11, 6+1,

    6+6+1, 6+6+2, 6+6+ 8, 6+6+ 7, 6+6+1,
    6+6+2, 6+6+3, 6+6+ 9, 6+6+ 8, 6+6+1,
    6+6+3, 6+6+4, 6+6+10, 6+6+ 9, 6+6+1,
    6+6+4, 6+6+5, 6+6+11, 6+6+10, 6+6+1,
    6+6+5, 6+6+6, 6+6+12, 6+6+11, 6+6+1

  };

static constexpr unsigned refine_steps = 4;


constexpr double kG = 8.8e-12;
constexpr double kS = kG/10;
