#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

enum{
        UNIFORM = 0,
        NONUNIFORM
};

enum{
        INCPOWERLAW = 0,
        DECPOWERLAW,
        GAUSS,
        CONVEX
};

#define DIST_TYPE UNIFORM
#define NU_DIST_MODEL INCPOWERLAW
#define RANGE 1

/* floating point radom number generator between 0 and 1 */
double randf();
/* non-uniform random distribution generator */
double * generate_rnd_vector(int type, int range, int numv, int dist);
