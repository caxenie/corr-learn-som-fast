#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

/* data distribution type */
enum{
        UNIFORM = 0,
        NONUNIFORM
};

/* non-uniform distribution type */
enum{
        INCPOWERLAW = 0,
        DECPOWERLAW,
        GAUSS,
        CONVEX
};

/* relation type */
enum{
	LINEAR = 1,
	ORDER2, 
	ORDER3,
	SINE,
	SQUARE,
	COMPLEX
};

#define DIST_TYPE UNIFORM
#define NU_DIST_MODEL GAUSS
#define RANGE 1

/* floating point radom number generator between 0 and 1 */
double randf();
/* non-uniform random distribution generator */
double * generate_rnd_vector(int type, int range, int numv, int dist);
/* Normal random variate vector generator, of mean m and std. dev. s */
double* gauss_rand(int num_vals, double mean, double std_dev);

