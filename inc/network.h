#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <float.h>

/* neural population definition */
typedef struct{
	short id;	/* id of the population */
	int size;	/* size of the population */
	double **Winput;/* sensory afferents synaptic connections */
	double **Wcross;/* cross modal afferents synaptic connections */
	double *s;	/* population tuning curves shapes */
	double *a;	/* population activation */
}population;

/* network definition */
typedef struct{
	short nsize;		/* number of populations in the net */
	population **pops;	/* populations in the net */
}network;

/* initialize a neural population */
population* init_population(short idx, int psize);
/* initialize the network */
network* init_network(int npop, int psize);

/* deallocate a neural population */
void deinit_populaion(population *pop);
/* deallocate a network */
void deinit_network(network *net);

