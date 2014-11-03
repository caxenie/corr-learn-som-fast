#include "network.h"

/* simulation constant parameters */
#define ALPHAI 0.1f
#define ALPHAF 0.001f
#define SIGMAF 1.0f
#define ETA   1.0f
#define XI    0.001f

/* simulation parameters */
typedef struct{
	int max_epochs; 	/* number of epochs to run the network */
	int t0;			/* initial time for simulation */
	int tf_lrn_in;		/* stop time for input learning */
	int tf_lrn_cross;	/* stop time for cross learning */
	double *alpha;		/* values of the learning rate */
	double *sigma;		/* neighborhood kernel size */
 	double *eta;		/* activity decay factor */
	double *xi;		/* cross modal learning rate */ 	
	network *n;		/* network to simulate */
}simulation;


/* initialize simulation */
simulation* init_simulation(int nepochs, network*net);

/* destroy the simulation */
void deinit_simulation(simulation* s);


