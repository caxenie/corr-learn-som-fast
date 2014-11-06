#include "network.h"
#include "data.h"

/* simulation constant parameters */
#define MAX_EPOCHS      1000
#define N_POP           2
#define POP_SIZE        100
#define DATASET_LEN     1500
#define ALPHAI 0.1f
#define ALPHAF 0.001f
#define SIGMAF 1.0f
#define ETA   1.0f
#define XI    0.001f

/* adaptive processes parametruization types */
enum{
	SIGMOID = 0, 
	INVTIME,
	EXP
};

/* cross-modal learning rules */
enum{
	HEBB = 0, 
	COVARIANCE,
	OJA
};

#define LEARNING_RULE COVARIANCE

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

/* output data struct */
typedef struct{
        simulation *sim;        /* simulation params */
        indata *in;             /* input data */
        network *n;             /* network */
}outdata;


/* initialize simulation */
simulation* init_simulation(int nepochs, network*net);
/* destroy the simulation */
void deinit_simulation(simulation* s);
/* run simulation and save runtime data struct */
outdata* run_simulation(indata *in, simulation *s);

/* parametrize adaptive parameters */
double* parametrize_process(double v0, double vf, int t0, int tf, short type);

