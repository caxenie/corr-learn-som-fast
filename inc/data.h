#include "simulation.h"
#include "tools.h"

/* input data struct */
typedef struct{
	int npop;		/* number of populations in the net */
	int popsize;		/* size of population: num of neurons */
	int len;		/* length of the training dataset */
	double **data;   	/* actual data */
}indata;

/* output data struct */
typedef struct{
	simulation *sim;	/* simulation params */
	indata *in;		/* input data */
	network *n;		/* network */
}outdata;

/* generate input data and populate struct */
indata* generate_input_data(int np, int psz, int l);

/* dump runtime data to file */
int write_output_data(outdata *odata);


