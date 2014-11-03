#include "simulation.h"

/* initialize simulation */
simulation* init_simulation(int nepochs, network*net)
{
	simulation* s = (simulation*)calloc(1, sizeof(simulation));
	double A = 0.0f, B = 0.0f; 

	s->max_epochs = nepochs;
	s->t0 = 0;
	s->tf_lrn_in = s->max_epochs/10;
	s->tf_lrn_cross = s->max_epochs;
	
	s->alpha = (double*)calloc(s->tf_lrn_in, sizeof(double));
	s->sigma = (double*)calloc(s->tf_lrn_in, sizeof(double));
	s->eta = (double*)calloc(s->tf_lrn_cross, sizeof(double));
	s->xi = (double*)calloc(s->tf_lrn_cross, sizeof(double));
	B = (ALPHAF*s->tf_lrn_in - ALPHAI*s->t0)/(ALPHAI - ALPHAF);
	A = ALPHAI*s->t0 + B*ALPHAI;
	for (int i = 0;i<s->tf_lrn_in;i++){
		s->alpha[i] = A/(i+B);
	}
	B = (SIGMAF*s->tf_lrn_in - ((net->pops[1]->size)/10)*s->t0)/(((net->pops[1]->size)/10)- SIGMAF);
	A =  ((net->pops[1]->size)/10)*s->t0 + B*((net->pops[1]->size)/10);
	for (int i = 0;i<s->tf_lrn_in;i++){
		s->sigma[i] = A/(i+B);
	}
	for (int i = 0;i<s->tf_lrn_cross;i++){
		s->eta[i] = ETA;
	}
	for (int i = 0;i<s->tf_lrn_cross;i++){
		s->xi[i] = XI;
	}
	s->n = net;
	return s;
}

/* destroy the simulation */
void deinit_simulation(simulation* s)
{
	free(s->alpha);
	free(s->sigma);	
	free(s->eta);
	free(s->xi);
	deinit_network(s->n);
}

