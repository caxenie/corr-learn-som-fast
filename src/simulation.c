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
	s->alpha = parametrize_process(ALPHAI, ALPHAF, s->t0, s->tf_lrn_in, INVTIME);
	s->sigma = parametrize_process((net->pops->size)/2, SIGMAF, s->t0, s->tf_lrn_in, INVTIME);
	for (int i = 0;i<s->tf_lrn_cross;i++){
		s->eta[i] = ETA;
	}
	for (int i = 0;i<s->tf_lrn_cross;i++){
		s->xi[i] = XI;
	}
	s->n = net;
	return s;
}

/* run simulation and save runtime data struct */
outdata* run_simulation(indata *in, simulation *s)
{
	double insample = 0.0f;
	double tot_act = 0.0f;
	double cur_act = 0.0f;
	double win_act = 0.0f;
	int win_idx = 0;
	double *hwi = (double*)calloc(s->n->pops[0].size, sizeof(double));
	simulation* runtime = (simulation*)calloc(1, sizeof(simulation));
	/* correlation learning loop */
        for(int tidx = s->t0; tidx<s->tf_lrn_cross; tidx++){	
		if(tidx<s->tf_lrn_in){
			/* input distribution learning loop */
			for(int didx = 0; didx < in->len; didx++){
				/* loop through populations */
				for(int pidx = 0; pidx < s->n->nsize; pidx++){
					tot_act = 0.0f;
					cur_act = 0.0f;
					win_act = 0.0f;
					win_idx = 0.0f;
					hwi = (double*)calloc(s->n->pops[0].size, sizeof(double));
					insample = in->data[didx][pidx];
					/* loop through neurons in current population */
					for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
						/* compute sensory elicited activation */
						cur_act = (1/(sqrt(2*M_PI)*s->n->pops[pidx].s[nidx]))*
							   exp(-pow((insample - s->n->pops[pidx].Winput[nidx]),2)/2*pow(s->n->pops[pidx].s[nidx], 2));
					}		
					/* normalize the activity vector of the population */
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++){
						tot_act	+= s->n->pops[pidx].a[snid];
					}	
					cur_act /= tot_act;
					/* update the activity for next iteration */
					for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
						s->n->pops[pidx].a[nidx] = (1-s->eta[tidx])*s->n->pops[pidx].a[nidx] + s->eta[tidx]*cur_act;
					}
					/* competition step - find the neuron with maximum activity */
					/* find the neuron with maximum activity and it's index in the population */
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++){ 
						if(s->n->pops[pidx].a[snid] > win_act){
							win_act = s->n->pops[pidx].a[snid];
							win_idx = snid;
						}
					}
					for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
						/* compute the neighborhood kernel */
						hwi[nidx] = exp(-pow(fabs(nidx -  win_idx), 2)/(2*pow(s->sigma[tidx], 2)));
						/* compute the sensory input synaptic weight */
						s->n->pops[pidx].Winput[nidx] += s->alpha[tidx]*hwi[nidx]*(insample - s->n->pops[pidx].Winput[nidx]); 

					}

				    }	
				}
			}
		}
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

/* parametrize adaptive parameters */
double* parametrize_process(double v0, double vf, int t0, int tf, short type)
{
	int len = tf-t0;
	double* out = (double*)calloc(len, sizeof(double));
	double s = 0.0f, p = 0.0f, A = 0.0f, B = 0.0f;
		
	switch(type){
		case SIGMOID:
			s = -floor(log10(tf))*pow(10, (-(floor(log10(tf)))));		
			p = abs(s*pow(10, (floor(log10(tf))+ floor(log10(tf)/2))));
			for(int i = 0;i<len;i++)
				out[i] = v0 - v0/(1+exp(s*(i-(tf/p)))) + vf;				
		break;
		case INVTIME:
			B = (vf*tf - v0*t0)/(v0-vf);
			A = v0*t0 + B*v0;
			for (int i=0;i<len;i++)
				out[i] = A/(i+B);
		break;
		case EXP:
			if(v0<1) p = -log(v0);
			else p = log(v0);
			for(int i=0;i<len;i++)
				out[i] = v0*exp(-i/(tf/p));
		break;
	}
	return out;
}
