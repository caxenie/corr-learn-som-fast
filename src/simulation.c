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
	double *cur_act = (double*)calloc(s->n->pops[0].size, sizeof(double));
	double **avg_act = (double**)calloc(s->n->nsize, sizeof(double*));
	for(int i=0; i<s->n->nsize; i++)
		avg_act[i] = (double*)calloc(s->n->pops[0].size, sizeof(double));
	double *sum_wcross = (double*)calloc(s->n->nsize, sizeof(double));
	double *max_wcross = (double*)calloc(s->n->nsize, sizeof(double));
	double win_act = 0.0f;
	int win_idx = 0;
	double omega = 0.0f;
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
					win_act = 0.0f;
					win_idx = 0.0f;
					hwi = (double*)calloc(s->n->pops[0].size, sizeof(double));
					insample = in->data[didx][pidx];
					/* loop through neurons in current population */
					for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
						/* compute sensory elicited activation */
						cur_act[nidx] = (1/(sqrt(2*M_PI)*s->n->pops[pidx].s[nidx]))*
							   exp(-pow((insample - s->n->pops[pidx].Winput[nidx]),2)/2*pow(s->n->pops[pidx].s[nidx], 2));
					}	
		                        printf("BP1 cur_act:%lf\n", cur_act[34]);
		                        printf("BP1 avg_act:%lf\n", avg_act[0][34]);
	
					/* normalize the activity vector of the population */
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++){
						tot_act	+= s->n->pops[pidx].a[snid];
					}	
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++){	
						cur_act[snid] /= tot_act;
					}
					printf("BP2 cur_act:%lf\n", cur_act[34]);
                                        printf("BP2 avg_act:%lf\n", avg_act[0][34]);
					/* update the activity for next iteration */
					for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
						s->n->pops[pidx].a[nidx] = (1-s->eta[tidx])*s->n->pops[pidx].a[nidx] + s->eta[tidx]*cur_act[nidx];
					}
					printf("BP1 cur_act:%lf\n", s->n->pops[0].a[34]);
                                        printf("BP1 avg_act:%lf\n", avg_act[0][34]);
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
						/* update the shape of the tuning curve for the current neuron */
						s->n->pops[pidx].s[nidx] += s->alpha[tidx]*(1/(sqrt(2*M_PI)*s->sigma[tidx]))*
									    exp(-pow(fabs(nidx -  win_idx), 2)/(2*pow(s->sigma[tidx], 2)))*
									    (pow((insample - s->n->pops[pidx].Winput[nidx]) , 2) - pow(s->n->pops[pidx].s[nidx], 2));
						
					}/* end for each neuron in the population */
				    }/* end loop through populations */	
				}/* end loop of sensory data presentation */
			}/* end loop for training input data distribution */
	
		        printf("BP2 cur_act:%lf\n", s->n->pops[0].a[34]);
		        printf("BP2 avg_act:%lf\n", avg_act[0][34]);

			/* cross-modal learning loop */
                        for(int didx = 0; didx < in->len; didx++){
				/* use the learned sensory elicited synaptic weights and compute activation for each population */
				/* loop through populations */
                                for(int pidx = 0; pidx < s->n->nsize; pidx++){
                                        tot_act = 0.0f;
                                        insample = in->data[didx][pidx];
                                        /* loop through neurons in current population */
                                        for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
                                                /* compute sensory elicited activation */
                                                cur_act[nidx] = (1/(sqrt(2*M_PI)*s->n->pops[pidx].s[nidx]))*
                                                           exp(-pow((insample - s->n->pops[pidx].Winput[nidx]),2)/2*pow(s->n->pops[pidx].s[nidx], 2));
                                        }
                                        /* normalize the activity vector of the population */
                                        for(int snid = 0; snid<s->n->pops[pidx].size; snid++){
                                                tot_act += s->n->pops[pidx].a[snid];
                                        }
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++){
	                                        cur_act[snid] /= tot_act;
					}
					/* update the activity for next iteration */
                                        for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
                                                s->n->pops[pidx].a[nidx] = (1-s->eta[tidx])*s->n->pops[pidx].a[nidx] + s->eta[tidx]*cur_act[nidx];
                                        }
				 }/* end loop for each population */
					/* apply the learning rule */
					int pidx = 0;
					switch(LEARNING_RULE){	
						case HEBB:
							/* cross-modal hebbian learning */
							for(int i=0;i<s->n->pops[pidx].size;i++){
								for(int j=0; j<s->n->pops[pidx].size; j++){
									s->n->pops[0].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[0].Wcross[i][j]+
													s->xi[tidx]*s->n->pops[0].a[i]*s->n->pops[1].a[j];
								}
							}
							for(int i=0;i<s->n->pops[pidx].size;i++){
                                                                for(int j=0; j<s->n->pops[pidx].size; j++){
                                                                        s->n->pops[1].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[1].Wcross[i][j]+     
                                                                                                        s->xi[tidx]*s->n->pops[1].a[i]*s->n->pops[0].a[j];
                                                                }
                                                        }
						break;
						case COVARIANCE:
							/* compute the mean value decay */
							omega = 0.002f + 0.998f/(tidx+2);
							/* compute the average activity */
							for(int pidx = 0; pidx < s->n->nsize; pidx++){
								for(int snid = 0; snid<s->n->pops[pidx].size; snid++){
                                                                        avg_act[pidx][snid] = (1-omega)*avg_act[pidx][snid] + omega*s->n->pops[pidx].a[snid];
                                                                }
							}
							/* cross-modal covariance learning rule */
							for(int i=0;i<s->n->pops[pidx].size;i++){
                                	                        for(int j=0; j<s->n->pops[pidx].size; j++){
                                                           	      printf("cur_act:%lf\n", s->n->pops[0].a[j]);
				                               	      printf("avg_act:%lf\n", avg_act[0][j]);
							              s->n->pops[0].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[0].Wcross[i][j]+
                                                                                                    s->xi[tidx]*
												    (s->n->pops[0].a[i] - avg_act[0][i])*
												    (s->n->pops[1].a[j] - avg_act[1][j]);
                        	                                }
                	                                }
        	                                        for(int i=0;i<s->n->pops[pidx].size;i++){
                                        	                for(int j=0; j<s->n->pops[pidx].size; j++){
                                                                      s->n->pops[1].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[1].Wcross[i][j]+
                                                                                                   s->xi[tidx]*
												   (s->n->pops[1].a[i] - avg_act[1][i])*
												   (s->n->pops[0].a[j] - avg_act[0][j]);
                                                       	        }
							}
						break;
						case OJA:
							/* compute the global synaptic strength in the likage */
							for(int i=0;i<s->n->pops[pidx].size;i++){
                                                                for(int j=0; j<s->n->pops[pidx].size; j++){
                                                                        s->n->pops[0].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[0].Wcross[i][j]+
                                                                                                        s->xi[tidx]*s->n->pops[0].a[i]*s->n->pops[1].a[j];
									sum_wcross[0] += s->n->pops[0].Wcross[i][j];
                                                                }
                                                        }
                                                        for(int i=0;i<s->n->pops[pidx].size;i++){
                                                                for(int j=0; j<s->n->pops[pidx].size; j++){
                                                                        s->n->pops[1].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[1].Wcross[i][j]+
                                                                                                        s->xi[tidx]*s->n->pops[1].a[i]*s->n->pops[0].a[j];
									sum_wcross[1] += s->n->pops[1].Wcross[i][j];
                                                                }
                                                        }
							/* compute the synaptic weights using Oja's local normalizing PCA rule */
							for(int i=0;i<s->n->pops[pidx].size;i++){
                                                                for(int j=0; j<s->n->pops[pidx].size; j++){
                                                                        s->n->pops[0].Wcross[i][j] = ((1-s->xi[tidx])*s->n->pops[0].Wcross[i][j]+
                                                                                                        s->xi[tidx]*s->n->pops[0].a[i]*s->n->pops[1].a[j])/
			                                                                        	sqrt(sum_wcross[0]);
                                                                }
                                                        }
                                                        for(int i=0;i<s->n->pops[pidx].size;i++){
                                                                for(int j=0; j<s->n->pops[pidx].size; j++){
                                                                        s->n->pops[1].Wcross[i][j] = ((1-s->xi[tidx])*s->n->pops[1].Wcross[i][j]+
                                                                                                        s->xi[tidx]*s->n->pops[1].a[i]*s->n->pops[0].a[j])/
                                                                        				sqrt(sum_wcross[1]);
                                                                }
                                                        }
						break;
					}/* and learning rule selector */
			}/* end dataset loop for cross-modal learning */
	}/* end cross-modal learning process */
	/* normalize cross-modal synaptic weights for visualization */
	for(int pidx = 0; pidx < s->n->nsize; pidx++){
		for(int i=0;i<s->n->pops[pidx].size;i++){
        	        for(int j=0; j<s->n->pops[pidx].size; j++){
				if(s->n->pops[pidx].Wcross[i][j]>max_wcross[pidx])
					max_wcross[pidx] = s->n->pops[pidx].Wcross[i][j];	
			}
		}
	}
	for(int pidx = 0; pidx < s->n->nsize; pidx++){
                for(int i=0;i<s->n->pops[pidx].size;i++){
                        for(int j=0; j<s->n->pops[pidx].size; j++){
				s->n->pops[pidx].Wcross[i][j] /= max_wcross[pidx];
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
