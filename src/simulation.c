#include "simulation.h"

/* initialize simulation */
simulation* init_simulation(int nepochs, network*net)
{
	simulation* s = (simulation*)calloc(1, sizeof(simulation));
	double A = 0.0f, B = 0.0f; 

	s->max_epochs = nepochs;
	s->t0 = 0;
	s->tf_lrn_in = s->max_epochs/4;
	s->tf_lrn_cross = s->max_epochs;
	
	s->alpha = (double*)calloc(s->tf_lrn_in, sizeof(double));
	s->sigma = (double*)calloc(s->tf_lrn_in, sizeof(double));
	s->eta = (double*)calloc(s->tf_lrn_cross, sizeof(double));
	s->xi = (double*)calloc(s->tf_lrn_cross, sizeof(double));
	s->alpha = parametrize_process(ALPHAI, ALPHAF, s->t0, s->tf_lrn_in, INVTIME);
	s->sigma = parametrize_process((net->pops->size)/3, SIGMAF, s->t0, s->tf_lrn_in, INVTIME);
	for (int i = 0;i<s->tf_lrn_cross;i++)
		s->eta[i] = ETA;
	for (int i = 0;i<s->tf_lrn_cross;i++)
		s->xi[i] = XI;
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
	outdata *runtime = (outdata*)calloc(1, sizeof(outdata));

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
					cur_act = (double*)calloc(s->n->pops[0].size, sizeof(double));
					/* loop through neurons in current population */
					for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
						/* compute sensory elicited activation */
						cur_act[nidx] = (1/(sqrt(2*M_PI)*s->n->pops[pidx].s[nidx]))*
							   exp(-pow((insample - s->n->pops[pidx].Winput[nidx]),2)/(2*pow(s->n->pops[pidx].s[nidx], 2)));
					}	
					/* normalize the activity vector of the population */
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++)
						tot_act	+= cur_act[snid];
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++)	
						cur_act[snid] /= tot_act;
					/* update the activity for next iteration */
					for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++)
						s->n->pops[pidx].a[nidx] = (1-s->eta[tidx])*s->n->pops[pidx].a[nidx] + s->eta[tidx]*cur_act[nidx];
					/* competition step - find the neuron with maximum activity */
					/* find the neuron with maximum activity and its index in the population */
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++){ 
						if(s->n->pops[pidx].a[snid] > win_act){
							win_act = s->n->pops[pidx].a[snid];
							win_idx = snid;
						}
					}
					for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
						/* compute the neighborhood kernel */
						if(!WRAP_POP)
							/* if we do not treat boundary effects use traditional neighborhood function */
							hwi[nidx] = exp(-pow(fabs(nidx -  win_idx), 2)/(2*pow(s->sigma[tidx], 2)));
						else
							/* wrap up the population to avoid boundary effects */
				                        /* dist = min{|i-j|, N - |i-j|} */
			        	                hwi[nidx] = exp(-pow(MIN(fabs(nidx-win_idx), s->n->pops[pidx].size-fabs(nidx-win_idx)), 2)/
								       (2*pow(s->sigma[tidx], 2)));
						/* compute the sensory input synaptic weight */
						s->n->pops[pidx].Winput[nidx] += s->alpha[tidx]*hwi[nidx]*(insample - s->n->pops[pidx].Winput[nidx]); 
						/* update the shape of the tuning curve for the current neuron */
						if(!ASYMM_FUNC)
							s->n->pops[pidx].s[nidx] += s->alpha[tidx]*
										    exp(-pow(fabs(nidx -  win_idx), 2)/(2*pow(s->sigma[tidx], 2)))*
										    (pow((insample - s->n->pops[pidx].Winput[nidx]) , 2) - pow(s->n->pops[pidx].s[nidx], 2));
						else
							s->n->pops[pidx].s[nidx] += s->alpha[tidx]*0.005*
                                                                        	    exp(-pow(fabs(nidx -  win_idx), 2)/(2*pow(s->sigma[tidx], 2)))*
                                        	                                    (pow((insample - s->n->pops[pidx].Winput[nidx]) , 2) - pow(s->n->pops[pidx].s[nidx], 2));
					}/* end for each neuron in the population */
				    }/* end loop through populations */	
				}/* end loop of sensory data presentation */
			}/* end loop for training input data distribution */

			/* cross-modal learning loop */
                        for(int didx = 0; didx < in->len; didx++){
				/* use the learned sensory elicited synaptic weights and compute activation for each population */
				/* loop through populations */
                                for(int pidx = 0; pidx < s->n->nsize; pidx++){
                                        tot_act = 0.0f;
					cur_act = (double*)calloc(s->n->pops[0].size, sizeof(double));
                                        insample = in->data[didx][pidx];
                                        /* loop through neurons in current population */
                                        for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++){
                                                /* compute sensory elicited activation */
                                                cur_act[nidx] = (1/(sqrt(2*M_PI)*s->n->pops[pidx].s[nidx]))*
                                                           exp(-pow((insample - s->n->pops[pidx].Winput[nidx]),2)/(2*pow(s->n->pops[pidx].s[nidx], 2)));
                                        }
                                        /* normalize the activity vector of the population */
                                        for(int snid = 0; snid<s->n->pops[pidx].size; snid++)
                                                tot_act += cur_act[snid];
					for(int snid = 0; snid<s->n->pops[pidx].size; snid++)
	                                        cur_act[snid] /= tot_act;
					/* update the activity for next iteration */
                                        for(int nidx = 0; nidx<s->n->pops[pidx].size;nidx++)
                                                s->n->pops[pidx].a[nidx] = (1-s->eta[tidx])*s->n->pops[pidx].a[nidx] + s->eta[tidx]*cur_act[nidx];
				 }/* end loop for each population */
				        
					/* apply the learning rule */
					int pidx = 0;
					/* utils for pops ids shuffling */
					unsigned int base_idx[] = {0,1,2};
					const unsigned int pre_post_pair = 2;
					switch(LEARNING_RULE){	
						case HEBB:
							/* cross-modal hebbian learning */
							do{
								/* shuffle the populations ids for updating */
								int* cur_ids = base_idx;
								printf("%d %d\n", cur_ids[0], cur_ids[1]);
								for(int i=0;i<s->n->pops[pidx].size;i++){
								  for(int j=0; j<s->n->pops[pidx].size; j++){
									s->n->pops[cur_ids[0]].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[cur_ids[0]].Wcross[i][j]+
													s->xi[tidx]*s->n->pops[cur_ids[0]].a[i]*s->n->pops[cur_ids[1]].a[j];
								  }
							        }
							}while(shuffle_pops_ids(base_idx, s->n->nsize, pre_post_pair));
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
						do{
                                                    /* shuffle the populations ids for updating */
                                                    int* cur_ids = base_idx;
						    printf("%d %d\n", cur_ids[0], cur_ids[1]);
						    for(int i=0;i<s->n->pops[pidx].size;i++){
                                	                        for(int j=0; j<s->n->pops[pidx].size; j++){
							              s->n->pops[cur_ids[0]].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[cur_ids[0]].Wcross[i][j]+
                                                                                                    s->xi[tidx]*
												    (s->n->pops[cur_ids[0]].a[i] - avg_act[cur_ids[0]][i])*
												    (s->n->pops[cur_ids[1]].a[j] - avg_act[cur_ids[1]][j]);
                        	                                }
                	                                }
						}while(shuffle_pops_ids(base_idx, s->n->nsize, pre_post_pair));
						break;
						case OJA:
						    /* compute the global synaptic strength in the likage */
                                                do{     
                                                    /* shuffle the populations ids for updating */
                                                    int* cur_ids = base_idx;
						    	printf("%d %d\n", cur_ids[0], cur_ids[1]);
							for(int i=0;i<s->n->pops[pidx].size;i++){
                                                                for(int j=0; j<s->n->pops[pidx].size; j++){
                                                                        s->n->pops[cur_ids[0]].Wcross[i][j] = (1-s->xi[tidx])*s->n->pops[cur_ids[0]].Wcross[i][j]+
                                                                                                        s->xi[tidx]*s->n->pops[cur_ids[0]].a[i]*s->n->pops[cur_ids[1]].a[j];
									sum_wcross[cur_ids[0]] += s->n->pops[cur_ids[0]].Wcross[i][j];
                                                                }
                                                        }
						}while(shuffle_pops_ids(base_idx, s->n->nsize, pre_post_pair));
				
						    /* compute the synaptic weights using Oja's local normalizing PCA rule */
                                                do{
                                                    /* shuffle the populations ids for updating */
                                                    int* cur_ids = base_idx;
							for(int i=0;i<s->n->pops[pidx].size;i++){
                                                                for(int j=0; j<s->n->pops[pidx].size; j++){
                                                                        s->n->pops[cur_ids[0]].Wcross[i][j] = ((1-s->xi[tidx])*s->n->pops[cur_ids[0]].Wcross[i][j]+
                                                                                                        s->xi[tidx]*s->n->pops[cur_ids[0]].a[i]*s->n->pops[cur_ids[1]].a[j])/
			                                                                        	sqrt(sum_wcross[cur_ids[0]]);
                                                                }
                                                        }
						}while(shuffle_pops_ids(base_idx, s->n->nsize, pre_post_pair));
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
	/* fill in the return struct */
	runtime->in = in;
	runtime->sim = s;
	
	return runtime;
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

/* dump the runtime data to file on disk */
char* dump_runtime_data(outdata *od)
{
	time_t rawt; time(&rawt);
        struct tm* tinfo = localtime(&rawt);
        FILE* fout;
        char* nfout = (char*)calloc(400, sizeof(char));
        char simparams[200];

        strftime(nfout, 150, "%Y-%m-%d__%H:%M:%S", tinfo);
        strcat(nfout, "_cln_runtime_data_");
        sprintf(simparams, "%d_epochs_%d_populations_%d_neurons",
                           od->sim->max_epochs,
                           od->sim->n->nsize,
                           od->sim->n->pops[od->sim->n->nsize-1].size);
        strcat(nfout, simparams);

        if((fout = fopen(nfout, "wb"))==NULL){
                printf("dump_runtime_data: Cannot create output file.\n");
                return NULL;
        }
        fwrite(od, sizeof(outdata), 1, fout);
        fclose(fout);
        return nfout;
}

/* dump the runtime data to file on disk - explicit sequential write */
char* dump_runtime_data_extended(outdata *od)
{
	time_t rawt; time(&rawt);
        struct tm* tinfo = localtime(&rawt);
        FILE* fout;
        char* nfout = (char*)calloc(400, sizeof(char));
        char simparams[200];

        strftime(nfout, 150, "%Y-%m-%d__%H:%M:%S", tinfo);
        strcat(nfout, "_cln_extended_runtime_data_");
        sprintf(simparams, "%d_epochs_%d_populations_%d_neurons",
                           od->sim->max_epochs,
                           od->sim->n->nsize,
                           od->sim->n->pops[od->sim->n->nsize-1].size);
        strcat(nfout, simparams);

        if((fout = fopen(nfout, "wb"))==NULL){
                printf("dump_runtime_data: Cannot create output file.\n");
                return NULL;
        }
	/* in order ot read the data in Matlab write explicitly every field */
        /* simulation data */
	fwrite(&(od->sim->max_epochs), sizeof(int), 1, fout);
	fwrite(&(od->sim->t0), sizeof(int), 1, fout);
	fwrite(&(od->sim->tf_lrn_in), sizeof(int), 1, fout);
	fwrite(&(od->sim->tf_lrn_cross), sizeof(int), 1, fout);
	for(int i = 0; i<od->sim->tf_lrn_in; i++)
		fwrite(&(od->sim->alpha[i]), sizeof(double), 1, fout);
	for(int i = 0; i<od->sim->tf_lrn_in; i++)
		fwrite(&(od->sim->sigma[i]), sizeof(double), 1, fout);
	for(int i = 0; i<od->sim->tf_lrn_cross; i++)
		fwrite(&(od->sim->eta[i]), sizeof(double), 1, fout);
	for(int i = 0; i<od->sim->tf_lrn_cross; i++)
		fwrite(&(od->sim->xi[i]), sizeof(double), 1, fout);
	/* network data */
	fwrite(&(od->sim->n->nsize), sizeof(short), 1, fout);
	for(int pidx = 0; pidx<od->sim->n->nsize;pidx++){
		fwrite(&(od->sim->n->pops[pidx].id), sizeof(short), 1, fout);
		fwrite(&(od->sim->n->pops[pidx].size), sizeof(int), 1, fout);
		for(int i=0; i<od->sim->n->pops[pidx].size;i++){
			fwrite(&(od->sim->n->pops[pidx].Winput[i]), sizeof(double), 1, fout);
			for(int j = 0;j < od->sim->n->pops[pidx].size; j++){
				fwrite(&(od->sim->n->pops[pidx].Wcross[i][j]), sizeof(double), 1, fout);
			}
			fwrite(&(od->sim->n->pops[pidx].s[i]), sizeof(double), 1, fout);
			fwrite(&(od->sim->n->pops[pidx].a[i]), sizeof(double), 1, fout);
		}
	}
	/* input data */
	fwrite(&(od->in->npop), sizeof(int), 1, fout);
	fwrite(&(od->in->popsize), sizeof(int), 1, fout);
	fwrite(&(od->in->len), sizeof(int), 1, fout);
	for(int i=0;i<od->in->len;i++){
		for(int j=0;j<od->in->npop;j++){
			fwrite(&(od->in->data[i][j]), sizeof(double), 1, fout);	
		}
	}
        fclose(fout);
        return nfout;
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

unsigned int shuffle_pops_ids(unsigned int *ar, size_t n, unsigned int k)
{
    unsigned int finished = 0;
    unsigned int changed = 0;
    unsigned int i;

    if (k > 0) {
        for (i = k - 1; !finished && !changed; i--) {
            if (ar[i] < (n - 1) - (k - 1) + i) {
                /* Increment this element */
                ar[i]++;
                if (i < k - 1) {
                    /* Turn the elements after it into a linear sequence */
                    unsigned int j;
                    for (j = i + 1; j < k; j++) {
                        ar[j] = ar[j - 1] + 1;
                    }
                }
                changed = 1;
            }
            finished = i == 0;
        }
        if (!changed) {
            /* Reset to first combination */
            for (i = 0; i < k; i++) {
                ar[i] = i;
            }
        }
    }
    return changed;
}


