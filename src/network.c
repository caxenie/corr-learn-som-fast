#include "simulation.h"

/* initialize a neural population */
population init_population(short idx, int psize)
{
	population p;
	double sigma_def = 0.045000f;
	double sumWcross = 0.0f;

	p.id = idx;
	p.size = psize;	

	p.Winput = (double*)malloc(p.size*sizeof(population));
	p.Wcross = (double**)calloc(p.size, sizeof(population*));
	for(int i = 0; i<p.size; i++)
		p.Wcross[i] = (double*)calloc(p.size, sizeof(population));
	for(int i =0; i<p.size; i++){
		for(int j = 0; j<p.size; j++){
			p.Wcross[i][j] = (double)rand()/(double)RAND_MAX;
			sumWcross += p.Wcross[i][j];
		}
	}
	for(int i =0; i<p.size; i++){
		for(int j = 0; j<p.size; j++){
			p.Wcross[i][j] /= sumWcross;
		}
	}

 	p.s = (double*)calloc(p.size, sizeof(double));
	for (int i=0; i<p.size;i++)
		p.s[i] = sigma_def;

	p.a = (double*)calloc(p.size, sizeof(double));

	return p;	
}


/* initialize the network */
network* init_network(int npop, int psz)
{
	network* n = (network*)calloc(1, sizeof(network));
	n->nsize = npop;
	n->pops = (population*)calloc(n->nsize, sizeof(population));
	for (int i = 0; i<n->nsize; i++)
		n->pops[i] = init_population((short)i, psz);

	return n;
}

/* deallocate a neural population */
void deinit_population(population *p)
{
	for(int i=0; i<p->size;i++){
		free(p->Wcross[i]);	
	}
	free(p->Winput);
	free(p->Wcross);	
	free(p->s);
	free(p->a);
	free(p);
}


/* deallocate a network */
void deinit_network(network *net)
{
	free(net->pops);
	free(net);
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

/* number of shuffles for maps ids for cross-modal circular permutation in Hebbian learning rule */
long num_shuffles(int n, int r)
{
    long f[n + 1];
    f[0]=1;
    for (int i=1;i<=n;i++)
        f[i]=i*f[i-1];
    return f[n]/f[r]/f[n-r];
}

/* shuffle the maps ids for cross-modal circular permutation in Hebbian learning rule */
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

/* decode population to real-world value */
double decode_population(network*n, int pre_id, int post_id, double init_cond)
{
	double dec_val = init_cond;
	int num_iters_opt = 2000;
	double precision = 1e-6;
	double fx = 0.0f, temp_fx = 0.0f, fx_ant = 0.0f;
	double dfx = 0.0f;
	double* dir_act = (double*)calloc(n->pops[post_id].size, sizeof(double));
	double* ind_act = n->pops[post_id].a;
	double tot_act = 0.0f;
	double* cur_act = (double*)calloc(n->pops[post_id].size, sizeof(double));

	printf("ai\n");
	for(int v = 0; v< n->pops[post_id].size; v++){	
			printf(" %lf  \n", ind_act[v]);
		}

	/* use Newton-Raphson optimization to find precise decoded value */
	for(int opt_iter = 0; opt_iter < num_iters_opt; opt_iter++){
		/* compute direct activation given the optimized variable */
		for(int i=0; i<n->pops[post_id].size; i++){
				cur_act[i] = (1/(sqrt(2*M_PI)*n->pops[post_id].s[i]))*
                                             exp(-pow((dec_val - n->pops[post_id].Winput[i]),2)/
					     (2*pow(n->pops[post_id].s[i], 2)));
				
		} 
		/* normalization routine in the optimization process */
		for(int snid = 0; snid<n->pops[post_id].size; snid++)
                              tot_act += cur_act[snid];
                for(int snid = 0; snid<n->pops[post_id].size; snid++)
                              cur_act[snid] /= tot_act;
                /* update the activity for next iteration */
		dir_act = cur_act;
		
		printf("ai | ad\n");
		for(int v = 0; v< n->pops[post_id].size; v++){	
			printf(" %lf  | % lf \n", ind_act[v], dir_act[v]);
		}
   	    	/* function to optimize is the error between the direct and indirect activation */
		for(int i=0;i<n->pops[post_id].size;i++)
			temp_fx += pow(ind_act[i] - dir_act[i], 2);
		fx = sqrt(temp_fx);
		/* ... and its derivative computed numerically */
		dfx = (fx - fx_ant)/2;
		/* update the estimate of the value */
		dec_val -= (fx/dfx);
		/* check convergence */
		if(fabs(fx/dfx) < precision){ printf("precision reached\n"); break;}
		/* update history */
		fx_ant = fx;
		temp_fx = 0.0f;
		tot_act = 0.0f;
	}
	return dec_val;
}
