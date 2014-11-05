#include "simulation.h"

/* initialize a neural population */
population* init_population(short idx, int psize)
{
	population* p = (population*)calloc(1, sizeof(population));
	double sigma_def = 0.045000f;
	double sumWcross = 0.0f;

	p->id = idx;
	p->size = psize;	

	p->Winput = (double*)malloc(p->size*sizeof(population));
	p->Wcross = (double**)calloc(p->size, sizeof(population*));
	for(int i = 0; i<p->size; i++)
		p->Wcross[i] = (double*)calloc(p->size, sizeof(population));
	for(int i =0; i<p->size; i++){
		for(int j = 0; j<p->size; j++){
			p->Wcross[i][j] = (double)rand()/(double)RAND_MAX;
			sumWcross += p->Wcross[i][j];
		}
	}
	for(int i =0; i<p->size; i++){
		for(int j = 0; j<p->size; j++){
			p->Wcross[i][j] /= sumWcross;
		}
	}

 	p->s = (double*)calloc(p->size, sizeof(double));
	for (int i=0; i<p->size;i++)
		p->s[i] = sigma_def;

	p->a = (double*)malloc(p->size*sizeof(double));

	return p;	
}


/* initialize the network */
network* init_network(int npop, int psz)
{
	network* n = (network*)calloc(1, sizeof(network));
	n->nsize = npop;
	n->pops = (population*)calloc(n->nsize, sizeof(population));
	for (int i = 0; i<n->nsize; i++){
		n->pops = init_population((short)i, psz);
	}
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

