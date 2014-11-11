#include "data.h"

/* read input data from file and populate struct */
indata* generate_input_data(int np, int psz, int l)
{
	indata *id = (indata*)calloc(1, sizeof(indata));
	double *base_var = generate_rnd_vector(DIST_TYPE, RANGE, l, NU_DIST_MODEL);
	double **rel_vars = (double**)calloc(l, sizeof(double));

	id->data = (double**)calloc(l, sizeof(double*));
	
  	id->npop = np;
	id->popsize = psz;
	id->len = l;
	
	for(int i = 0;i<id->len; i++)
		id->data[i] = (double*)calloc(id->npop, sizeof(double));	

	for(int i=0;i<id->len;i++)
		rel_vars[i] = (double*)calloc(id->npop-1, sizeof(double));

	/* HERE EMBED THE RELATIONS */
	for(int i = 0;i<id->len;i++){
		for(int j = 0;j<id->npop;j++){
			rel_vars[i][j] = pow(base_var[i], 3);
		}
	}
        for(int i = 0;i<id->len;i++){
               for(int j = 0;j<id->npop;j++){
			if(j==0) id->data[i][j] = base_var[i];
			else id->data[i][j] = rel_vars[i][j];	
		}
	}
	
	return id;
}	

