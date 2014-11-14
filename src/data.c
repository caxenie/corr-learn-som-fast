#include "data.h"

/* read input data from file and populate struct */
indata* generate_input_data(int np, int psz, int l, int rtype)
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
	/* check how many variables we haev in the network */
	switch(id->npop){
		case 2:
			for(int i = 0;i<id->len;i++){
				for(int j = 0;j<id->npop;j++){
					switch(rtype){
						case LINEAR:	
							rel_vars[i][j] = RANGE*pow(base_var[i], 1);
						break;
						case ORDER2:
							rel_vars[i][j] = pow(base_var[i], 2);
						break;
						case ORDER3:
							rel_vars[i][j] = pow(base_var[i], 3);
						break;
						case SINE:
							rel_vars[i][j] = sin(base_var[i]);
						break;
					}	
				}		
			}
		break;
		case 3:
		    switch(rtype){
			case COMPLEX:
			for(int i = 0;i<id->len;i++){
				rel_vars[i][0] = RANGE*pow(base_var[i], 1);
				rel_vars[i][1] = pow(base_var[i], 2); 
			}
		     }
		break;
		case 4:
		    switch(rtype){
			case COMPLEX:
			for(int i = 0;i<id->len;i++){
				rel_vars[i][0] = RANGE*pow(base_var[i], 1);
				rel_vars[i][1] = pow(base_var[i], 3);
				rel_vars[i][2] = pow(base_var[i], 2);
			}
		    }	
		break;
	}

        for(int i = 0;i<id->len;i++){
               for(int j = 0;j<id->npop;j++){
			if(j==0) id->data[i][j] = base_var[i];
			else id->data[i][j] = rel_vars[i][j];	
		}
	}

	return id;
}	

