#include "tools.h"

/* floating point radom number generator between 0 and 1 */
double randf()
{
	return (double) rand () / (double) RAND_MAX;
}

/* non-uniform random distribution generator */
double * generate_rnd_vector(int type, int range, int numv, int dist)
{	
	double *out = (double*)calloc(numv, sizeof(double));
	double *var = (double*)calloc(numv, sizeof(double));
	int expn = 3;
	int vmin = 0.0;
	int vmax = range;
	
	for(int i=0;i<numv; i++){
		var[i] = randf()*range;
	}
	
	switch(type){
		case UNIFORM:
			for(int i=0;i<numv; i++){
				out[i] = -range + randf()*(2*range);
			}
		break;
		case NONUNIFORM:
			switch(dist){
				case INCPOWERLAW:
					for(int i=0;i<numv; i++){
						out[i] = exp(log(var[i]*(pow(-vmin,(expn+1)) + pow(vmax, (expn+1))))/(expn+1));
					}	
				break;
				case DECPOWERLAW:
					for(int i=0;i<numv; i++){
						out[i] = -exp(log(var[i]*(pow(-vmin,(expn+1)) + pow(vmax,(expn+1))))/(expn+1));
					}	
				break;
				case GAUSS:
					for(int i=0;i<numv; i++){
                                                out[i] = randf()*(vmax/4);
                                        }
				break;
				case CONVEX:
					for(int li = 0;li<numv;li++){
						if(li<=numv/2){
							vmin = 0.00000005f;
							out[li] = -exp(log(var[li]*(pow(-vmin, (expn+1)) + pow(vmax, (expn+1))))/(expn+1));
						}		
						else{
							vmin = -vmin;
							out[li] = exp(log(var[li]*(pow(-vmin, (expn+1)) + pow(vmax, (expn+1))))/(expn+1));
						}
					}
				break;
			}
	}
	return out;	
}

