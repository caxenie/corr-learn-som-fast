#include "tools.h"

/* floating point radom number generator between 0 and 1 */
double randf()
{
	return (double) rand () / (double) RAND_MAX;
}

/* Normal random variate vector generator, of mean m and std. dev. s */
double* gauss_rand(int num_vals, double mean, double std_dev)
{
        double x1 = 0.0f, x2 = 0.0f, w = 0.0f, y1 = 0.0f;
	double *out = (double*)calloc(num_vals, sizeof(double));
	for(int vidx = 0; vidx < num_vals; vidx++){
	        do
        	{
                	x1 = 2.0 * (double)rand()/RAND_MAX - 1.0;
	                x2 = 2.0 * (double)rand()/RAND_MAX - 1.0;
        	        w = x1 * x1 + x2 * x2;
	        } while ( w >= 1.0 );
	
        	w = sqrt( (-2.0 * log( w ) ) / w );
	        y1 = x1 * w;
		out[vidx] = y1;
	}
        return out;
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
                                       out = gauss_rand(numv, 0.0, vmax/(numv/2));
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

