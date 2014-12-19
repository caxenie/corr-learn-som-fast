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
/* sort generated data vector */
void sort_input_data(int n, double *ra)
{
        int l, j, ir, i;
        double rra;

        l = (n >> 1) + 1;
        ir = n;

        /* The index l will be decremented from its initial value down to 1 during
         * the "hiring" (heap creation) phase.  Once it reaches 1, the index ir will
         * be decremented from its initial value down to 1 during the "retirement-
         * and-promotion" (heap selection) phase.
        */
        for ( ; ; )
        {
                if (l > 1)                                      /* still in hiring phase */
                        rra = ra[--l];
                else                                            /* in retirement-and-promotion phase */
                {
                        rra = ra[ir];           /* clear a space at end of array */
                        ra[ir]=ra[1];                   /* retire the top of the heap into it */
                        if (--ir == 1)                  /* done with last promotion */
                        {
                                ra[1] = rra;
	                        return;
                        }                                               /* end if */
                }                                                       /* end else */
                i = l;                                          /* whether we are in the hiring phase */
                j = l << 1;                                     /* or promotion phase, we here set up */
                while ( j <= ir )
                {
                        if ( (j < ir) && (ra[j] < ra[j + 1]) )
                                ++j;                            /* compare to the better underling */
                                if ( rra < ra[j] )      /* demote rra */
                                {
                                        ra[i] = ra[j];
                                        j += (i = j);
                                }
                                else
                                        j = ir + 1;             /* this is rra's level; set j to */
                }                           /* terminate the sift-down */
                ra[i] = rra;                            /* put rra into its slot */
        }
}

/* swap data points */
void swap_data_points(double *in1, double *in2)
{
	double tmp = *in1; *in1 = *in2; *in2 = tmp;	
}

/* shuffle generated data vector */
void shuffle_input_data(double *arr, int n)
{
	/* Use a different seed value so that we don't get same result each time we run this program */
	srand ( time(NULL) );
	for (int i = n-1; i > 0; i--){
		int j = rand() % (i+1);
		swap_data_points(&arr[i], &arr[j]);	
	}
}

/* non-uniform random distribution generator */
double * generate_rnd_vector(int type, int range, int numv, int dist, int dtype)
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
			/* random data */
			// for(int i=0;i<numv; i++)
			//	out[i] = -range + randf()*(2*range);
			/* sequential data */
			out[0] = -range;
			for(int i=1; i<numv; i++){
				out[i] = out[i-1] + (double)(2*range)/numv; // non-random values
			}
			if(dtype==TRAINING)
				/* randomize the data */
			        shuffle_input_data(out, numv);	
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

