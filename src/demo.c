#include "simulation.h"

#define MAX_EPOCHS	1000
#define N_POP		2
#define POP_SIZE	100
#define DATASET_LEN	1500

int main(int argc, char* argv[])
{
	indata *indata = generate_input_data(N_POP, POP_SIZE, DATASET_LEN);
	network *n = init_network(indata->npop, indata->popsize); 
	simulation *s = init_simulation(MAX_EPOCHS, n);

	
	deinit_simulation(s);
	return EXIT_SUCCESS;
}
