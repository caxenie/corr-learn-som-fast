#include "simulation.h"

#define MAX_EPOCHS	1000
#define N_POP		2
#define POP_SIZE	100
#define DATASET_LEN	1500

int main(int argc, char* argv[])
{
	indata *indata = generate_input_data(N_POP, POP_SIZE, DATASET_LEN);
	network *net = init_network(indata->npop, indata->popsize); 
	simulation *sim = init_simulation(MAX_EPOCHS, net);
	outdata *rutime = run_simulation(indata, sim);
	deinit_simulation(sim);
	return EXIT_SUCCESS;
}
