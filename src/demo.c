#include "simulation.h"

int main(int argc, char* argv[])
{
	indata *indata = generate_input_data(N_POP, POP_SIZE, DATASET_LEN, ORDER2);
	network *net = init_network(indata->npop, indata->popsize); 
	simulation *sim = init_simulation(MAX_EPOCHS, net);
	outdata *runtime = run_simulation(indata, sim);
	char* dump_file = dump_runtime_data_extended(runtime);
	printf("CORR_LEARN_NET: Runtime data dumped on disk in: %s\n", dumnp_file);
	deinit_simulation(sim);
	return EXIT_SUCCESS;
}
