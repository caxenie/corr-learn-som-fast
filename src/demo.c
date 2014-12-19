#include "simulation.h"

int main(int argc, char* argv[])
{
	indata *indata = generate_input_data(N_POP, POP_SIZE, DATASET_LEN, ORDER2, TRAINING);
	network *net = init_network(indata->npop, indata->popsize); 
	simulation *sim = init_simulation(MAX_EPOCHS, net);
	outdata *runtime = run_simulation(indata, sim);
	char* dump_file = dump_runtime_data_extended(runtime, BASIC);
	printf("CORR_LEARN_NET: Runtime data dumped on disk in: %s\n", dump_file);
#ifdef TESTS_ON
	runtime->in = generate_input_data(N_POP, POP_SIZE, DATASET_LEN, ORDER2, TESTING);
	outdata *test_runtime = test_inference(runtime);
        char* dump_file_test_runtime = dump_runtime_data_extended(test_runtime, TESTS);
        printf("CORR_LEARN_NET_TESTS: Inference test runtime data dumped on disk in: %s\n", dump_file_test_runtime);
#endif 
	deinit_simulation(sim);
	return EXIT_SUCCESS;
}
