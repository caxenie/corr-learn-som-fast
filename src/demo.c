#include "demo.h"

#define MAX_EPOCHS	1000
#define N_POP		2
#define POP_SIZE	200

int main(int argc, char* argv[])
{
	network *n = init_network(N_POP, POP_SIZE); 
	simulation *s = init_simulation(MAX_EPOCHS, n);
	
	deinit_simulation(s);
	return EXIT_SUCCESS;
}
