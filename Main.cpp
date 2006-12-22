
#include "HMC.h"


int main(int argc, char *argv[])
{
	srand(time(NULL));

	HMC hmc(argc, argv);
	hmc.run();
	return 0;
}
