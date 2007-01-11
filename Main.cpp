
#include "HMC.h"

#include "MemLeak.h"


int main(int argc, char *argv[])
{
	EnableMemLeakCheck();
//	_CrtSetBreakAlloc(1682);

	srand(time(NULL));

	HMC hmc(argc, argv);
	hmc.run();
	return 0;
}
