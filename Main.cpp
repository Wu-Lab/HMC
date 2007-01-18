
#include "HMC.h"

#include "MemLeak.h"


int main(int argc, char *argv[])
{
	EnableMemLeakCheck();
//	_CrtSetBreakAlloc(1682);

	srand(time(NULL));

	try {
		HMC hmc(argc, argv);
		hmc.run();
	}
    catch (exception &e)
    {
        cout << e.what() << "\n";
        return 1;
	}
	return 0;
}
