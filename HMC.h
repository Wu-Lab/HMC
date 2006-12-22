
#ifndef __HMC_H
#define __HMC_H


#include "HaploData.h"
#include "HaploFile.h"
#include "HaploModel.h"
#include "HaploComp.h"


class HMC {
	ArgParser m_args;
	HaploFile *m_haplo_file;
	HaploModel m_builder;
	HaploData m_genos;
	HaploData m_resolutions;

	static const char *m_version;
	static const char *m_year;

public:
	HMC(int argc = 0, char *argv[] = NULL);

	void usage();
	void copyright();
	void defineOptions();
	void parseOptions(int argc = 0, char *argv[] = NULL);

	void run();

	void resolve();
	void compareWith(const char *filename);

	const char *version() { return m_version; }
};


#endif // __HMC_H
