
#ifndef __HMC_H
#define __HMC_H


#include <iostream>
#include <boost/program_options.hpp>

#include "Utils.h"
#include "HaploData.h"
#include "HaploModel.h"


namespace po = boost::program_options;


class HaploFile;


class HMC {
	po::options_description m_visible_options;
	po::variables_map m_args;
	string m_input_format;
	vector<string> m_input_files;
	HaploFile *m_haplo_file;
	HaploModel m_builder;
	HaploData m_genos;
	HaploData m_resolutions;

	static const char *m_version;
	static const char *m_year;

public:
	explicit HMC(int argc = 0, char *argv[] = NULL);
	~HMC();

	void usage();
	void copyright();
	void defineOptions();
	void parseOptions(int argc = 0, char *argv[] = NULL);

	void run();

	void resolve();
	void compareWith(const string &filename);

	const char *version() { return m_version; }
};


#endif // __HMC_H
