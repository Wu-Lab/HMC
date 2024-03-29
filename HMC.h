
#ifndef __HMC_H
#define __HMC_H


#include <boost/program_options.hpp>

#include "Utils.h"
#include "GenoData.h"
#include "HaploModel.h"


namespace po = ::boost::program_options;


class HaploFile;


class HMC {
	po::options_description m_options;
	po::options_description m_visible_options;
	po::variables_map m_args;
	vector<string> m_filenames;
	string m_input_format, m_convert_format;
	tr1::shared_ptr<HaploFile> m_input_file, m_target_file;

	HaploModel m_builder;
	GenoData m_genos;
	GenoData m_resolutions;

	static const char *m_version;
	static const char *m_year;

public:
	explicit HMC(int argc = 0, char *argv[] = NULL);
	~HMC();

	void usage();
	void copyright();

	void printOptions();
	void parseOptions();
	void parseFileNames();

	void run();

	void resolve();

	void convert();
	void compare();
};


#endif // __HMC_H
