
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdexcept>

#include "HMC.h"
#include "HaploFile.h"
#include "HaploComp.h"
#include "Options.h"


const char *HMC::m_version = "0.8.0";
const char *HMC::m_year = "2007";


HMC::HMC(int argc, char *argv[])
{
	po::options_description generics("Generic options");
	generics.add_options()
		("help,h", "Show help message")
		("config,c", po::value<string>()->default_value("HMC.cfg"), "Configuration file")
		;

	po::options_description configs("Configurations");
	configs.add_options()
		("nologo", "Suppress logo and copyright information")
		("debug,d", po::value<int>()->default_value(4), "Set debug level")
		("input-format,f", po::value<string>(&m_input_format)->default_value("PHASE"), "Set input file format")
		("output-patterns", po::value<string>(), "")
		;

	po::options_description parameters("Model parameters");
	parameters.add_options()
		("model,m", po::value<string>()->default_value("MV"), "Set the inference model")
		("min-freq,r", po::value<double>(&m_builder.min_freq)->default_value(0.02), "Minimum relative frequency of patterns")
		("num-patterns,n", po::value<int>(&m_builder.num_patterns), "Maximum number of patterns")
		("min-pattern-len", po::value<int>(&m_builder.min_pattern_len)->default_value(1), "Minimum length of patterns")
		("max-pattern-len", po::value<int>(&m_builder.max_pattern_len)->default_value(30), "Maximum length of patterns")
		("mc-order,o", po::value<int>(&m_builder.mc_order)->default_value(1), "Markov chain order")
		("iteration,i", po::value<int>(&m_builder.max_iteration)->default_value(1), "")
		;

	po::options_description utilities("Utility options");
	utilities.add_options()
		("compare-with,w", po::value<string>(), "Compare input data with target data")
		("convert-to,t", po::value<string>(&m_convert_format), "Convert input data to specified format")
		("randomize", "Randomize genotype phases when converting format")
		("simplify", "Simplify allele symbols when converting format")
		;

	po::options_description hidden;
	hidden.add_options()
		("filename", po::value<vector<string> >(&m_filenames), "input/output files")
		;

	m_options.add(generics).add(configs).add(parameters).add(utilities).add(hidden);
	m_visible_options.add(generics).add(configs).add(parameters).add(utilities);

	po::positional_options_description p;
	p.add("filename", -1);

	po::options_description cmdline_options, file_options;
	cmdline_options.add(generics).add(configs).add(parameters).add(utilities).add(hidden);
	file_options.add(configs).add(parameters).add(utilities);

	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), m_args);

	ifstream ifs(m_args["config"].as<string>().c_str());
	po::store(po::parse_config_file(ifs, file_options), m_args);

	po::notify(m_args);

	conflicting_options(m_args, "min-freq", "num-patterns");

	parseOptions();
	parseFileNames();
}

HMC::~HMC()
{
}

void HMC::usage()
{
	Logger::info("Usage: HMC [option ...] datafiles");
	if (Logger::log_level() >= Logger::log_level_info()) {
		cout << m_visible_options << "\n";
	}
}

void HMC::copyright()
{
	Logger::info("HMC - Haplotype inference tool based on Markov Chain [Version %s]", m_version);
	Logger::info("Copyright (C) ZHANGroup@BIOINFOAMSS.ORG 2004-%s. All rights reserved.\n", m_year);
}

void HMC::printOptions()
{
	cout << "Used options:" << endl;
	print_options(m_args, cout, DisplayOption::defaulted);
	print_options(m_args, cout, DisplayOption::specified);
	cout << endl;
}

void HMC::parseOptions()
{
	Logger::setLogLevel(m_args["debug"].as<int>());
	if (!m_args.count("nologo")) {
		copyright();
	}
	if (m_args.count("help") || !m_filenames.size()) {
		usage();
		exit(0);
	}

	if (Logger::isDebug()) {
		printOptions();
	}

	//////////////////////////////////////////////////////////////////////////
	// model parameters

	if (m_args.count("min_pattern_len") && m_args.count("max_pattern_len")
		&& m_args["min_pattern_len"].as<int>() >= m_args["max_pattern_len"].as<int>()) {
		Logger::error("The value of option -max_pattern_len must greater than that of option -min_pattern_len!");
		exit(1);
	}

	m_builder.setModel(m_args["model"].as<string>());
}

void HMC::parseFileNames()
{
	int ni = 0;
	int nc = 0;

	ni = HaploFile::getFileNameNum(m_input_format);
	if (ni == 0) {
		Logger::error("Unknown input format %s!", m_input_format.c_str());
		exit(1);
	}
	if (m_args.count("convert-to")) {
		nc = HaploFile::getFileNameNum(m_convert_format);
		if (nc == 0) {
			Logger::error("Unknown convert format %s!", m_convert_format.c_str());
			exit(1);
		}
	}
	if (m_filenames.size() != (ni + nc)) {
		Logger::error("Input format %s require %d filenames!", m_input_format.c_str(), ni);
		if (m_args.count("convert-to")) {
			Logger::error("Convert format %s require %d filenames!", m_convert_format.c_str(), nc);
		}
		exit(1);
	}

	m_input_file.reset(HaploFile::getHaploFile(m_input_format, m_filenames.begin()));
	if (m_args.count("convert-to")) {
		m_convert_file.reset(HaploFile::getHaploFile(m_convert_format, m_filenames.begin()+ni));
	}
}

void HMC::run()
{
	// read input file
	Logger::info("Reading genotype file ...");
	m_input_file->readHaploData(m_genos);
	Logger::info("Succesfully read Haplotype file with %d markers and %d genotypes.",
					m_genos.genotype_len(), m_genos.genotype_num());

	if (m_args.count("compare-with")) {
		compareWith(m_args["compare-with"].as<string>());
	}
	else if (m_args.count("convert-to")) {
		convertFormat();
	}
	else {
		resolve();
	}
}

void HMC::resolve()
{
	Logger::debug("");
	Logger::beginTimer(3, "Find bottleneck");
	Logger::beginTimer(4, "Find bottleneck");

	Logger::pauseTimer(3);
	Logger::pauseTimer(4);

	Logger::debug("");
	Logger::beginTimer(2, "Resolve Genotype");

	m_builder.run(m_genos, m_resolutions);

	Logger::debug("");
	Logger::endTimer(4);
	Logger::debug("");
	Logger::endTimer(3);
	Logger::debug("");
	Logger::endTimer(2);
	Logger::info("Solving Time = %f", Logger::timer(1).time()+Logger::timer(2).time());

// 	for (i=0; i<m_genos.unphased_num(); i++) {
// 		double w1, w2, w3;
// 		w1 = m_builder.getLikelihood(m_genos[i]);
// 		w2 = m_builder.getLikelihood(m_resolutions[i]);
// 		w3 = m_resolutions[i].weight();
// 		Logger::debug("Genotype[%d] %g, %g, %g", i, w1, w2, w3);
// 	}

	m_input_file->writeHaploData(m_resolutions, ".reconstructed");
	if (m_args.count("output-patterns"))
	{
		m_input_file->writePattern(m_builder, ".patterns");
	}
}

void HMC::convertFormat()
{
	if (m_args.count("simplify")) {
		m_genos.simplify();
	}
	if (m_args.count("randomize")) {
		m_genos.randomizePhase();
	}
	m_convert_file->writeHaploData(m_genos);
}

void HMC::compareWith(const string &filename)
{
	HaploFile comp_file;
	HaploData comp_genos;

	comp_file.setFileName(filename);
	comp_file.readHaploData(comp_genos);
	Logger::info("Succesfully read Haplotype file with %d markers and %d genotypes.",
					comp_genos.genotype_len(), comp_genos.genotype_num());
	HaploComp compare(&comp_genos, &m_genos);
	Logger::info("");
	Logger::info("Switch Error = %f, IHP = %f, IGP = %f",
		compare.switch_error(), compare.incorrect_haplotype_percentage(), compare.incorrect_genotype_percentage());
}
