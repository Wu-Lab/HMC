
#include <iostream>
#include <fstream>
#include <iterator>
#include <stdexcept>

#include "HMC.h"
#include "HaploFile.h"
#include "HaploComp.h"
#include "Options.h"


const char *HMC::m_version = "0.7.4";
const char *HMC::m_year = "2007";


HMC::HMC(int argc, char *argv[])
: m_haplo_file(0)
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
		("model,m", po::value<string>()->default_value("MC-v"), "Set the inference model")
		("mc-order,o", po::value<int>(), "Markov chain order")
		("min-freq-rel,r", po::value<double>(), "Minimum relative frequency of patterns")
		("min-freq-abs,a", po::value<double>()->default_value(2.0), "Minimum absolute frequency of patterns")
		("min-pattern-len", po::value<int>(), "Minimum length of patterns")
		("max-pattern-len", po::value<int>()->default_value(30), "Maximum length of patterns")
		("num-patterns", po::value<int>(), "Maximum number of patterns")
		("iteration,i", po::value<int>()->default_value(1), "")
		;

	po::options_description utilities("Utility options");
	utilities.add_options()
		("compare-with,w", po::value<string>(), "Compare input data with target data")
		("convert-to,t", po::value<string>(), "Convert input data to default format")
		("randomize", "Randomize genotype phases when convert")
		("simplify", "Simplify allele symbols when convert")
		;

	po::options_description hidden;
	hidden.add_options()
		("input-file", po::value<vector<string> >(&m_input_files), "input file")
		;

	m_options.add(generics).add(configs).add(parameters).add(utilities).add(hidden);
	m_visible_options.add(generics).add(configs).add(parameters).add(utilities);

	po::positional_options_description p;
	p.add("input-file", -1);

	po::options_description cmdline_options, file_options;
	cmdline_options.add(generics).add(configs).add(parameters).add(utilities).add(hidden);
	file_options.add(configs).add(parameters).add(utilities).add(hidden);

	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), m_args);

	ifstream ifs(m_args["config"].as<string>().c_str());
	po::store(po::parse_config_file(ifs, file_options), m_args);

	po::notify(m_args);

	conflicting_options(m_args, "min-freq-rel", "min-freq-abs");
	conflicting_options(m_args, "min-freq-rel", "num-patterns");
	conflicting_options(m_args, "min-freq-abs", "num-patterns");

	parseOptions();
}

HMC::~HMC()
{
	delete m_haplo_file;
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
	if (m_args.count("help") || !m_input_files.size()) {
		usage();
		exit(0);
	}

	if (Logger::isDebug()) {
		printOptions();
	}

	//////////////////////////////////////////////////////////////////////////
	// input file format

	if (m_input_format == "PHASE") {
		if (m_input_files.size() != 1) {
			Logger::error("Input format PHASE only require 1 input file!");
			exit(1);
		}
		delete m_haplo_file;
		m_haplo_file = new HaploFile(m_input_files[0]);
	}
	else if (m_input_format == "HPM") {
		if (m_input_files.size() != 1) {
			Logger::error("Input format HPM only require 1 input file!");
			exit(1);
		}
		delete m_haplo_file;
		m_haplo_file = new HaploFileHPM(m_input_files[0]);
	}
	else if (m_input_format == "BENCH2") {
		if (m_input_files.size() != 2) {
			Logger::error("Input format BENCH2 require 2 input files!");
			exit(1);
		}
		delete m_haplo_file;
		m_haplo_file = new HaploFileBench(m_input_files[0], m_input_files[1]);
	}
	else if (m_input_format == "BENCH3") {
		if (m_input_files.size() != 3) {
			Logger::error("Input format BENCH3 require 3 input files!");
			exit(1);
		}
		delete m_haplo_file;
		m_haplo_file = new HaploFileBench(m_input_files[0], m_input_files[1], m_input_files[2]);
	}
	else {
		Logger::error("Unknown input format %s!", m_input_format.c_str());
		exit(1);
	}

	//////////////////////////////////////////////////////////////////////////
	// model parameters

	if (m_args.count("min_pattern_len") && m_args.count("max_pattern_len")
		&& m_args["min_pattern_len"].as<int>() >= m_args["max_pattern_len"].as<int>()) {
		Logger::error("The value of option -max_pattern_len must greater than that of option -min_pattern_len!");
		exit(1);
	}

	if (m_args["model"].as<string>() == "MC-v") {
		m_builder.setModel(MC_v);
	}
	else if (m_args["model"].as<string>() == "MC-d") {
		m_builder.setModel(MC_d);
	}
	else if (m_args["model"].as<string>() == "MC-b") {
		m_builder.setModel(MC_b);
	}
	else {
		Logger::error("Unknown model %s!", m_args["model"].as<string>().c_str());
		exit(1);
	}

	if (m_args.count("min-pattern-len")) {
		m_builder.setMinPatternLen(m_args["min-pattern-len"].as<int>());
	}
	else {
		m_builder.setMinPatternLen(1);
	}

	if (m_args.count("max-pattern-len")) {
		m_builder.setMaxPatternLen(m_args["max-pattern-len"].as<int>());
	}
	else {
		m_builder.setMaxPatternLen(-1);
	}

	if (m_args.count("num-patterns")) {
		m_builder.setNumPatterns(m_args["num-patterns"].as<int>());
	}
	else if (m_args.count("min-freq-abs")) {
		m_builder.setMinFreqAbs(m_args["min-freq-abs"].as<double>());
		m_builder.setNumPatterns(-1);
	}
	else if (m_args.count("min-freq-rel")) {
		m_builder.setMinFreqRel(m_args["min-freq-rel"].as<double>());
		m_builder.setMinFreqAbs(-1);
		m_builder.setNumPatterns(-1);
	}

	if (m_args.count("mc-order")) {
		m_builder.setMinPatternLen(m_args["mc-order"].as<int>());
		m_builder.setMaxPatternLen(m_args["mc-order"].as<int>()+1);
	}
}

void HMC::run()
{
	// read input file
	Logger::info("Reading genotype file ...");
	m_haplo_file->readHaploData(m_genos);
	Logger::info("Succesfully read Haplotype file with %d markers and %d genotypes.",
					m_genos.genotype_len(), m_genos.genotype_num());

	if (m_args.count("compare-with")) {
		compareWith(m_args["compare-with"].as<string>());
	}
	else if (m_args.count("convert-to")) {
		HaploFile newfile;
		if (m_args.count("simplify")) {
			m_genos.simplify();
		}
		if (m_args.count("randomize")) {
			m_genos.randomizePhase();
		}
		newfile.setFileName(m_args["convert-to"].as<string>());
		newfile.writeHaploData(m_genos);
	}
	else {
		resolve();
	}
}

void HMC::resolve()
{
	int i, iter, max_iter;
	HaploData unphased_genos;
	vector<HaploPair*> res_list;

	Logger::verbose("");
	Logger::beginTimer(3, "Find bottleneck");
	Logger::beginTimer(4, "Find bottleneck");

	Logger::pauseTimer(3);
	Logger::pauseTimer(4);

	unphased_genos = m_genos;
	unphased_genos.randomizePhase();
	HaploComp compare(&m_genos, &unphased_genos);
	Logger::info("");
	Logger::info("  Switch Error = %f, IHP = %f, IGP = %f",
		compare.switch_error(), compare.incorrect_haplotype_percentage(), compare.incorrect_genotype_percentage());
	m_resolutions = unphased_genos;
	m_builder.build(unphased_genos);

	Logger::verbose("");
	Logger::beginTimer(2, "Resolve Genotype");

	max_iter = min(m_genos.unphased_num(), m_args["iteration"].as<int>());
	for (iter=0; iter<max_iter; iter++) {
		double ll = 0;
		for (i=0; i<m_genos.unphased_num(); i++) {
//			if (!unphased_genos[i].isPhased()) {
				Logger::status("Iteration %d: Resolving Genotype[%d] %s ...", iter, i, m_genos[i].id().c_str());
				m_builder.resolve(unphased_genos[i], m_resolutions[i], res_list);
				m_resolutions[i].setID(m_genos[i].id());
				if (res_list.size() == 0) {
					Logger::warning("Unable to resolve Genotype[%d]: %s!", i, m_genos[i].id().c_str());
				}
				unphased_genos[i].setLikelihood(m_resolutions[i].likelihood());
//			}
			ll += log(m_resolutions[i].likelihood());
		}

		HaploComp compare(&m_genos, &m_resolutions);
		Logger::info("");
		Logger::info("  Switch Error = %f, IHP = %f, IGP = %f, LL = %f",
			compare.switch_error(), compare.incorrect_haplotype_percentage(), compare.incorrect_genotype_percentage(), ll);

		if (iter < max_iter - 1) m_builder.adjust(2);
	}
	Logger::verbose("");
	Logger::endTimer(4);
	Logger::verbose("");
	Logger::endTimer(3);
	Logger::verbose("");
	Logger::endTimer(2);
	Logger::info("Solving Time = %f", Logger::timer(1).time()+Logger::timer(2).time());

// 	for (i=0; i<m_genos.unphased_num(); i++) {
// 		double w1, w2, w3;
// 		w1 = m_builder.getLikelihood(m_genos[i]);
// 		w2 = m_builder.getLikelihood(m_resolutions[i]);
// 		w3 = m_resolutions[i].weight();
// 		Logger::debug("Genotype[%d] %g, %g, %g", i, w1, w2, w3);
// 	}

	m_haplo_file->writeHaploData(m_resolutions, ".reconstructed");
	if (m_args.count("output-patterns"))
	{
		m_haplo_file->writePattern(m_builder, ".patterns");
	}
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
