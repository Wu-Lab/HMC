
#include "HMC.h"
#include "HaploFile.h"
#include "HaploComp.h"


const char *HMC::m_version = "0.7.2";
const char *HMC::m_year = "2007";


HMC::HMC(int argc, char *argv[])
{
	m_haplo_file = NULL;
	defineOptions();
	parseOptions(argc, argv);
}

HMC::~HMC()
{
	delete m_haplo_file;
}

void HMC::usage()
{
	char buf[10000];
	Logger::info("Usage: HMC [option ...] datafiles\n");
	Logger::info(m_args.printDef(buf));
}

void HMC::copyright()
{
	Logger::info("HMC - Haplotype inference tool based on Markov Chain [Version %s]", m_version);
	Logger::info("Copyright (C) ZHANGroup@BIOINFOAMSS.ORG 2004-%s. All rights reserved.\n", m_year);
}

void HMC::defineOptions()
{
	ArgOption *option;

	// debug settings
	option = m_args.addOption("debug_level", "d", true, "4", 5, "1", "2", "3", "4", "5");
	option->addHelpInfo("Set debug level, larger level for more debug information.");
	option = m_args.addOption("no_logo", "n");
	option->addHelpInfo("Suppress the logo and copyright information.");
	option = m_args.addOption("help", "h");
	option->addHelpInfo("Print the help message.");

	// comparison
	option = m_args.addOption("compare_with", "w", true);
	option->addHelpInfo("Compare the input data with the target data.");

	// convert data file format
	option = m_args.addOption("convert_to", "t", true);
	option->addHelpInfo("Convert data file format to default format.");
	option = m_args.addOption("randomize", "tr");
	option->addHelpInfo("Used with option -t, randomize the phase.");
	option = m_args.addOption("simplify", "ts");
	option->addHelpInfo("Used with option -t, simplify the allele symbols.");

	// input data file format
	option = m_args.addOption("file_format_0", "f0");	// PHASE like format (default)
	option->addHelpInfo("Read input file as PHASE like format (default).");
	option = m_args.addOption("file_format_1", "f1");	// HPM format
	option->addHelpInfo("Read input file as HaploRec like format.");
	option = m_args.addOption("file_format_2", "f2");	// Benchmark format (unrelated individuals)
	option->addHelpInfo("Read input file as Benchmark format (unrelated individuals).");
	option = m_args.addOption("file_format_3", "f3");	// Benchmark format (trios)
	option->addHelpInfo("Read input file as Benchmark format (trios).");

	// engine models
	option = m_args.addOption("model", "m", true, "MC-v", 3, "MC-v", "MC-d", "MC-b");
	option->addHelpInfo("Set the inference model.");

	// model options
	option = m_args.addOption("mc_order", "mo", true);
	option = m_args.addOption("min_freq_rel", "mfr", true);
	option = m_args.addOption("min_freq_abs", "mfa", true, "2.0");
	option = m_args.addOption("min_pattern_len", "minpl", true);
	option = m_args.addOption("max_pattern_len", "maxpl", true, "30");
	option = m_args.addOption("num_patterns", "np", true);

	option = m_args.addOption("iteration_number", "i", true, "1");
	option = m_args.addOption("output_patterns", "op", true);

	// for poor machine
	option = m_args.addOption("ram_limit", "ram", true, "10000000");
	option->addHelpInfo("Set memory limitation, used when the memory of the machine is poor.");

	// conflicts between options
	m_args.addConflictOptions("mc_order", "min_freq_rel");
	m_args.addConflictOptions("mc_order", "min_freq_abs");
	m_args.addConflictOptions("mc_order", "min_pattern_len");
	m_args.addConflictOptions("mc_order", "max_pattern_len");
	m_args.addConflictOptions("mc_order", "num_patterns");
	m_args.addConflictOptions("min_freq_rel", "min_freq_abs");
	m_args.addConflictOptions("min_freq_rel", "num_patterns");
	m_args.addConflictOptions("min_freq_abs", "num_patterns");
}

void HMC::parseOptions(int argc, char *argv[])
{
	char buf[10000];
	m_args.parseArguments(argc, argv);

	Logger::setLogLevel(m_args.getOption("debug_level")->getValueAsInt());
	if (!m_args.isDefined("no_logo")) {
		copyright();
	}
	if (m_args.isDefined("help") || m_args.getNonOptionArgNum() < 1) {
		usage();
		exit(0);
	}
	Constant::setRAMLimit(m_args.getOption("ram_limit")->getValueAsInt());

	//////////////////////////////////////////////////////////////////////////
	// input file format

	if (m_args.isDefined("file_format_1")) {	
		if (m_args.getNonOptionArgNum() != 1) {
			Logger::error("Option file_format_1 require 1 input files!");
			exit(1);
		}
		delete m_haplo_file;
		m_haplo_file = new HaploFileHPM(m_args.getNonOptionArgument(0));
	}
	else if (m_args.isDefined("file_format_2")) {
		if (m_args.getNonOptionArgNum() != 2) {
			Logger::error("Option file_format_2 require 2 input files!");
			exit(1);
		}
		delete m_haplo_file;
		m_haplo_file = new HaploFileBench(m_args.getNonOptionArgument(0), m_args.getNonOptionArgument(1));
	}
	else if (m_args.isDefined("file_format_3")) {
		if (m_args.getNonOptionArgNum() != 3) {
			Logger::error("Option file_format_3 require 3 input files!");
			exit(1);
		}
		delete m_haplo_file;
		m_haplo_file = new HaploFileBench(m_args.getNonOptionArgument(0), m_args.getNonOptionArgument(1), m_args.getNonOptionArgument(2));
	}
	else {
		if (m_args.getNonOptionArgNum() != 1) {
			Logger::error("Require 1 input files!");
			exit(1);
		}
		delete m_haplo_file;
		m_haplo_file = new HaploFile(m_args.getNonOptionArgument(0));
	}

	//////////////////////////////////////////////////////////////////////////
	// model parameters

	if (m_args.isDefined("min_pattern_len") && m_args.isDefined("max_pattern_len") && m_args.getOption("min_pattern_len")->getValueAsInt() >= m_args.getOption("max_pattern_len")->getValueAsInt()) {
		Logger::error("The value of option -max_pattern_len must greater than that of option -min_pattern_len!");
		exit(1);
	}
	if (m_args.getOption("model")->equal("MC-d") && !m_args.isDefined("mc_order")) {
		Logger::error("The option -mc_order is needed for MC-d model!");
		exit(1);
	}

	if (m_args.getOption("model")->equal("MC-v")) {
		m_builder.setModel(MC_v);
	}
	else if (m_args.getOption("model")->equal("MC-d")) {
		m_builder.setModel(MC_d);
	}
	else if (m_args.getOption("model")->equal("MC-b")) {
		m_builder.setModel(MC_b);
	}

	if (m_args.isDefined("min_pattern_len")) {
		m_builder.setMinPatternLen(m_args.getOption("min_pattern_len")->getValueAsInt());
	}
	else {
		m_builder.setMinPatternLen(1);
	}

	if (m_args.isDefined("max_pattern_len")) {
		m_builder.setMaxPatternLen(m_args.getOption("max_pattern_len")->getValueAsInt());
	}
	else {
		m_builder.setMaxPatternLen(-1);
	}

	if (m_args.isDefined("num_patterns")) {
		m_builder.setNumPatterns(m_args.getOption("num_patterns")->getValueAsInt());
	}
	else if (m_args.isDefined("min_freq_abs")) {
		m_builder.setMinFreqAbs(m_args.getOption("min_freq_abs")->getValueAsDouble());
		m_builder.setNumPatterns(-1);
	}
	else {
		m_builder.setMinFreqRel(m_args.getOption("min_freq_rel")->getValueAsDouble());
		m_builder.setMinFreqAbs(-1);
		m_builder.setNumPatterns(-1);
	}

	if (m_args.isDefined("mc_order")) {
		m_builder.setMinPatternLen(m_args.getOption("mc_order")->getValueAsInt());
		m_builder.setMaxPatternLen(m_args.getOption("mc_order")->getValueAsInt()+1);
	}

	Logger::debug(m_args.printVal(buf));
}

void HMC::run()
{
	// read input file
	Logger::info("Reading genotype file ...");
	m_haplo_file->readHaploData(m_genos);
	Logger::info("Succesfully read Haplotype file with %d markers and %d genotypes.",
					m_genos.genotype_len(), m_genos.genotype_num());

	if (m_args.isDefined("compare_with")) {
		compareWith(m_args.getOption("compare_with")->getValue());
	}
	else if (m_args.isDefined("convert_to")) {
		HaploFile newfile;
		if (m_args.isDefined("simplify")) {
			m_genos.simplify();
		}
		if (m_args.isDefined("randomize")) {
			m_genos.randomizePhase();
		}
		newfile.setFileName(m_args.getOption("convert_to")->getValue());
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

	max_iter = min(m_genos.unphased_num(), m_args.getOption("iteration_number")->getValueAsInt());
	for (iter=0; iter<max_iter; iter++) {
		double ll = 0;
		for (i=0; i<m_genos.unphased_num(); i++) {
//			if (!unphased_genos[i].isPhased()) {
				Logger::status("Iteration %d: Resolving Genotype[%d] %s ...", iter, i, m_genos[i].id());
				m_builder.resolve(unphased_genos[i], m_resolutions.genotype(i), res_list);
				m_resolutions.genotype(i).setID(m_genos[i].id());
				if (res_list.size() == 0) {
					Logger::warning("Unable to resolve Genotype[%d]: %s!", i, m_genos[i].id());
				}
				unphased_genos.genotype(i).setLikelihood(m_resolutions[i].likelihood());
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
	if (m_args.isDefined("output_patterns"))
	{
		m_haplo_file->writePattern(m_builder, ".patterns");
	}
}

void HMC::compareWith(const char *filename)
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
