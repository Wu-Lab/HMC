
#include "HaploModel.h"
#include "GenoData.h"
#include "HaploComp.h"

#include <cfloat>

#include "MemLeak.h"


HaploModel::HaploModel()
{
	m_model = "MV";
	min_freq = -1;
	num_patterns = -1;
	min_pattern_len = 1;
	max_pattern_len = -1;
	mc_order = 1;
	max_iteration = -1;
	min_iteration = -1;
	sample_size = 1;
	max_sample_size = 1;
	final_sample_size = 1;
	exact_estimate = false;
}

void HaploModel::setModel(string model)
{
	if (model == "MV" || model == "MC" || model == "MA") {
		m_model = model;
	}
	else {
		Logger::error("Unknown model %s!", model.c_str());
		exit(1);
	}
}

void HaploModel::build(GenoData &genos)
{
	setGenoData(genos);

	Logger::info("");
	Logger::info("Running HMC engine %s ...", m_model.c_str());
	Logger::verbose("");
	Logger::beginTimer(1, "Search Haplotype pattern");

	findPatterns();

	Logger::endTimer(1);
}

void HaploModel::findPatterns()
{
	if (min_freq_abs > 0) {
		min_freq = min_freq_abs / (2.0 * genos()->genotype_num());
	}
	if (m_model == "MV") {
		if (num_patterns > 0) {
			m_patterns.findPatternByNum(num_patterns, min_pattern_len, max_pattern_len);
		}
		else {
			m_patterns.findPatternByFreq(min_freq, min_pattern_len, max_pattern_len);
		}
	}
	else if (m_model == "MC") {
		m_patterns.findPatternBlock(mc_order+1);
	}
	else if (m_model == "MA") {
		if (num_patterns > 0) {
			m_patterns.findPatternByNum(num_patterns, min_pattern_len, max_pattern_len);
		}
		else {
			m_patterns.findPatternByFreq(min_freq, min_pattern_len, max_pattern_len);
		}
		m_patterns.adjustFrequency();
	}
}

double HaploModel::resolveAll(GenoData &genos, GenoData &resolutions)
{
	int i, j, n;
	vector<Genotype> res_list;
	double sampling_coverage;
	double log_likelihood = 0;
	samples()->clear();
	for (i=0; i<genos.genotype_num(); ++i) {
		if (!genos[i].isPhased()) {
			Logger::status("  Resolving Genotype[%d] %s ...", i, genos[i].id().c_str());
			sampling_coverage = resolve(genos[i], resolutions[i], res_list, sample_size);
			resolutions[i].setID(genos[i].id());
			if (res_list.empty()) {
				Logger::warning("Unable to resolve Genotype[%d]: %s!", i, genos[i].id().c_str());
			}
			else {
				n = res_list.size();
				for (j=0; j<n; ++j) {
					res_list[j](0).setWeight(res_list[j].posterior_probability() / sampling_coverage);
					res_list[j](1).setWeight(res_list[j].posterior_probability() / sampling_coverage);
// 					for (int k=0; k<genos.genotype_len(); ++k) {
// 						if (genos[i].isMissing(k)) {
// 							res_list[j](0)[k] = -1;
// 							res_list[j](1)[k] = -1;
// 						}
// 					}
					samples()->addHaplotype(res_list[j](0));
					samples()->addHaplotype(res_list[j](1));
				}
			}
			genos[i].setGenotypeProbability(resolutions[i].genotype_probability());
			log_likelihood += log(resolutions[i].genotype_probability());
		}
	}
	samples()->checkTotalWeight();
	return log_likelihood;
}

void HaploModel::run(const GenoData &genos, GenoData &resolutions)
{
	int iter;
	double ll, old_ll;
	GenoData unphased, resolved;
	vector<Genotype> res_list;

	unphased = resolved = genos;
//	unphased.randomizePhase();
	build(unphased);
	resolutions = unphased;

	old_ll = -DBL_MAX;
	for (iter=1; iter<=max_iteration; ++iter) {
		ll = resolveAll(unphased, resolved);
		if (ll >= old_ll) resolutions = resolved;

		HaploComp compare(&genos, &resolutions);
		Logger::info("");
		Logger::info("  Switch Error = %f, IHP = %f, IGP = %f, LL = %f",
			compare.switch_error(), compare.incorrect_haplotype_percentage(), compare.incorrect_genotype_percentage(), ll);

		if (iter < max_iteration && ll >= old_ll && (old_ll - ll) / old_ll > 0.0001) {
 			if (exact_estimate) {
				m_patterns.estimatePatterns();
			}
			else {
				findPatterns();
			}
			old_ll = ll;
			if (m_model == "MA") {
				m_patterns.adjustFrequency();
			}
		}
		else {
			break;
		}
	}
}
