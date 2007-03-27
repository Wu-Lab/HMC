
#include "HaploModel.h"
#include "HaploData.h"
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

void HaploModel::build(HaploData &hd)
{
	setHaploData(hd);
	findPatterns();
	initialize();
}

void HaploModel::findPatterns()
{
	Logger::info("");
	Logger::info("Running HMC engine %s ...", m_model.c_str());
	Logger::verbose("");
	Logger::beginTimer(1, "Search Haplotype pattern");

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

	Logger::endTimer(1);
}

void HaploModel::run(const HaploData &genos, HaploData &resolutions)
{
	int i, iter;
	double ll, old_ll;
	HaploData unphased;
	vector<HaploPair*> res_list;

	unphased = genos;
	unphased.randomizePhase();
	build(unphased);
	resolutions = unphased;

	old_ll = -DBL_MAX;
	for (iter=1; iter<=max_iteration; ++iter) {
		ll = 0;
		for (i=0; i<genos.unphased_num(); ++i) {
//			if (!unphased[i].isPhased()) {
				Logger::status("Iteration %d: Resolving Genotype[%d] %s ...", iter, i, genos[i].id().c_str());
				resolve(unphased[i], resolutions[i], res_list);
				resolutions[i].setID(genos[i].id());
				if (res_list.size() == 0) {
					Logger::warning("Unable to resolve Genotype[%d]: %s!", i, genos[i].id().c_str());
				}
				unphased[i].setLikelihood(resolutions[i].likelihood());
//			}
			ll += log(resolutions[i].likelihood());
		}

		HaploComp compare(&genos, &resolutions);
		Logger::info("");
		Logger::info("  Switch Error = %f, IHP = %f, IGP = %f, LL = %f",
			compare.switch_error(), compare.incorrect_haplotype_percentage(), compare.incorrect_genotype_percentage(), ll);

		if (iter < max_iteration) {
// 			if (ll >= old_ll && (old_ll - ll) / old_ll < 0.01) {
 				m_patterns.estimatePatterns();
// 				old_ll = -DBL_MAX;
// 			}
// 			else {
//				m_patterns.estimateFrequency();
// 				old_ll = ll;
// 			}
			if (m_model == "MA") {
				m_patterns.adjustFrequency();
			}
		}
	}
}
