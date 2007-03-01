
#include "HaploModel.h"
#include "HaploData.h"
#include "HaploComp.h"

#include <cfloat>

#include "MemLeak.h"


HaploModel::HaploModel()
{
	m_model = MC_v;
	min_freq = -1;
	num_patterns = -1;
	min_pattern_len = 2;
	max_pattern_len = -1;
	mc_order = 1;
}

void HaploModel::build(HaploData &hd)
{
	setHaploData(hd);
	switch (m_model) {
	case MC_v:
		findPatterns_mc_v();
		break;
	case MC_d:
		findPatterns_mc_d();
		break;
	case MC_b:
		findPatterns_mc_b();
		break;
	}
	initialize();
}

void HaploModel::findPatterns_mc_v()
{
	Logger::info("");
	Logger::info("Running HMC engine MC-v ...");
	Logger::verbose("");
	Logger::beginTimer(1, "Search Haplotype pattern");
	if (num_patterns > 0) {
		m_patterns.findPatternByNum(num_patterns, min_pattern_len, max_pattern_len);
	}
	else {
		m_patterns.findPatternByFreq(min_freq, min_pattern_len, max_pattern_len);
	}
	Logger::endTimer(1);
}

void HaploModel::findPatterns_mc_d()
{
	Logger::info("");
	Logger::info("Run HMC engine MC-d ...");
	Logger::verbose("");
	Logger::beginTimer(1, "Search Haplotype pattern");
	m_patterns.findPatternBlock(mc_order+1);
	Logger::endTimer(1);
}

void HaploModel::findPatterns_mc_b()
{
	Logger::info("");
	Logger::info("Running HMC engine MC-b ...");
	Logger::verbose("");
	Logger::beginTimer(1, "Search Haplotype pattern");
	if (num_patterns > 0) {
		m_patterns.findPatternByNum(num_patterns, min_pattern_len, max_pattern_len);
	}
	else {
		m_patterns.findPatternByFreq(min_freq, min_pattern_len, max_pattern_len);
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
	resolutions = genos;

	unphased.randomizePhase();
	build(unphased);

	old_ll = -DBL_MAX;
	for (iter=1; iter<=max_iteration; ++iter) {
		ll = 0;
		for (i=0; i<genos.unphased_num(); ++i) {
//			if (!unphased[i].isPhased()) {
				Logger::status("Iteration %d: Resolving Genotype[%d] %s ...", iter, i, genos[i].id().c_str());
				resolve(genos[i], resolutions[i], res_list);
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
			if (ll >= old_ll && (old_ll - ll) / old_ll < 0.001) {
				m_patterns.adjustPatterns();
				old_ll = -DBL_MAX;
			}
			else {
				m_patterns.adjustFrequency();
				old_ll = ll;
			}
		}
	}
}
