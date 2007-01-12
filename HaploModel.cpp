
#include "HaploModel.h"
#include "HaploData.h"

#include "MemLeak.h"


HaploModel::HaploModel()
{
	m_model = MC_v;
	m_min_pattern_len = 1;
	m_max_pattern_len = 2;
	m_num_patterns = -1;
	m_min_freq_abs = -1;
	m_min_freq_rel = -1;
}

double HaploModel::getMinFreq() const
{
	if (m_min_freq_abs > 0) {
		return m_min_freq_abs;
	}
	else {
		return m_min_freq_rel * m_haplo_data->genotype_num();
	}
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
	if (m_num_patterns > 0) {
		findHaploPatternByNum(m_num_patterns, m_min_pattern_len, m_max_pattern_len);
	}
	else if (m_min_freq_abs > 0) {
		findHaploPatternByFreq(m_min_freq_abs, m_min_pattern_len, m_max_pattern_len);
	}
	else {
		findHaploPatternByFreq(m_min_freq_rel * m_haplo_data->genotype_num(), m_min_pattern_len, m_max_pattern_len);
	}
	Logger::endTimer(1);
}

void HaploModel::findPatterns_mc_d()
{
	Logger::info("");
	Logger::info("Run HMC engine MC-d ...");
	Logger::verbose("");
	Logger::beginTimer(1, "Search Haplotype pattern");
	findHaploPatternByFreq(0, m_min_pattern_len, m_max_pattern_len);
	Logger::endTimer(1);
}

void HaploModel::findPatterns_mc_b()
{
	Logger::info("");
	Logger::info("Running HMC engine MC-b ...");
	Logger::verbose("");
	Logger::beginTimer(1, "Search Haplotype pattern");
	if (m_num_patterns > 0) {
		findHaploPatternBlockByNum(m_num_patterns, m_min_pattern_len, m_max_pattern_len);
	}
	else if (m_min_freq_abs > 0) {
		findHaploPatternBlockByFreq(m_min_freq_abs, m_min_pattern_len, m_max_pattern_len);
	}
	else {
		findHaploPatternBlockByFreq(m_min_freq_rel * m_haplo_data->genotype_num(), m_min_pattern_len, m_max_pattern_len);
	}
	Logger::endTimer(1);
}
