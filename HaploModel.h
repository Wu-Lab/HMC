
#ifndef __HAPLOMODEL_H
#define __HAPLOMODEL_H


#include "Utils.h"
#include "HaploBuilder.h"


typedef enum { MC_v, MC_d, MC_b } HaploModelTypes;


class HaploModel : public HaploBuilder {
protected:
	HaploModelTypes m_model;
	int m_min_pattern_len;
	int m_max_pattern_len;
	int m_num_patterns;
	double m_min_freq_abs;
	double m_min_freq_rel;

public:
	HaploModel();

	int min_pattern_len() const { return m_min_pattern_len; }
	int max_pattern_len() const { return m_max_pattern_len; }
	int num_patterns() const { return m_num_patterns; }
	double min_freq_abs() const { return m_min_freq_abs; }
	double min_freq_rel() const { return m_min_freq_rel; }

	void setModel(HaploModelTypes model) { m_model = model; }
	void setMinPatternLen(int len) { m_min_pattern_len = len; }
	void setMaxPatternLen(int len) { m_max_pattern_len = len; }
	void setNumPatterns(int num) { m_num_patterns = num; }
	void setMinFreqAbs(double freq) { m_min_freq_abs = freq; }
	void setMinFreqRel(double freq) { m_min_freq_rel = freq; }

	double getMinFreq() const;

	void build(HaploData &hd);


protected:
	void findPatterns_mc_v();
	void findPatterns_mc_d();
	void findPatterns_mc_b();
};


#endif // __HAPLOMODEL_H