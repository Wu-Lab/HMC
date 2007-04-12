
#ifndef __HAPLOMODEL_H
#define __HAPLOMODEL_H


#include "Utils.h"
#include "HaploBuilder.h"


class HaploModel : public HaploBuilder {
protected:
	string m_model;

public:
	double min_freq;
	int num_patterns;
	int min_pattern_len;
	int max_pattern_len;
	int mc_order;
	int max_iteration;
	int min_iteration;
	int sample_size;
	int max_sample_size;
	int final_sample_size;
	bool re_estimate_without_sampling;

public:
	HaploModel();

	void setModel(string model);

	void run(const GenoData &genos, GenoData &resolutions);

protected:
	void build(GenoData &genos);
	void findPatterns();

	double resolveAll(GenoData &genos, GenoData &resolutions);
};


#endif // __HAPLOMODEL_H