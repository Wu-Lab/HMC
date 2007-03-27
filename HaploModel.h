
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

public:
	HaploModel();

	void setModel(string model);

	void build(HaploData &hd);

	void run(const HaploData &genos, HaploData &resolutions);

protected:
	void findPatterns();
};


#endif // __HAPLOMODEL_H