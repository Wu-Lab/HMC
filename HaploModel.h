
#ifndef __HAPLOMODEL_H
#define __HAPLOMODEL_H


#include "Utils.h"
#include "HaploBuilder.h"


typedef enum { MC_v, MC_d, MC_b } HaploModelTypes;


class HaploModel : public HaploBuilder {
protected:
	HaploModelTypes m_model;

public:
	double min_freq;
	int num_patterns;
	int min_pattern_len;
	int max_pattern_len;
	int mc_order;
	int max_iteration;

public:
	HaploModel();

	void setModel(HaploModelTypes model) { m_model = model; }

	void build(HaploData &hd);

	void run(const HaploData &genos, HaploData &resolutions);

protected:
	void findPatterns_mc_v();
	void findPatterns_mc_d();
	void findPatterns_mc_b();
};


#endif // __HAPLOMODEL_H