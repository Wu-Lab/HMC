
#ifndef __HAPLOPAIR_H
#define __HAPLOPAIR_H


#include "Genotype.h"
#include "HaploPattern.h"


class HaploPair {
	int m_end;
	int m_genotype_len;
	HaploPattern **m_patterns[2];
	double m_weight;
	double m_likelihood;
	double m_total_likelihood;
	bool m_homogenous;

public:
	HaploPair(HaploPattern *hpa, HaploPattern *hpb);
	HaploPair(const HaploPair &hp);
	~HaploPair();

	int end() const { return m_end; }
	HaploPattern *patterns(int i) const { return m_patterns[i][m_end]; }
	double likelihood() const { return m_likelihood; }

	Genotype getGenotype();

	void getID(int id[2]);
	bool extendable(int a, int b);
	void extend(int a, int b, HaploPattern *hpa = NULL, HaploPattern *hpb = NULL);

	friend class HaploBuilder;
};


#endif // __HAPLOPAIR_H
