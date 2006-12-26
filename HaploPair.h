
#ifndef __HAPLOPAIR_H
#define __HAPLOPAIR_H


#include "Genotype.h"
#include "HaploPattern.h"


class HaploPair {
	int m_end;
	int m_genotype_len;
	HaploPattern **m_patterns[2];
	double m_likelihood;
	double m_total_likelihood;
	bool m_homogenous;
	bool m_half, m_match_a, m_match_b, m_match_next_a, m_match_next_b;

public:
	HaploPair(HaploPattern *hpa, HaploPattern *hpb, HaploPattern *target_pattern = NULL);
	HaploPair(const HaploPair &hp);
	~HaploPair();

	int end() const { return m_end; }
	HaploPattern *patterns(int i) const { return m_patterns[i][m_end]; }
	double likelihood() const { return m_likelihood; }

	Genotype getGenotype();

	void getID(int id[2]);
	bool extendable(int a, int b, HaploPattern *target_pattern = NULL);
	void extend(int a, int b);

	friend class HaploBuilder;
};


#endif // __HAPLOPAIR_H
