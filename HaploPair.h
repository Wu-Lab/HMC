
#ifndef __HAPLOPAIR_H
#define __HAPLOPAIR_H


#include "Utils.h"
#include "Genotype.h"
#include "HaploPattern.h"


struct PatternPair {
	HaploPattern &pattern_a, &pattern_b;
	tr1::shared_ptr<PatternPair> prev;

	PatternPair(HaploPattern *hpa, HaploPattern *hpb) : pattern_a(*hpa), pattern_b(*hpb) {};
};


class HaploPair {
	tr1::shared_ptr<PatternPair> m_pair;
	int m_end, m_id[2];
	double m_likelihood;
	double m_total_likelihood;
	bool m_homogenous;
	bool m_half, m_match_a, m_match_b, m_match_next_a, m_match_next_b;

public:
	HaploPair(HaploPattern *hpa, HaploPattern *hpb, HaploPattern *target_pattern = NULL);

	int end() const { return m_end; }
	double likelihood() const { return m_likelihood; }

	Genotype getGenotype();

	bool extendable(int a, int b, HaploPattern *target_pattern = NULL);
	void extend_trial(int a, int b, int id[2], double &likelihood, double &total_likelihood);
	void extend(int a, int b);

	static void getID(int id[2], HaploPattern *hpa, HaploPattern *hpb);

	friend class HaploBuilder;
	friend bool greater_likelihood(const HaploPair *hp1, const HaploPair *hp2);
};


inline bool greater_likelihood(const HaploPair *hp1, const HaploPair *hp2)
{
	return hp1->m_likelihood > hp2->m_likelihood;
}


#endif // __HAPLOPAIR_H
