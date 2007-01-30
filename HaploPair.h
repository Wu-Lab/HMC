
#ifndef __HAPLOPAIR_H
#define __HAPLOPAIR_H


#include <boost/pool/pool.hpp>

#include "Utils.h"
#include "HaploPattern.h"


class HaploPair {
	static boost::pool<> m_pool;

	const HaploPattern &m_pattern_a, &m_pattern_b;
	const Allele &m_allele_a, &m_allele_b;
	vector<HaploPair*> m_forward_links;
	HaploPair *m_backward_link;
	double m_transition_prob;
	double m_forward_likelihood;
	double m_backward_likelihood;
	double m_best_likelihood;

public:
	explicit HaploPair(const HaploPattern *hpa, const HaploPattern *hpb);
	explicit HaploPair(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb);
	void add(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb);

	const HaploPattern &pattern_a() const { return m_pattern_a; }
	const HaploPattern &pattern_b() const { return m_pattern_b; }
	const Allele &allele_a() const { return m_allele_a; }
	const Allele &allele_b() const { return m_allele_b; }
	int id_a() const { return m_pattern_a.id(); }
	int id_b() const { return m_pattern_b.id(); }
	int end() const { return m_pattern_a.end(); }

	double transition_prob() const { return m_transition_prob; }
	double forward_likelihood() const { return m_forward_likelihood; }
	double backward_likelihood() const { return m_backward_likelihood; }
	double best_likelihood() const { return m_best_likelihood; }

	Genotype getGenotype();

	template <typename T>
	const HaploPattern *successor_a(T &i) const { return m_pattern_a.successors(i); }
	template <typename T>
	const HaploPattern *successor_b(T &i) const { return m_pattern_b.successors(i); }

	void calcBackwardLikelihood();
	static double evaluatePattern(const HaploPattern *hp, vector<HaploPair*> &hp_list);

	static void *operator new(std::size_t) { return m_pool.malloc(); }
	static void operator delete(void *rawMemory) { m_pool.free(rawMemory); }

#ifdef _DEBUG
	static void *operator new(unsigned int, int, const char *, int) { return m_pool.malloc(); }
	static void operator delete(void *rawMemory, int, const char *, int) { m_pool.free(rawMemory); }
#endif // _DEBUG

	struct greater_likelihood {
		bool operator()(const HaploPair *hp1, const HaploPair *hp2) const
		{
			return hp1->best_likelihood() > hp2->best_likelihood();
		}
	};

	friend class HaploBuilder;

private:
	HaploPair(const HaploPair &);
	HaploPair &operator=(const HaploPair &);
};


#endif // __HAPLOPAIR_H
