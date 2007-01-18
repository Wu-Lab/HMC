
#ifndef __HAPLOPAIR_H
#define __HAPLOPAIR_H


#include <boost/shared_ptr.hpp>
#include <boost/pool/pool.hpp>

#include "Utils.h"
#include "HaploPattern.h"


class PatternPair;
typedef tr1::shared_ptr<PatternPair> SP_PatternPair;

class PatternPair {
	static boost::pool<> m_pool;
	const HaploPattern &m_pattern_a, &m_pattern_b;
	SP_PatternPair m_prev;

public:
	explicit PatternPair(const HaploPattern *hpa, const HaploPattern *hpb) : m_pattern_a(*hpa), m_pattern_b(*hpb) {};
	explicit PatternPair(const HaploPattern *hpa, const HaploPattern *hpb, SP_PatternPair &pp) : m_pattern_a(*hpa), m_pattern_b(*hpb), m_prev(pp) {};

	static void *operator new(std::size_t) { return m_pool.malloc(); }
	static void operator delete(void *rawMemory) { m_pool.free(rawMemory); }

#ifdef _DEBUG
	static void *operator new(unsigned int, int, const char *, int) { return m_pool.malloc(); }
	static void operator delete(void *rawMemory, int, const char *, int) { m_pool.free(rawMemory); }
#endif // _DEBUG

	friend class HaploPair;
};


class HaploPair {
	static boost::pool<> m_pool;
	SP_PatternPair m_pair;
	int m_end, m_id[2];
	double m_likelihood;
	double m_total_likelihood;
	bool m_homogenous;
	bool m_half, m_match_a, m_match_b, m_match_next_a, m_match_next_b;

public:
	explicit HaploPair(const HaploPattern *hpa, const HaploPattern *hpb, HaploPattern *target_pattern = NULL);
	explicit HaploPair(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb);

	int end() const { return m_end; }
	double likelihood() const { return m_likelihood; }

	Genotype getGenotype();

	bool extendable(int a, int b, HaploPattern *target_pattern = NULL);
	void extend_trial(int a, int b, int id[2], double &likelihood, double &total_likelihood);

	const HaploPattern *successor_a(int i) const { return m_pair->m_pattern_a.successors(i); }
	const HaploPattern *successor_b(int i) const { return m_pair->m_pattern_b.successors(i); }

	static void getID(int id[2], const HaploPattern *hpa, const HaploPattern *hpb);

	static void *operator new(std::size_t) { return m_pool.malloc(); }
	static void operator delete(void *rawMemory) { m_pool.free(rawMemory); }

#ifdef _DEBUG
	static void *operator new(unsigned int, int, const char *, int) { return m_pool.malloc(); }
	static void operator delete(void *rawMemory, int, const char *, int) { m_pool.free(rawMemory); }
#endif // _DEBUG

	struct greater_likelihood {
		bool operator()(const HaploPair *hp1, const HaploPair *hp2) const
		{
			return hp1->likelihood() > hp2->likelihood();
		}
	};

	friend class HaploBuilder;
};


#endif // __HAPLOPAIR_H
