
#ifndef __HAPLOPAIR_H
#define __HAPLOPAIR_H


#include <boost/pool/pool.hpp>

#include "Utils.h"
#include "HaploPattern.h"


class HaploPair;

struct HaploPairLink {
	HaploPair *link;
	int index;
	bool reversed;
	double likelihood;

	HaploPairLink() : link(0), index(0), reversed(false), likelihood(0) {}
	HaploPairLink(HaploPair *h, int i, bool r, double l);

	void set(HaploPair *h, int i, bool r, double l);

	friend bool operator <(const HaploPairLink &lhs, const HaploPairLink &rhs);
	friend bool operator >(const HaploPairLink &lhs, const HaploPairLink &rhs);
};

inline HaploPairLink::HaploPairLink(HaploPair *h, int i, bool r, double l)
: link(h), index(i), reversed(r), likelihood(l)
{
}

inline void HaploPairLink::set(HaploPair *h, int i, bool r, double l)
{
	link = h;
	index = i;
	reversed = r;
	likelihood = l;
}

inline bool operator <(const HaploPairLink &lhs, const HaploPairLink &rhs)
{
	return (lhs.likelihood < rhs.likelihood);
}

inline bool operator >(const HaploPairLink &lhs, const HaploPairLink &rhs)
{
	return (lhs.likelihood > rhs.likelihood);
}


class HaploPair : public NoThrowNewDelete {
	static boost::pool<> m_pool;

	const HaploPattern &m_pattern_a, &m_pattern_b;
	const Allele &m_allele_a, &m_allele_b;
	vector<HaploPair*> m_forward_links[2];
	vector<HaploPairLink> m_best_links;
	double m_transition_prob;
	double m_forward_likelihood;
	double m_backward_likelihood;

public:
	explicit HaploPair(const HaploPattern *hpa, const HaploPattern *hpb);
	explicit HaploPair(const HaploPattern *hpa, const HaploPattern *hpb, HaploPair *hp, bool reversed = false);
	void add(HaploPair *hp, bool reversed = false, int best_num = 1);

	const HaploPattern &pattern_a() const { return m_pattern_a; }
	const HaploPattern &pattern_b() const { return m_pattern_b; }
	const Allele &allele_a() const { return m_allele_a; }
	const Allele &allele_b() const { return m_allele_b; }
	int id_a() const { return m_pattern_a.id(); }
	int id_b() const { return m_pattern_b.id(); }
	int end() const { return m_pattern_a.end(); }

	const vector<HaploPair*> &forward_links(int i) const { return m_forward_links[i]; }
	double transition_prob() const { return m_transition_prob; }
	double forward_likelihood() const { return m_forward_likelihood; }
	double backward_likelihood() const { return m_backward_likelihood; }

	Genotype getGenotype(int index = 0) const;
	double getLikelihood(int index = 0) const;

	template <typename T>
	const HaploPattern *successor_a(T &i) const { return m_pattern_a.successors(i); }
	template <typename T>
	const HaploPattern *successor_b(T &i) const { return m_pattern_b.successors(i); }

	void sortBestLinks();
	void calcBackwardLikelihood();

	using NoThrowNewDelete::operator new;
	using NoThrowNewDelete::operator delete;

	static void *operator new(std::size_t) { return m_pool.malloc(); }
	static void operator delete(void *pMemory) { m_pool.free(pMemory); }

	struct greater_likelihood {
		bool operator()(const HaploPair *hp1, const HaploPair *hp2) const
		{
			return hp1->getLikelihood() > hp2->getLikelihood();
		}
	};

private:
	HaploPair(const HaploPair &);
	HaploPair &operator=(const HaploPair &);
};

inline double HaploPair::getLikelihood(int index) const
{
	if ((index >= 0) && (index < m_best_links.size())) {
		return m_best_links[index].likelihood;
	}
	return 0;
}

inline void HaploPair::sortBestLinks()
{
	sort(m_best_links.begin(), m_best_links.end(), greater<HaploPairLink>());
}


#endif // __HAPLOPAIR_H
