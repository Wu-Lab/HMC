
#ifndef __HAPLOPATTERN_H
#define __HAPLOPATTERN_H


#include <list>
#include <boost/pool/pool.hpp>

#include "Utils.h"
#include "Allele.h"
#include "Haplotype.h"
#include "Genotype.h"
#include "GenoData.h"


class HaploPattern : public AlleleSequence, public NoThrowNewDelete {
protected:
	static boost::pool<> m_pool;
	const GenoData &m_genos;
	int m_start, m_end;
	unsigned int m_id;
	double m_frequency, m_prefix_freq;
	double m_transition_prob;
	vector<const HaploPattern*> m_successors;

public:
	explicit HaploPattern(const GenoData &genos, int start = 0);

	const GenoData &genos() const { return m_genos; }
	unsigned int id() const { return m_id; }
	int start() const { return m_start; }
	int end() const { return m_end; }
	double frequency() const { return m_frequency; }
	double prefix_freq() const { return m_prefix_freq; }
	double transition_prob() const { return m_transition_prob; }
	const HaploPattern *successors(int i) const { return i < m_successors.size() ? m_successors[i] : 0; }
	const HaploPattern *successors(Allele &a) const { return successors(m_genos.getAlleleIndex(m_end, a)); }

	int getAlleleIndex(int local_locus) const;
	int getGlobalLocus(int local_locus) const { return m_start+local_locus; }

	void setPattern(const Allele &a, int start = 0);
	void setPattern(const AlleleSequence &as, int start = 0);
	void setID(int i) { m_id = i; }
	void setFrequency(double f) { m_frequency = f; }
	void setPrefixFreq(double f) { m_prefix_freq = f; }
	void setTransitionProb(double p) { m_transition_prob = p < 1.0 ? p : 1.0; }
	void setSuccessor(size_t i, const HaploPattern *pn);

	void repack();

	bool isMatch(const Haplotype &h) const;
	bool isMatch(const Genotype &g) const;
	bool isMatch(const HaploPattern &hp) const;
	bool isMatch(const Haplotype &h, int locus) const;
	bool isMatch(const Genotype &g, int locus) const;
	bool isMatch(const HaploPattern &hp, int locus) const;
	bool isMatch(const Haplotype &h, int start, int len) const;
	bool isMatch(const Genotype &g, int start, int len) const;
	bool isMatch(const HaploPattern &hp, int start, int len) const;

	char *read(char *buffer, int len = 0);
	char *write(char *buffer, bool long_format = true) const;

	HaploPattern &assign(const HaploPattern &hp, const Allele &a);

	HaploPattern &operator +=(const Allele &a);

	using NoThrowNewDelete::operator new;
	using NoThrowNewDelete::operator delete;

	static void *operator new(std::size_t) { return m_pool.malloc(); }
	static void operator delete(void *pMemory) { m_pool.free(pMemory); }

	struct greater_frequency {
		bool operator()(const HaploPattern *hp1, const HaploPattern *hp2) const
		{
			return hp1->frequency() > hp2->frequency();
		}
	};
};

inline HaploPattern::HaploPattern(const GenoData &genos, int start)
: m_genos(genos),
  m_start(start),
  m_end(start),
  m_id(0),
  m_frequency(1.0),
  m_prefix_freq(1.0),
  m_transition_prob(1.0)
{
}

inline void HaploPattern::setSuccessor(size_t i, const HaploPattern *pn)
{
	m_successors.resize(i+1); // if (i < m_successors.size()) resize do nothing
	m_successors[i] = pn;
}

inline bool HaploPattern::isMatch(const Haplotype &h) const
{
	return h.isMatch(*this, m_start, 0, length());
}

inline bool HaploPattern::isMatch(const Genotype &g) const
{
	return g.isMatch(*this, m_start, 0, length());	
}

inline bool HaploPattern::isMatch(const HaploPattern &hp) const
{
	int start, len;
	start = m_start > hp.m_start ? m_start : hp.m_start;
	len = m_end < hp.m_end ? m_end : hp.m_end;
	len -= start;
	if (len <= 0) {
		return true;
	}
	else {
		return isMatch(hp, start, len);
	}
}

inline bool HaploPattern::isMatch(const Haplotype &h, int locus) const
{
	return m_alleles[locus-m_start].isMatch(h[locus]);	
}

inline bool HaploPattern::isMatch(const Genotype &g, int locus) const
{
	return g.isMatch(m_alleles[locus-m_start], locus);	
}

inline bool HaploPattern::isMatch(const HaploPattern &hp, int locus) const
{
	if (locus < m_start || locus >= m_end || locus < hp.m_start || locus >= hp.m_end) {
		return true;
	}
	else {
		return m_alleles[locus-m_start].isMatch(hp[locus-hp.m_start]);
	}
}

inline bool HaploPattern::isMatch(const Haplotype &h, int start, int len) const
{
	return h.isMatch(*this, start, start-m_start, len);	
}

inline bool HaploPattern::isMatch(const Genotype &g, int start, int len) const
{
	return g.isMatch(*this, start, start-m_start, len);	
}

inline bool HaploPattern::isMatch(const HaploPattern &hp, int start, int len) const
{
	return AlleleSequence::isMatch(hp, start-m_start, start-hp.m_start, len);
}

inline HaploPattern &HaploPattern::assign(const HaploPattern &hp, const Allele &a)
{
	m_start = hp.m_start;
	m_end = hp.m_end;
	m_frequency = hp.m_frequency;
	AlleleSequence::assign(hp, a);
	++m_end;
// 	checkFrequencyWithExtension(hp.m_match_frequency, m_end-1);
// 	if (length() > 1) m_transition_prob = m_frequency / hp.m_frequency;
// 	else m_transition_prob = 0;
	return *this;
}

inline HaploPattern &HaploPattern::operator +=(const Allele &a)
{
	double old_freq = m_frequency;
	AlleleSequence::operator +=(a);
	++m_end;
// 	list<pair<int, double> > mf;
// 	m_match_frequency.swap(mf);
// 	checkFrequencyWithExtension(mf, m_end-1);
// 	if (length() > 1) m_transition_prob = m_frequency / old_freq;
// 	else m_transition_prob = 0;
	return *this;
}


#endif // __HAPLOPATTERN_H
