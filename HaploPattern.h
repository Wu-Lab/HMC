
#ifndef __HAPLOPATTERN_H
#define __HAPLOPATTERN_H


#include <list>
#include <boost/pool/pool.hpp>

#include "Utils.h"
#include "Allele.h"
#include "Haplotype.h"
#include "Genotype.h"


class HaploData;


class HaploPattern : public AlleleSequence {
protected:
	static boost::pool<> m_pool;
	const HaploData &m_haplodata;
	int m_start, m_end;
	unsigned int m_id;
	double m_frequency;
	double m_transition_prob;
	const HaploPattern *m_prefix;
	vector<const HaploPattern*> m_successors;
	list<pair<int, double> > m_match_frequency;

public:
	explicit HaploPattern(const HaploData &hd, int start);

	const HaploData &haplodata() const { return m_haplodata; }
	unsigned int id() const { return m_id; }
	int start() const { return m_start; }
	int end() const { return m_end; }
	double frequency() const { return m_frequency; }
	double transition_prob() const { return m_transition_prob; }
	const HaploPattern *prefix() const { return m_prefix; }
	const HaploPattern *successors(int i) const { return m_successors[i]; }

	int getAlleleIndex(int local_locus) const;
	int getGlobalLocus(int local_locus) const { return m_start+local_locus; }

	void setPattern(const Allele &a, int start = 0);
	void setPattern(const AlleleSequence &as, int start = 0);
	void setID(int i) { m_id = i; }
	void setFrequency(double f) { m_frequency = f; }
	void setTransitionProb(double p) { m_transition_prob = p; }
	void setPrefix(const HaploPattern *p) { m_prefix = p; }
	void setSuccessor(int i, const HaploPattern *pn) { m_successors.resize(i+1); m_successors[i] = pn; }

	void repack();

	void releaseMatchGenotype() { m_match_frequency.clear(); }
	double checkFrequency();
	double checkFrequencyWithExtension(int ext, int len = 1);

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

	HaploPattern &operator +=(const AlleleSequence &as);
	HaploPattern &operator +=(const Allele &a);

	static void *operator new(std::size_t) { return m_pool.malloc(); }
	static void operator delete(void *rawMemory) { m_pool.free(rawMemory); }

#ifdef _DEBUG
	static void *operator new(unsigned int, int, const char *, int) { return m_pool.malloc(); }
	static void operator delete(void *rawMemory, int, const char *, int) { m_pool.free(rawMemory); }
#endif // _DEBUG

protected:
	double getMatchingFrequency(const Genotype &g, const Allele *pa, int start, int len) const;
};

inline HaploPattern::HaploPattern(const HaploData &hd, int start)
: m_haplodata(hd),
  m_start(start),
  m_end(start),
  m_id(0),
  m_frequency(0),
  m_transition_prob(1.0),
  m_prefix(0)
{
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
	m_match_frequency = hp.m_match_frequency;
	AlleleSequence::assign(hp, a);
	++m_end;
	checkFrequencyWithExtension(m_end-1);
	return *this;
}

inline HaploPattern &HaploPattern::operator +=(const AlleleSequence &as)
{
	AlleleSequence::operator +=(as);
	m_end += as.length();
	checkFrequencyWithExtension(m_end-as.length(), as.length());
	return *this;
}

inline HaploPattern &HaploPattern::operator +=(const Allele &a)
{
	AlleleSequence::operator +=(a);
	++m_end;
	checkFrequencyWithExtension(m_end-1);
	return *this;
}


#endif // __HAPLOPATTERN_H
