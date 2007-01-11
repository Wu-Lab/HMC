
#ifndef __HAPLOPATTERN_H
#define __HAPLOPATTERN_H


#include "Utils.h"
#include "Allele.h"
#include "Haplotype.h"
#include "Genotype.h"


class HaploData;


class HaploPattern : public AlleleSequence {
protected:
	const HaploData *m_haplo_data;
	unsigned int m_id;
	int m_start, m_end;
	double m_frequency;
	double m_weight;
	int m_genotype_num;
	bool *m_match_genotype;
	double *m_match_frequency;
	int m_successors_num;
	HaploPattern **m_successors;
	double m_transition_prob;
	HaploPattern *m_prefix;

public:
	HaploPattern();
	HaploPattern(const HaploPattern &hp);
	explicit HaploPattern(const HaploData *haplo, int start = 0);
	explicit HaploPattern(const HaploData *haplo, const Allele &a, int start = 0);
	explicit HaploPattern(const HaploData *haplo, const AlleleSequence &as, int start = 0);
	~HaploPattern();

	const HaploData *haplo_data() { return m_haplo_data; }
	unsigned int id() { return m_id; }
	int start() const { return m_start; }
	int end() const { return m_end; }
	double frequency() const { return m_frequency; }
	double weight() const { return m_weight; }
	HaploPattern *successors(int i) const { return m_successors[i]; }

	int getAlleleIndex(int local_locus) const;
	int getGlobalLocus(int local_locus) const { return m_start+local_locus; }

	void setHaploData(const HaploData *haplo);
	void setPattern(const Allele &a, int start = 0);
	void setPattern(const AlleleSequence &as, int start = 0);
	void setSuccessor(int i, HaploPattern *pn) { m_successors[i] = pn; }

	void repack();

	void copyMatchGenotype(const HaploPattern &hp);
	void releaseMatchGenotype();
	double checkFrequency();
	double checkFrequencyWithExtension(int ext);
	double checkFrequencyWithExtension(int ext, int len);

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

	HaploPattern &operator =(const HaploPattern &hp);
	HaploPattern &operator +=(const HaploPattern &hp);
	HaploPattern &operator +=(const AlleleSequence &as);
	HaploPattern &operator +=(const Allele &a);

protected:
	double getMatchingFrequency(int h, const Allele *a, int start, int len) const;

	HaploPattern &assign(const HaploPattern &hp);
	HaploPattern &concatenate(const HaploPattern &hp);
	HaploPattern &concatenate(const AlleleSequence &as);
	HaploPattern &concatenate(const Allele &a);

public:
	friend HaploPattern operator +(const HaploPattern &hp1, const HaploPattern &hp2);
	friend HaploPattern operator +(const HaploPattern &hp, const AlleleSequence &as);
	friend HaploPattern operator +(const AlleleSequence &as, const HaploPattern &hp);
	friend HaploPattern operator +(const HaploPattern &hp, const Allele &a);
	friend HaploPattern operator +(const Allele &a, const HaploPattern &hp);

	HaploPattern &concatenate(const HaploPattern &hp1, const HaploPattern &hp2);
	HaploPattern &concatenate(const HaploPattern &hp, const AlleleSequence &as);
	HaploPattern &concatenate(const AlleleSequence &as, const HaploPattern &hp);
	HaploPattern &concatenate(const HaploPattern &hp, const Allele &a);
	HaploPattern &concatenate(const Allele &a, const HaploPattern &hp);

	friend class HaploPair;
	friend class PatternTree;
	friend class HaploBuilder;
	friend class HaploData;
};

inline bool HaploPattern::isMatch(const Haplotype &h) const
{
	return h.isMatch(*this, m_start, 0, m_length);
}

inline bool HaploPattern::isMatch(const Genotype &g) const
{
	return g.isMatch(*this, m_start, 0, m_length);	
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
	return h.isMatch(m_alleles[locus-m_start], locus);	
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
		return AlleleSequence::isMatch(hp.m_alleles[locus-hp.m_start], locus-m_start);
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

inline double HaploPattern::checkFrequencyWithExtension(int ext)
{
	return checkFrequencyWithExtension(ext, 1);
}

inline HaploPattern &HaploPattern::operator =(const HaploPattern &hp)
{
	return assign(hp);
}

inline HaploPattern &HaploPattern::operator +=(const HaploPattern &hp)
{
	return concatenate(hp);
}

inline HaploPattern &HaploPattern::operator +=(const AlleleSequence &as)
{
	return concatenate(as);
}

inline HaploPattern &HaploPattern::operator +=(const Allele &a)
{
	return concatenate(a);
}

inline HaploPattern operator +(const HaploPattern &hp1, const HaploPattern &hp2)
{
	return HaploPattern().concatenate(hp1, hp2);
}

inline HaploPattern operator +(const HaploPattern &hp, const AlleleSequence &as)
{
	return HaploPattern().concatenate(hp, as);
}

inline HaploPattern operator +(const AlleleSequence &as, const HaploPattern &hp)
{
	return HaploPattern().concatenate(as, hp);
}

inline HaploPattern operator +(const HaploPattern &hp, const Allele &a)
{
	return HaploPattern().concatenate(hp, a);
}

inline HaploPattern operator +(const Allele &a, const HaploPattern &hp)
{
	return HaploPattern().concatenate(a, hp);
}


#endif // __HAPLOPATTERN_H
