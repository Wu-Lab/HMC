
#ifndef __GENOTYPE_H
#define __GENOTYPE_H


#include <string>

#include "Utils.h"
#include "Haplotype.h"


class Genotype {
protected:
	Haplotype m_haplotypes[2];
	string m_id;
	int m_heterozygous_num;
	int m_missing_num;
	int m_missing_allele_num;
	double m_prior_probability;
	double m_posterior_probability;
	double m_genotype_probability;
	bool m_is_phased;

public:
	Genotype();
	explicit Genotype(int len);
	explicit Genotype(const Haplotype &h1, const Haplotype &h2);

	Haplotype &operator ()(int i) { return m_haplotypes[i]; }
	const Haplotype &operator ()(int i) const { return m_haplotypes[i]; }
	const string &id() const { return m_id; }
	int length() const { return m_haplotypes[0].length(); }
	int heterozygous_num() const { return m_heterozygous_num; }
	int missing_num() const { return m_missing_num; }
	int missing_allele_num() const { return m_missing_allele_num; }
	double prior_probability() const { return m_prior_probability; }
	double posterior_probability() const { return m_posterior_probability; }
	double genotype_probability() const { return m_genotype_probability; }
	bool isPhased() const { return m_is_phased; }

	void setID(const string &id);
	void setLength(int len) { m_haplotypes[0].setLength(len); m_haplotypes[1].setLength(len); }
	void setPriorProbability(double p) { m_prior_probability = p; }
	void setPosteriorProbability(double p) { m_posterior_probability = p; }
	void setGenotypeProbability(double p) { m_genotype_probability = p; }
	void setIsPhased(bool state) { m_is_phased = state; }
	void setHaplotypes(Haplotype &h1, Haplotype &h2);
	void checkGenotype();
	void randomizePhase();

	bool isHeterozygous(int i) const { return !m_haplotypes[0][i].isMatch(m_haplotypes[1][i]); }
	int getHeterozygousNum(int start, int len) const;

	bool isMissing(int locus) const;
	bool hasMissing(int locus) const;
	bool hasAllele(int locus, Allele a) const;
	bool isMatch(const Allele &a, int locus) const;
	bool isMatch(const AlleleSequence &as, int start1, int start2, int len) const;
	bool isMatch(const Genotype &g, int locus, bool reversed = false) const;
	bool isMatch(const Genotype &g) const;
	bool isMatchUnphased(const Genotype &g) const;
	bool isMatchIgnoreMissing(const Genotype &g) const;
	int getDiffNum(const Genotype &g) const;
	int getDiffNumIgnoreMissing(const Genotype &g) const;
	int getSwitchDistance(const Genotype &g) const;
	int getSwitchDistanceIgnoreMissing(const Genotype &g) const;

	struct greater_posterior_probability {
		bool operator()(const Genotype &g1, const Genotype &g2) const
		{
			return g1.posterior_probability() > g2.posterior_probability();
		}
	};

	struct greater_prior_probability {
		bool operator()(const Genotype &g1, const Genotype &g2) const
		{
			return g1.prior_probability() > g2.prior_probability();
		}
	};
};

inline Genotype::Genotype()
: m_heterozygous_num(0),
  m_missing_num(0),
  m_missing_allele_num(0),
  m_prior_probability(0),
  m_posterior_probability(0),
  m_genotype_probability(0),
  m_is_phased(false)
{
}

inline Genotype::Genotype(int len)
: m_heterozygous_num(0),
  m_missing_num(0),
  m_missing_allele_num(0),
  m_prior_probability(0),
  m_posterior_probability(0),
  m_genotype_probability(0),
  m_is_phased(false)
{
  m_haplotypes[0].setLength(len);
  m_haplotypes[1].setLength(len);
}

inline void Genotype::setID(const string &id)
{
	m_id = id;
	string_replace(m_id, " ", "_");
	m_haplotypes[0].setID(m_id);
	m_haplotypes[1].setID(m_id);
}

inline bool Genotype::isMissing(int locus) const
{
	return (m_haplotypes[0][locus].isMissing() && m_haplotypes[1][locus].isMissing());
}

inline bool Genotype::hasMissing(int locus) const
{
	return (m_haplotypes[0][locus].isMissing() || m_haplotypes[1][locus].isMissing());
}

inline bool Genotype::hasAllele(int locus, Allele a) const
{
	return (m_haplotypes[0][locus] == a || m_haplotypes[1][locus] == a);
}

inline bool Genotype::isMatch(const Allele &a, int locus) const
{
	return (m_haplotypes[0][locus].isMatch(a) || m_haplotypes[1][locus].isMatch(a));
}

inline bool Genotype::isMatchIgnoreMissing(const Genotype &g) const
{
	return (getSwitchDistanceIgnoreMissing(g) == 0);
}


#endif // __GENOTYPE_H
