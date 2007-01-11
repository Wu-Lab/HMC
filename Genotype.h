
#ifndef __GENOTYPE_H
#define __GENOTYPE_H


#include "Utils.h"
#include "Haplotype.h"


class Genotype {
protected:
	Haplotype m_haplotypes[2];
	char m_id[STR_LEN_GENOTYPE_ID];
	int m_length;
	int m_heterozygous_num;
	int m_missing_num;
	int m_missing_allele_num;
	double m_likelihood;
	double m_weight;
	bool m_is_phased;

public:
	Genotype();
	Genotype(const Genotype &g);
	explicit Genotype(int len);
	explicit Genotype(const Haplotype &h1, const Haplotype &h2);

	Haplotype &operator ()(int i) { return m_haplotypes[i]; }
	const Haplotype &operator ()(int i) const { return m_haplotypes[i]; }
	char *id() const { return (char *) m_id; }
	int length() const { return m_length; }
	int heterozygous_num() const { return m_heterozygous_num; }
	int missing_num() const { return m_missing_num; }
	int missing_allele_num() const { return m_missing_allele_num; }
	double likelihood() const { return m_likelihood; }
	double weight() const { return m_weight; }
	bool isPhased() const { return m_is_phased; }

	void setID(const char *id) { strcpy(m_id, id); }
	void setLength(int len) { m_length = len; }
	void setLikelihood(double likelihood) { m_likelihood = likelihood; }
	void setWeight(double weight) { m_weight = weight; }
	void setIsPhased(bool state) { m_is_phased = state; }
	void setHaplotypes(Haplotype &h1, Haplotype &h2);
	void checkGenotype();
	void randomizePhase();

	bool isHeterozygous(int i) const { return !m_haplotypes[0].isMatch(m_haplotypes[1][i], i); }
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

	Genotype &operator =(const Genotype &g);
	Genotype &operator +=(const Genotype &g);

protected:
	Genotype &assign(const Genotype &g);
	Genotype &concatenate(const Genotype &g);

public:
	friend Genotype operator +(const Genotype &g1, const Genotype &g2);

	Genotype &concatenate(const Genotype &g1, const Genotype &g2);

	friend class HaploPattern;
	friend class HaploData;
};

inline bool Genotype::isMissing(int locus) const
{
	return (m_haplotypes[0].isMissing(locus) && m_haplotypes[1].isMissing(locus));
}

inline bool Genotype::hasMissing(int locus) const
{
	return (m_haplotypes[0].isMissing(locus) || m_haplotypes[1].isMissing(locus));
}

inline bool Genotype::hasAllele(int locus, Allele a) const
{
	return (m_haplotypes[0][locus] == a || m_haplotypes[1][locus] == a);
}

inline bool Genotype::isMatch(const Allele &a, int locus) const
{
	return (m_haplotypes[0].isMatch(a, locus) || m_haplotypes[1].isMatch(a, locus));
}

inline bool Genotype::isMatchIgnoreMissing(const Genotype &g) const
{
	return (getSwitchDistanceIgnoreMissing(g) == 0);
}

inline Genotype &Genotype::operator =(const Genotype &g)
{
	return assign(g);
}

inline Genotype &Genotype::operator +=(const Genotype &g)
{
	return concatenate(g);
}

inline Genotype operator +(const Genotype &g1, const Genotype &g2)
{
	return Genotype().concatenate(g1, g2);
}


#endif // __GENOTYPE_H
