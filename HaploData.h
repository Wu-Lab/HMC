
#ifndef __HAPLODATA_H
#define __HAPLODATA_H


#include "Utils.h"
#include "Constant.h"
#include "Allele.h"
#include "Genotype.h"


class HaploData {
protected:
	Genotype *m_genotypes;
	int m_genotype_num;
	int m_genotype_len;
	int m_unphased_num;

	char *m_allele_type;
	int *m_allele_postition;
	int *m_allele_num;
	char (*m_allele_name)[STR_LEN_ALLELE_NAME];
	Allele **m_allele_symbol;
	double **m_allele_frequency;

public:
	HaploData();
	HaploData(const HaploData &hd);
	explicit HaploData(int num, int len);
	~HaploData();

	Genotype &operator [](int i) { return m_genotypes[i]; }
	const Genotype &operator [](int i) const { return m_genotypes[i]; }
	int genotype_num() const { return m_genotype_num; }
	int genotype_len() const { return m_genotype_len; }
	int unphased_num() const { return m_unphased_num; }
	const char *allele_types() const { return m_allele_type; }
	char allele_type(int locus) const { return m_allele_type[locus]; }
	int allele_postition(int locus) const { return m_allele_postition[locus]; }
	int allele_num(int locus) const { return m_allele_num[locus]; }
	const char *allele_name(int locus) const { return m_allele_name[locus]; }
	Allele allele_symbol(int locus, int index) const { return m_allele_symbol[locus][index]; }
	double allele_frequency(int locus, int index) const { return m_allele_frequency[locus][index]; }
	double allele_frequency(int locus, Allele a) const { return m_allele_frequency[locus][getAlleleIndex(locus, a)]; }

	int getAlleleIndex(int locus, Allele a) const;

	void setGenotypeNum(int i);
	void setGenotypeLen(int i);
	void setUnphasedNum(int i) { m_unphased_num = i; }
	void setAlleleType(int locus, char type) { m_allele_type[locus] = type; }
	void setAllelePosition(int locus, int position) { m_allele_postition[locus] = position; }
	void setAlleleName(int locus, const char *name) { strcpy(m_allele_name[locus], name); }

	void checkAlleleNum();
	void checkAlleleFrequency();
	void randomizePhase();
	void simplify();
	
	HaploData &operator =(const HaploData &hd);

protected:
	HaploData &assign(const HaploData &hd);
};

inline HaploData &HaploData::operator =(const HaploData &hd)
{
	return assign(hd);
}


#endif // __HAPLODATA_H
