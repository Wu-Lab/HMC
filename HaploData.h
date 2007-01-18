
#ifndef __HAPLODATA_H
#define __HAPLODATA_H


#include <vector>

#include "Utils.h"
#include "Constant.h"
#include "Allele.h"
#include "Genotype.h"


class HaploData {
protected:
	vector<Genotype> m_genotypes;
	int m_genotype_num;
	int m_genotype_len;
	int m_unphased_num;

	string m_allele_type;
	vector<int> m_allele_postition;
	vector<string> m_allele_name;
	vector<vector<pair<Allele, double> > > m_allele_symbol;

public:
	HaploData();
	explicit HaploData(int num, int len);

	Genotype &operator [](int i) { return m_genotypes[i]; }
	const Genotype &operator [](int i) const { return m_genotypes[i]; }
	int genotype_num() const { return m_genotype_num; }
	int genotype_len() const { return m_genotype_len; }
	int unphased_num() const { return m_unphased_num; }
	const string &allele_type() const { return m_allele_type; }
	char allele_type(int locus) const { return m_allele_type[locus]; }
	int allele_postition(int locus) const { return m_allele_postition[locus]; }
	const string &allele_name(int locus) const { return m_allele_name[locus]; }

	int allele_num(int locus) const { return m_allele_symbol[locus].size(); }
	int max_allele_num() const;
	Allele allele_symbol(int locus, int index) const { return m_allele_symbol[locus][index].first; }
	double allele_frequency(int locus, int index) const { return m_allele_symbol[locus][index].second; }
	double allele_frequency(int locus, Allele a) const { return allele_frequency(locus, getAlleleIndex(locus, a)); }

	int getAlleleIndex(int locus, Allele a) const;

	void setGenotypeNum(int i);
	void setGenotypeLen(int i);
	void setUnphasedNum(int i) { m_unphased_num = i; }
	void setAlleleType(int locus, char type) { m_allele_type[locus] = type; }
	void setAllelePosition(int locus, int position) { m_allele_postition[locus] = position; }
	void setAlleleName(int locus, const string &name) { m_allele_name[locus] = name; }

	void checkAlleleSymbol();
	void randomizePhase();
	void simplify();
};


#endif // __HAPLODATA_H
