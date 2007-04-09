
#ifndef __HAPLODATA_H
#define __HAPLODATA_H


#include <vector>

#include "Utils.h"
#include "Constant.h"
#include "Allele.h"
#include "Haplotype.h"
#include "GenoData.h"


class HaploData {
protected:
	vector<Haplotype> m_haplotypes;
	int m_haplotype_num;
	int m_haplotype_len;
	double m_total_weight;
	
	string m_allele_type;
	vector<int> m_allele_postition;
	vector<string> m_allele_name;
	vector<vector<pair<Allele, double> > > m_allele_symbol;

public:
	HaploData();
	explicit HaploData(int num, int len);
	explicit HaploData(const GenoData &genos);
	HaploData &operator =(const GenoData &genos);

	Haplotype &operator [](int i) { return m_haplotypes[i]; }
	const Haplotype &operator [](int i) const { return m_haplotypes[i]; }
	int haplotype_num() const { return m_haplotype_num; }
	int haplotype_len() const { return m_haplotype_len; }
	double total_weight() const { return m_total_weight; }
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

	void setHaplotypeNum(int i);
	void setHaplotypeLen(int i);
	void setAlleleType(int locus, char type) { m_allele_type[locus] = type; }
	void setAllelePosition(int locus, int position) { m_allele_postition[locus] = position; }
	void setAlleleName(int locus, const string &name) { m_allele_name[locus] = name; string_replace(m_allele_name[locus], " ", "_"); }

	void addHaplotype(const Haplotype &haplo);
	void addHaplotype(const vector<Haplotype> &haplos);

	void checkTotalWeight();
	void checkAlleleSymbol();
	void simplify();
};


#endif // __HAPLODATA_H
