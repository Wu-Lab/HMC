
#include "HaploData.h"

#include "MemLeak.h"


////////////////////////////////
//
// class HaploData

HaploData::HaploData()
: m_genotype_num(0),
  m_genotype_len(0),
  m_unphased_num(0)
{
}

HaploData::HaploData(int num, int len)
: m_genotype_num(0),
  m_genotype_len(0),
  m_unphased_num(0)
{
	setGenotypeNum(num);
	setGenotypeLen(len);
}

int HaploData::getAlleleIndex(int locus, Allele a) const
{
	int i;
	bool found;
	found = false;
	for (i=0; i<allele_num(locus); ++i) {
		if (a == allele_symbol(locus, i)) {
			found = true;
			break;
		}
	}
	if (found) {
		return i;
	}
	else {
		return -1;
	}
}

void HaploData::setGenotypeNum(int num)
{
	int i;
	if (m_genotype_num != num) {
		m_genotype_num = num > 0 ? num : 0;
		m_unphased_num = m_genotype_num;
		m_genotypes.resize(m_genotype_num);
		for (i=0; i<m_genotype_num; ++i) {
			m_genotypes[i].setLength(m_genotype_len);
		}
	}
}

void HaploData::setGenotypeLen(int len)
{
	int i;
	if (m_genotype_len != len) {
		m_genotype_len = len > 0 ? len : 0;
		for (i=0; i<m_genotype_num; ++i) {
			m_genotypes[i].setLength(m_genotype_len);
		}
		m_allele_type.resize(m_genotype_len);
		m_allele_postition.resize(m_genotype_len);
		m_allele_name.resize(m_genotype_len);
		m_allele_symbol.resize(m_genotype_len);
		for (i=0; i<m_genotype_len; ++i) {
			m_allele_type[i] = 'M';
			m_allele_postition[i] = i*Constant::average_marker_distance();
		}
	}
}

void HaploData::checkAlleleSymbol()
{
	int i, j, k, l;
	for (i=0; i<m_genotype_len; ++i) {
		m_allele_symbol[i].clear();
	}
	for (i=0; i<m_genotype_num; ++i) {
		for (j=0; j<2; ++j) {
			const Haplotype &h = m_genotypes[i](j);
			for (k=0; k<m_genotype_len; ++k) {
				if (!h[k].isMissing()) {
					if (getAlleleIndex(k, h[k]) < 0) {			// not found
						m_allele_symbol[k].push_back(make_pair(h[k], 0));
					}
				}
			}
		}
	}
	for (i=0; i<m_genotype_len; ++i) {
		sort(m_allele_symbol[i].begin(), m_allele_symbol[i].end());
	}
	vector<double> total_weight;
	total_weight.resize(m_genotype_len, 0);
	for (i=0; i<m_genotype_num; ++i) {
		for (j=0; j<2; ++j) {
			const Haplotype &h = m_genotypes[i](j);
			for (k=0; k<m_genotype_len; ++k) {
				if (!h[k].isMissing()) {
					l = getAlleleIndex(k, h[k]);
					m_allele_symbol[k][l].second += h.weight();
					total_weight[k] += h.weight();
				}
			}
		}
	}
	for (i=0; i<m_genotype_len; ++i) {
		for (j=0; j<allele_num(i); ++j) {
			m_allele_symbol[i][j].second /= total_weight[i];
		}
	}
}

void HaploData::randomizePhase()
{
	int i;
	for (i=0; i<m_unphased_num; ++i) {
		m_genotypes[i].randomizePhase();
	}
}

void HaploData::simplify()
{
	int i, j, k;
	for (i=0; i<m_genotype_num; ++i) {
		for (j=0; j<2; ++j) {
			Haplotype &h = m_genotypes[i](j);
			for (k=0; k<m_genotype_len; ++k) {
				if (!h[k].isMissing())
				{
					if (m_allele_type[k] == 'S') {
						h[k] = getAlleleIndex(k, h[k]) + '1';
					}
					else {
						h[k] = getAlleleIndex(k, h[k]) + 1;
					}
				}
			}
		}
	}
}

int HaploData::max_allele_num() const
{
	int i, m;
	m = 0;
	for (i=0; i<m_genotype_len; ++i) {
		if (allele_num(i) > m) m = allele_num(i);
	}
	return m;
}
