
#include <string.h>

#include "HaploData.h"


////////////////////////////////
//
// class HaploData

HaploData::HaploData()
{
	m_genotypes = NULL;
	m_genotype_num = 0;
	m_genotype_len = 0;
	m_unphased_num = 0;
	m_allele_type = NULL;
	m_allele_postition = NULL;
	m_allele_num = NULL;
	m_allele_name = NULL;
	m_allele_symbol = NULL;
	m_allele_frequency = NULL;
}

HaploData::HaploData(int num, int len)
{
	m_genotypes = NULL;
	m_genotype_num = 0;
	m_genotype_len = 0;
	m_unphased_num = 0;
	m_allele_type = NULL;
	m_allele_postition = NULL;
	m_allele_num = NULL;
	m_allele_name = NULL;
	m_allele_symbol = NULL;
	m_allele_frequency = NULL;
	setGenotypeNum(num);
	setGenotypeLen(len);
}

HaploData::HaploData(const HaploData &hd)
{
	m_genotypes = NULL;
	m_genotype_num = 0;
	m_genotype_len = 0;
	m_unphased_num = 0;
	m_allele_type = NULL;
	m_allele_postition = NULL;
	m_allele_num = NULL;
	m_allele_name = NULL;
	m_allele_symbol = NULL;
	m_allele_frequency = NULL;
	assign(hd);
}

HaploData::~HaploData()
{
	int i;
	delete[] m_genotypes;
	delete[] m_allele_type;
	delete[] m_allele_postition;
	delete[] m_allele_num;
	delete[] m_allele_name;
	for (i=0; i<m_genotype_len; i++) {
		delete[] m_allele_symbol[i];
		delete[] m_allele_frequency[i];
	}
	delete[] m_allele_symbol;
	delete[] m_allele_frequency;
}

int HaploData::getAlleleIndex(int locus, Allele a) const
{
	int i;
	bool found;
	found = false;
	for (i=0; i<m_allele_num[locus]; i++) {
		if (a == m_allele_symbol[locus][i]) {
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
		delete[] m_genotypes;
		m_genotype_num = num;
		m_unphased_num = m_genotype_num;
		m_genotypes = new Genotype [m_genotype_num];
	}
	for (i=0; i<m_genotype_num; i++) {
		m_genotypes[i].setLength(m_genotype_len);
	}
}

void HaploData::setGenotypeLen(int len)
{
	int i, j;
	if (m_genotype_len != len) {
		if (m_genotype_len > 0) {
			delete[] m_allele_type;
			delete[] m_allele_postition;
			delete[] m_allele_num;
			delete[] m_allele_name;
			for (i=0; i<m_genotype_len; i++) {
				delete[] m_allele_symbol[i];
				delete[] m_allele_frequency[i];
			}
			delete[] m_allele_symbol;
			delete[] m_allele_frequency;
		}
		m_genotype_len = len;
		for (i=0; i<m_genotype_num; i++) {
			m_genotypes[i].setLength(m_genotype_len);
		}
		if (m_genotype_len > 0) {
			m_allele_type = new char [m_genotype_len+1];
			m_allele_postition = new int [m_genotype_len];
			m_allele_num = new int [m_genotype_len];
			m_allele_name = new char [m_genotype_len][STR_LEN_ALLELE_NAME];
			m_allele_symbol = new Allele * [m_genotype_len];
			m_allele_frequency = new double * [m_genotype_len];
			for (i=0; i<m_genotype_len; i++) {
				m_allele_symbol[i] = new Allele [Constant::max_allele_num()];
				m_allele_frequency[i] = new double [Constant::max_allele_num()];
			}
		}
	}
	for (i=0; i<m_genotype_len; i++) {
		m_allele_type[i] = 'M';
		m_allele_postition[i] = i*Constant::average_marker_distance();
		m_allele_num[i] = 2;
		m_allele_name[i][0] = 0;
		for (j=0; j<Constant::max_allele_num(); j++) {
			m_allele_symbol[i][j] = -1;
			m_allele_frequency[i][j] = 0;
		}
	}
	m_allele_type[m_genotype_len] = 0;
}

void HaploData::checkAlleleNum()
{
	int i, j, k;
	for (i=0; i<m_genotype_len; i++) {
		m_allele_num[i] = 0;
	}
	for (i=0; i<m_genotype_num; i++) {
		for (j=0; j<2; j++) {
			const Haplotype &h = m_genotypes[i](j);
			for (k=0; k<m_genotype_len; k++) {
				if (h[k] >= 0) {
					if (getAlleleIndex(k, h[k]) < 0) {			// not found
						if (m_allele_num[k] < Constant::max_allele_num()) {
							m_allele_symbol[k][m_allele_num[k]] = h[k];
							m_allele_num[k]++;
						}
						else {
							Logger::error("The number of alleles at locus %d exceeds the capacity!", k);
							exit(1);
						}
					}
				}
			}
		}
	}
}

void HaploData::checkAlleleFrequency()
{
	int i, j, k, l;
	double *total_weight;
	total_weight = new double [m_genotype_len];
	for (i=0; i<m_genotype_len; i++) {
		total_weight[i] = 0;
		for (j=0; j<Constant::max_allele_num(); j++) {
			m_allele_frequency[i][j] = 0;
		}
	}
	for (i=0; i<m_genotype_num; i++) {
		for (j=0; j<2; j++) {
			const Haplotype &h = m_genotypes[i](j);
			for (k=0; k<m_genotype_len; k++) {
				if (h[k] >= 0) {
					l = getAlleleIndex(k, h[k]);
					m_allele_frequency[k][l] += h.m_weight;
					total_weight[k] += h.m_weight;
				}
			}
		}
	}
	for (i=0; i<m_genotype_len; i++) {
		for (j=0; j<m_allele_num[i]; j++) {
			m_allele_frequency[i][j] /= total_weight[i];
		}
	}
	delete[] total_weight;
}

void HaploData::randomizePhase()
{
	int i;
	for (i=0; i<m_unphased_num; i++) {
		m_genotypes[i].randomizePhase();
	}
}

void HaploData::simplify()
{
	int i, j, k;
	for (i=0; i<m_genotype_num; i++) {
		for (j=0; j<2; j++) {
			Haplotype &h = m_genotypes[i].haplotypes(j);
			for (k=0; k<m_genotype_len; k++) {
				if (!h.isMissing(k))
				{
					if (m_allele_type[k] == 'S') {
						h.allele(k) = getAlleleIndex(k, h[k]) + '1';
					}
					else {
						h.allele(k) = getAlleleIndex(k, h[k]) + 1;
					}
				}
			}
		}
	}
}

HaploData &HaploData::assign(const HaploData &hd)
{
	int i, j;
	setGenotypeNum(hd.m_genotype_num);
	setGenotypeLen(hd.m_genotype_len);
	for (i=0; i<m_genotype_num; i++) {
		m_genotypes[i] = hd.m_genotypes[i];
	}
	for (i=0; i<m_genotype_len; i++) {
		m_allele_type[i] = hd.m_allele_type[i];
		m_allele_postition[i] = hd.m_allele_postition[i];
		m_allele_num[i] = hd.m_allele_num[i];
		strcpy(m_allele_name[i], hd.m_allele_name[i]);
		for (j=0; j<Constant::max_allele_num(); j++) {
			m_allele_symbol[i][j] = hd.m_allele_symbol[i][j];
			m_allele_frequency[i][j] = hd.m_allele_frequency[i][j];
		}
	}
	return *this;
}
