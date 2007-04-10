
#include "HaploData.h"

#include "MemLeak.h"


////////////////////////////////
//
// class HaploData

HaploData::HaploData()
: m_haplotype_num(0),
  m_haplotype_len(0)
{
}

HaploData::HaploData(int num, int len)
: m_haplotype_num(0),
  m_haplotype_len(0)
{
	setHaplotypeNum(num);
	setHaplotypeLen(len);
}

HaploData::HaploData(const GenoData &genos)
: m_haplotype_num(genos.genotype_num() * 2),
  m_haplotype_len(genos.genotype_len()),
  m_allele_type(genos.m_allele_type),
  m_allele_postition(genos.m_allele_postition),
  m_allele_name(genos.m_allele_name),
  m_allele_symbol(genos.m_allele_symbol)
{
	m_haplotypes.resize(m_haplotype_num);
	for (int i=0; i<genos.genotype_num(); ++i) {
		m_haplotypes[2*i] = genos[i](0);
		m_haplotypes[2*i+1] = genos[i](1);
	}
}

HaploData &HaploData::operator =(const GenoData &genos)
{
	m_haplotype_num = genos.genotype_num() * 2;
	m_haplotype_len = genos.genotype_len();
	m_allele_type = genos.m_allele_type;
	m_allele_postition = genos.m_allele_postition;
	m_allele_name = genos.m_allele_name;
	m_allele_symbol = genos.m_allele_symbol;
	m_haplotypes.resize(m_haplotype_num);
	for (int i=0; i<genos.genotype_num(); ++i) {
		m_haplotypes[2*i] = genos[i](0);
		m_haplotypes[2*i+1] = genos[i](1);
	}
	return *this;
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

void HaploData::setHaplotypeNum(int num)
{
	int i;
	if (m_haplotype_num != num) {
		m_haplotype_num = num > 0 ? num : 0;
		m_haplotypes.resize(m_haplotype_num);
		for (i=0; i<m_haplotype_num; ++i) {
			m_haplotypes[i].setLength(m_haplotype_len);
		}
	}
}

void HaploData::setHaplotypeLen(int len)
{
	int i;
	if (m_haplotype_len != len) {
		m_haplotype_len = len > 0 ? len : 0;
		for (i=0; i<m_haplotype_num; ++i) {
			m_haplotypes[i].setLength(m_haplotype_len);
		}
		m_allele_type.resize(m_haplotype_len);
		m_allele_postition.resize(m_haplotype_len);
		m_allele_name.resize(m_haplotype_len);
		m_allele_symbol.resize(m_haplotype_len);
		for (i=0; i<m_haplotype_len; ++i) {
			m_allele_type[i] = 'M';
			m_allele_postition[i] = i*Constant::average_marker_distance();
		}
	}
}

void HaploData::addHaplotype(const Haplotype &haplo)
{
	m_haplotypes.push_back(haplo);
	m_haplotype_num++;
}

void HaploData::addHaplotype(const vector<Haplotype> &haplos)
{
	m_haplotypes.insert(m_haplotypes.end(), haplos.begin(), haplos.end());
	m_haplotype_num += haplos.size();
}

void HaploData::checkTotalWeight()
{
	m_total_weight = 0;
	for (int i=0; i<m_haplotype_num; ++i) {
		m_total_weight += m_haplotypes[i].weight();
	}
}

void HaploData::checkAlleleSymbol()
{
	int i, j, k, l;
	for (i=0; i<m_haplotype_len; ++i) {
		m_allele_symbol[i].clear();
	}
	for (i=0; i<m_haplotype_num; ++i) {
		const Haplotype &h = m_haplotypes[i];
		for (k=0; k<m_haplotype_len; ++k) {
			if (!h[k].isMissing()) {
				if (getAlleleIndex(k, h[k]) < 0) {			// not found
					m_allele_symbol[k].push_back(make_pair(h[k], 0));
				}
			}
		}
	}
	for (i=0; i<m_haplotype_len; ++i) {
		sort(m_allele_symbol[i].begin(), m_allele_symbol[i].end());
	}
	vector<double> total_weight;
	total_weight.resize(m_haplotype_len, 0);
	for (i=0; i<m_haplotype_num; ++i) {
		const Haplotype &h = m_haplotypes[i];
		for (k=0; k<m_haplotype_len; ++k) {
			if (!h[k].isMissing()) {
				l = getAlleleIndex(k, h[k]);
				m_allele_symbol[k][l].second += h.weight();
				total_weight[k] += h.weight();
			}
		}
	}
	for (i=0; i<m_haplotype_len; ++i) {
		for (j=0; j<allele_num(i); ++j) {
			m_allele_symbol[i][j].second /= total_weight[i];
		}
	}
}

void HaploData::simplify()
{
	int i, k;
	for (i=0; i<m_haplotype_num; ++i) {
		Haplotype &h = m_haplotypes[i];
		for (k=0; k<m_haplotype_len; ++k) {
			if (!h[k].isMissing()) {
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

int HaploData::max_allele_num() const
{
	int i, m;
	m = 0;
	for (i=0; i<m_haplotype_len; ++i) {
		if (allele_num(i) > m) m = allele_num(i);
	}
	return m;
}
