
#include "HaploPattern.h"
#include "HaploData.h"

#include "MemLeak.h"


boost::pool<> HaploPattern::m_pool(sizeof(HaploPattern));

////////////////////////////////
//
// class HaploPattern

int HaploPattern::getAlleleIndex(int local_locus) const
{
	int index;
	if (m_alleles[local_locus].isMissing()) {
		index = -1;
	}
	else {
		index = m_haplodata.getAlleleIndex(m_start+local_locus, m_alleles[local_locus]);
	}
	return index;
}

void HaploPattern::setPattern(const Allele &a, int start)
{
	m_alleles.clear();
	m_alleles.push_back(a);
	m_start = start;
	m_end = m_start + length();
	checkFrequency();
}

void HaploPattern::setPattern(const AlleleSequence &as, int start)
{
	AlleleSequence::operator =(as);
	m_start = start;
	m_end = m_start + length();
	checkFrequency();
}

void HaploPattern::repack()
{
	int i, start, len;
	vector<Allele> temp;
	bool repacked;
	start = 0;
	repacked = false;
	for (i=0; i<length(); i++) {
		if (m_alleles[i].isMissing()) {
			m_start++;
			start++;
			repacked = true;
		}
		else {
			break;
		}
	}
	for (i=length()-1; i>=m_start; i--) {
		if (m_alleles[i].isMissing()) {
			m_end--;
			repacked = true;
		}
		else {
			break;
		}
	}
	if (repacked) {
		len = m_end - m_start;
		temp = m_alleles;
		if (len > 0) {
			m_alleles.assign(&temp[start], &temp[start+len]);
		}
		else {
			m_alleles.clear();
		}
	}
}

double HaploPattern::checkFrequency()
{
	int i;
	double w, freq;
	if (length() == 0) {
		m_frequency = m_haplodata.genotype_num();
	}
	else
	{
		m_match_frequency.clear();
		m_frequency = 0;
		for (i=0; i<m_haplodata.genotype_num(); i++) {
			const Genotype &g = m_haplodata[i];
			if (isMatch(g)) {
				freq = 1.0;
				if (g.isPhased()) {
					if (g.weight() < 1) {
						freq = getMatchingFrequency(g, &m_alleles[0], m_start, length());
						m_frequency += freq * (1 - g.weight());
					}
					w = 0;
					if (isMatch(g(0))) {
						w += 0.5;
					}
					if (isMatch(g(1))) {
						w += 0.5;
					}
					m_frequency += w * g.weight();
				}
				else {
					freq = getMatchingFrequency(g, &m_alleles[0], m_start, length());
					m_frequency += freq;
				}
				m_match_frequency.push_back(make_pair(i, freq));
			}
		}
	}
	return m_frequency;
}

double HaploPattern::checkFrequencyWithExtension(int ext, int len)
{
	double w;
	list<pair<int, double> >::iterator i_mf;
	if (length() == 0) {
		m_frequency = m_haplodata.genotype_num();
	}
	else if (m_match_frequency.empty()) {
		checkFrequency();
	}
	else
	{
		m_frequency = 0;
		i_mf = m_match_frequency.begin();
		while (i_mf != m_match_frequency.end()) {
			const Genotype &g = m_haplodata[i_mf->first];
			if (isMatch(g, ext, len)) {
				if (g.isPhased()) {
					if (g.weight() < 1) {
						i_mf->second *= getMatchingFrequency(g, &m_alleles[ext-m_start], ext, len);
						m_frequency += i_mf->second * (1 - g.weight());
					}
					w = 0;
					if (isMatch(g(0))) {
						w += 0.5;
					}
					if (isMatch(g(1))) {
						w += 0.5;
					}
					m_frequency += w * g.weight();
				}
				else {
					i_mf->second *= getMatchingFrequency(g, &m_alleles[ext-m_start], ext, len);
					m_frequency += i_mf->second;
				}
				++i_mf;
			}
			else {
				i_mf = m_match_frequency.erase(i_mf);
			}
		}
	}
	return m_frequency;
}

double HaploPattern::getMatchingFrequency(const Genotype &g, const Allele *pa, int start, int len) const
{
	double total_freq, freq;
	int i, j;
	Allele b;
	total_freq = 1.0;
	for (i=0; i<len; i++) {
		if (!pa[i].isMissing()) {			// pa[i] is not missing
			freq = 0;
			for (j=0; j<2; j++) {
				b = g(j)[start+i];
				if (b.isMissing()) {		// b is missing
					freq += m_haplodata.allele_frequency(start+i, pa[i]) > 0 ? (1.0/m_haplodata.allele_num(start+i)) : 0;
				}
				else if (b == pa[i]) {
					freq += 1;
				}
			}
			total_freq *= 0.5 * freq;
		}
	}
	return total_freq;
}

char *HaploPattern::read(char *buffer, int len)
{
	buffer = AlleleSequence::read(NULL, buffer, len);
	m_start = 0;
	m_end = length();
	repack();
	releaseMatchGenotype();
	checkFrequency();
	return buffer;
}

char *HaploPattern::write(char *buffer, bool long_format) const
{
	char *s = buffer;
	if (long_format) {
		AlleleSequence(m_start).write(NULL, s);
		s += strlen(s);
	}
	AlleleSequence::write(NULL, s);
	s += strlen(s);
	if (long_format) {
		AlleleSequence(m_haplodata.genotype_len()-m_end).write(NULL, s);
	}
	return buffer;
}
