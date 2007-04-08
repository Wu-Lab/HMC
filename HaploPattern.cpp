
#include "HaploPattern.h"
#include "GenoData.h"

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
		index = m_genos.getAlleleIndex(m_start+local_locus, m_alleles[local_locus]);
	}
	return index;
}

void HaploPattern::setPattern(const Allele &a, int start)
{
	m_alleles.clear();
	m_alleles.push_back(a);
	m_start = start;
	m_end = m_start + length();
}

void HaploPattern::setPattern(const AlleleSequence &as, int start)
{
	AlleleSequence::operator =(as);
	m_start = start;
	m_end = m_start + length();
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
	for (i=length()-1; i>=m_start; --i) {
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

char *HaploPattern::read(char *buffer, int len)
{
	buffer = AlleleSequence::read(NULL, buffer, len);
	m_start = 0;
	m_end = length();
	repack();
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
		AlleleSequence(m_genos.genotype_len()-m_end).write(NULL, s);
	}
	return buffer;
}
