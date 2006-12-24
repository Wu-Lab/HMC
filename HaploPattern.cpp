
#include <string.h>

#include "HaploPattern.h"


////////////////////////////////
//
// class HaploPattern

HaploPattern::HaploPattern()
{
	int i;
	m_haplo_data = NULL;
	m_id = 0;
	m_start = m_end = 0;
	m_frequency = 0;
	m_weight = 0;
	m_genotype_num = 0;
	m_match_genotype = NULL;
	m_match_frequency = NULL;
	m_successors_num = 0;
	m_successors = new HaploPattern * [Constant::max_allele_num()];
	for (i=0; i<Constant::max_allele_num(); i++) {
		m_successors[i] = NULL;
	}
	m_transition_prob = 0;
	m_prefix = NULL;
}

HaploPattern::HaploPattern(const HaploData *haplo, int start)
{
	int i;
	m_haplo_data = haplo;
	m_id = 0;
	m_start = m_end = start;
	m_frequency = 0;
	m_genotype_num = m_haplo_data->m_genotype_num;
	m_match_genotype = NULL;
	m_match_frequency = NULL;
	m_successors_num = m_haplo_data->m_allele_num[m_end];
	m_successors = new HaploPattern * [Constant::max_allele_num()];
	for (i=0; i<Constant::max_allele_num(); i++) {
		m_successors[i] = NULL;
	}
	m_transition_prob = 0;
	m_prefix = NULL;
	checkFrequency();
}

HaploPattern::HaploPattern(const HaploData *haplo, const Allele &a, int start)
	: AlleleSequence(a)
{
	int i;
	m_haplo_data = haplo;
	m_id = 0;
	m_start = start;
	m_end = m_start + m_length;
	m_frequency = 0;
	m_genotype_num = m_haplo_data->m_genotype_num;
	m_match_genotype = NULL;
	m_match_frequency = NULL;
	m_successors_num = m_haplo_data->m_allele_num[m_end];
	m_successors = new HaploPattern * [Constant::max_allele_num()];
	for (i=0; i<Constant::max_allele_num(); i++) {
		m_successors[i] = NULL;
	}
	m_transition_prob = 0;
	m_prefix = NULL;
	checkFrequency();
}

HaploPattern::HaploPattern(const HaploData *haplo, const AlleleSequence &as, int start)
	: AlleleSequence(as)
{
	int i;
	m_haplo_data = haplo;
	m_id = 0;
	m_start = start;
	m_end = m_start + m_length;
	m_frequency = 0;
	m_genotype_num = m_haplo_data->m_genotype_num;
	m_match_genotype = NULL;
	m_match_frequency = NULL;
	m_successors_num = m_haplo_data->m_allele_num[m_end];
	m_successors = new HaploPattern * [Constant::max_allele_num()];
	for (i=0; i<Constant::max_allele_num(); i++) {
		m_successors[i] = NULL;
	}
	m_transition_prob = 0;
	m_prefix = NULL;
	checkFrequency();
}

HaploPattern::HaploPattern(const HaploPattern &hp)
	: AlleleSequence(hp)
{
	int i;
	m_haplo_data = NULL;
	m_match_genotype = NULL;
	m_match_frequency = NULL;
	m_successors = new HaploPattern * [Constant::max_allele_num()];
	for (i=0; i<Constant::max_allele_num(); i++) {
		m_successors[i] = NULL;
	}
	m_prefix = NULL;
	assign(hp);
}

HaploPattern::~HaploPattern()
{
	delete[] m_match_genotype;
	delete[] m_match_frequency;
	delete[] m_successors;
}

int HaploPattern::getAlleleIndex(int local_locus) const
{
	int index;
	if (isMissing(local_locus)) {
		index = -1;
	}
	else {
		index = m_haplo_data->getAlleleIndex(m_start+local_locus, m_alleles[local_locus]);
	}
	return index;
}

void HaploPattern::setHaploData(const HaploData *haplo)
{
	m_haplo_data = haplo;
	if (m_genotype_num != m_haplo_data->m_genotype_num) {
		delete[] m_match_genotype;
		delete[] m_match_frequency;
		m_genotype_num = m_haplo_data->m_genotype_num;
		m_match_genotype = NULL;
		m_match_frequency = NULL;
	}
}

void HaploPattern::setPattern(const Allele &a, int start)
{
	AlleleSequence::assign(a);
	m_start = start;
	m_end = m_start + m_length;
	checkFrequency();
}

void HaploPattern::setPattern(const AlleleSequence &as, int start)
{
	AlleleSequence::assign(as);
	m_start = start;
	m_end = m_start + m_length;
	checkFrequency();
}

void HaploPattern::repack()
{
	int i, start;
	Allele *temp;
	bool repacked;
	start = 0;
	repacked = false;
	for (i=0; i<m_length; i++) {
		if (isMissing(i)) {
			m_start++;
			start++;
			repacked = true;
		}
		else {
			break;
		}
	}
	for (i=m_length-1; i>=m_start; i--) {
		if (isMissing(i)) {
			m_end--;
			repacked = true;
		}
		else {
			break;
		}
	}
	if (repacked) {
		m_length = m_end - m_start;
		temp = m_alleles;
		if (m_length > 0) {
			m_alleles = new Allele [m_length];
			for (i=0; i<m_length; i++) {
				m_alleles[i] = temp[start+i];
			}
		}
		else {
			m_alleles = NULL;
		}
		delete[] temp;
	}
}

void HaploPattern::copyMatchGenotype(const HaploPattern &hp)
{
	int i;
	m_haplo_data = hp.m_haplo_data;
	if (hp.m_match_genotype == NULL) {
		m_genotype_num = hp.m_genotype_num;
		delete[] m_match_genotype;
		delete[] m_match_frequency;
		m_match_genotype = NULL;
		m_match_frequency = NULL;
	}
	else {
		if (m_genotype_num != hp.m_genotype_num) {
			m_genotype_num = hp.m_genotype_num;
			delete[] m_match_genotype;
			delete[] m_match_frequency;
			m_match_genotype = new bool [m_genotype_num];
			m_match_frequency = new double [m_genotype_num];
		}
		for (i=0; i<m_genotype_num; i++) {
			m_match_genotype[i] = hp.m_match_genotype[i];
			m_match_frequency[i] = hp.m_match_frequency[i];
		}
	}
}

void HaploPattern::releaseMatchGenotype()
{
	delete[] m_match_genotype;
	delete[] m_match_frequency;
	m_match_genotype = NULL;
	m_match_frequency = NULL;
}

double HaploPattern::checkFrequency()
{
	int i;
	double w;
	if (m_haplo_data != NULL) {
		if (m_length == 0) {
			m_frequency = m_genotype_num;
		}
		else
		{
			if (m_match_genotype == NULL) {
				m_match_genotype = new bool [m_genotype_num];
				m_match_frequency = new double [m_genotype_num];
			}
			m_frequency = 0;
			m_weight = 0;
			for (i=0; i<m_genotype_num; i++) {
				m_match_genotype[i] = isMatch(m_haplo_data->m_genotypes[i]);
				m_match_frequency[i] = 0;
				if (m_match_genotype[i]) {
					if (m_haplo_data->m_genotypes[i].isPhased()) {
						if (m_haplo_data->m_genotypes[i].weight() < 1) {
							m_match_frequency[i] = getMatchingFrequency(i, m_alleles, m_start, m_length);
							m_frequency += m_match_frequency[i] * (1 - m_haplo_data->m_genotypes[i].weight());
						}
						w = 0;
						if (isMatch(m_haplo_data->m_genotypes[i](0))) {
							w += 0.5;
						}
						if (isMatch(m_haplo_data->m_genotypes[i](1))) {
							w += 0.5;
						}
						m_frequency += w * m_haplo_data->m_genotypes[i].weight();
					}
					else {
						m_match_frequency[i] = getMatchingFrequency(i, m_alleles, m_start, m_length);
						m_frequency += m_match_frequency[i];
					}
				}
			}
		}
	}
	return m_frequency;
}

double HaploPattern::checkFrequencyWithExtension(int ext, int len)
{
	int i;
	double w;
	if (m_haplo_data != NULL) {
		if (m_length == 0) {
			m_frequency = m_genotype_num;
		}
		else if (m_match_genotype == NULL) {
			checkFrequency();
		}
		else
		{
			m_frequency = 0;
			m_weight = 0;
			for (i=0; i<m_genotype_num; i++) {
				if (m_match_genotype[i]) {
					m_match_genotype[i] = isMatch(m_haplo_data->m_genotypes[i], ext, len);
					if (m_match_genotype[i]) {
						if (m_haplo_data->m_genotypes[i].isPhased()) {
							if (m_haplo_data->m_genotypes[i].weight() < 1) {
								m_match_frequency[i] *= getMatchingFrequency(i, m_alleles+ext-m_start, ext, len);
								m_frequency += m_match_frequency[i] * (1 - m_haplo_data->m_genotypes[i].weight());
							}
							w = 0;
							if (isMatch(m_haplo_data->m_genotypes[i](0))) {
								w += 0.5;
							}
							if (isMatch(m_haplo_data->m_genotypes[i](1))) {
								w += 0.5;
							}
							m_frequency += w * m_haplo_data->m_genotypes[i].weight();
						}
						else {
							m_match_frequency[i] *= getMatchingFrequency(i, m_alleles+ext-m_start, ext, len);
							m_frequency += m_match_frequency[i];
						}
					}
				}
			}
		}
	}
	return m_frequency;
}

double HaploPattern::getMatchingFrequency(int g, const Allele *a, int start, int len) const
{
	double total_freq, freq;
	int i, j;
	Allele b;
	total_freq = 1.0;
	for (i=0; i<len; i++) {
		if (a[i] >= 0) {			// a[i] is not missing
			freq = 0;
			for (j=0; j<2; j++) {
				b = m_haplo_data->m_genotypes[g](j)[start+i];
				if (b < 0) {			// b is missing
					freq += m_haplo_data->allele_frequency(start+i, a[i]) > 0 ? (1.0/m_haplo_data->allele_num(start+i)) : 0;
				}
				else if (b == a[i]) {
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
	m_end = m_length;
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
	if (long_format && m_haplo_data != NULL) {
		AlleleSequence(m_haplo_data->m_genotype_len-m_end).write(NULL, s);
	}
	return buffer;
}

HaploPattern &HaploPattern::assign(const HaploPattern &hp)
{
	int i;
	AlleleSequence::assign(hp);
	m_id = hp.m_id;
	m_start = hp.m_start;
	m_end = hp.m_end;
	m_frequency = hp.m_frequency;
	copyMatchGenotype(hp);
	m_successors_num = hp.m_successors_num;
	for (i=0; i<Constant::max_allele_num(); i++) {
		m_successors[i] = hp.m_successors[i];
	}
	m_transition_prob = 0;
	m_prefix = NULL;
	return *this;
}

HaploPattern &HaploPattern::concatenate(const HaploPattern &hp)
{
	if (m_end != hp.m_start) {
		Logger::warning("[HaploPattern::concatenate] Attempt to concatenate incontinuous HaploPatterns!");
	}
	AlleleSequence::concatenate(hp);
	m_end = m_start + m_length;
	checkFrequencyWithExtension(m_end-hp.m_length, hp.m_length);
	return *this;
}

HaploPattern &HaploPattern::concatenate(const AlleleSequence &as)
{
	AlleleSequence::concatenate(as);
	m_end = m_start + m_length;
	checkFrequencyWithExtension(m_end-as.length(), as.length());
	return *this;
}

HaploPattern &HaploPattern::concatenate(const Allele &a)
{
	AlleleSequence::concatenate(a);
	m_end = m_start + m_length;
	checkFrequencyWithExtension(m_end-1);
	return *this;
}

HaploPattern &HaploPattern::concatenate(const HaploPattern &hp1, const HaploPattern &hp2)
{
	if (hp1.m_end != hp2.m_start) {
		Logger::warning("[HaploPattern::concatenate] Attempt to concatenate incontinuous HaploPatterns!");
	}
	AlleleSequence::concatenate(hp1, hp2);
	m_start = hp1.m_start;
	m_end = m_start + m_length;
	copyMatchGenotype(hp1);
	checkFrequencyWithExtension(m_end-hp2.m_length, hp2.m_length);
	return *this;
}

HaploPattern &HaploPattern::concatenate(const HaploPattern &hp, const AlleleSequence &as)
{
	AlleleSequence::concatenate(hp, as);
	m_start = hp.m_start;
	m_end = m_start + m_length;
	copyMatchGenotype(hp);
	checkFrequencyWithExtension(m_end-as.length(), as.length());
	return *this;
}

HaploPattern &HaploPattern::concatenate(const AlleleSequence &as, const HaploPattern &hp)
{
	AlleleSequence::concatenate(as, hp);
	m_end = hp.m_end;
	m_start = m_end - m_length;
	copyMatchGenotype(hp);
	checkFrequencyWithExtension(m_start, as.length());
	return *this;
}

HaploPattern &HaploPattern::concatenate(const HaploPattern &hp, const Allele &a)
{
	AlleleSequence::concatenate(hp, a);
	m_start = hp.m_start;
	m_end = m_start + m_length;
	copyMatchGenotype(hp);
	checkFrequencyWithExtension(m_end-1);
	return *this;
}

HaploPattern &HaploPattern::concatenate(const Allele &a, const HaploPattern &hp)
{
	AlleleSequence::concatenate(a, hp);
	m_end = hp.m_end;
	m_start = m_end - m_length;
	copyMatchGenotype(hp);
	checkFrequencyWithExtension(m_start);
	return *this;
}
