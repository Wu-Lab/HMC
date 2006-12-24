
#include <string.h>

#include "Genotype.h"


////////////////////////////////
//
// class Haplotype

Genotype::Genotype()
{
	m_id[0] = 0;
	m_length = 0;
	m_heterozygous_num = 0;
	m_weight = 1.0;
	m_is_phased = false;
}

Genotype::Genotype(int len)
{
	m_haplotypes[0].setLength(len);
	m_haplotypes[1].setLength(len);
	m_id[0] = 0;
	m_length = len;
	m_heterozygous_num = 0;
	m_weight = 1.0;
	m_is_phased = false;
}

Genotype::Genotype(const Haplotype &h1, const Haplotype &h2)
{
	if (h1.length() == h2.length()) {
		m_haplotypes[0] = h1;
		m_haplotypes[1] = h2;
		m_id[0] = 0;
		m_length = h1.length();
		checkGenotype();
		m_weight = 1.0;
		m_is_phased = false;
	}
	else {
		Logger::error("Attempt to combine two inconsistent haplotypes!");
		exit(1);
	}
}

Genotype::Genotype(const Genotype &g)
{
	m_haplotypes[0] = g.m_haplotypes[0];
	m_haplotypes[1] = g.m_haplotypes[1];
	strcpy(m_id, g.m_id);
	m_length = g.m_length;
	m_heterozygous_num = g.m_heterozygous_num;
	m_weight = g.m_weight;
	m_is_phased = g.m_is_phased;
}

void Genotype::setHaplotypes(Haplotype &h1, Haplotype &h2)
{
	m_haplotypes[0] = h1;
	m_haplotypes[1] = h2;
	m_length = h1.length();
	m_weight = 1.0;
	m_is_phased = false;
	checkGenotype();
}

void Genotype::checkGenotype()
{
	int i;
	m_heterozygous_num = 0;
	m_missing_num = 0;
	m_missing_allele_num = 0;
	for (i=0; i<m_length; i++) {
		if (isHeterozygous(i)) m_heterozygous_num++;
		if (hasMissing(i)) m_missing_num++;
		if (m_haplotypes[0].isMissing(i)) m_missing_allele_num++;
		if (m_haplotypes[1].isMissing(i)) m_missing_allele_num++;
	}
}

void Genotype::randomizePhase()
{
	int i;
	for (i=0; i<m_length; i++) {
		if (2.0 * rand() < RAND_MAX ) Utils::swap(m_haplotypes[0].allele(i), m_haplotypes[1].allele(i));
	}
	m_is_phased = false;
}

int Genotype::getHeterozygousNum(int start, int len) const
{
	int i;
	int num = 0;
	for (i=start; i<start+len; i++) {
		if (isHeterozygous(i)) {
			num++;
		}
	}
	return num;
}

Genotype &Genotype::assign(const Genotype &g)
{
	m_haplotypes[0] = g.m_haplotypes[0];
	m_haplotypes[1] = g.m_haplotypes[1];
	strcpy(m_id, g.m_id);
	m_length = g.m_length;
	m_heterozygous_num = g.m_heterozygous_num;
	m_weight = g.m_weight;
	m_is_phased = g.m_is_phased;
	return *this;
}

Genotype &Genotype::concatenate(const Genotype &g)
{
	m_haplotypes[0] += g.m_haplotypes[0];
	m_haplotypes[1] += g.m_haplotypes[1];
	m_length += g.m_length;
	m_heterozygous_num += g.m_heterozygous_num;
	return *this;
}

Genotype &Genotype::concatenate(const Genotype &g1, const Genotype &g2)
{
	m_haplotypes[0] = g1.m_haplotypes[0] + g2.m_haplotypes[0];
	m_haplotypes[1] = g1.m_haplotypes[1] + g2.m_haplotypes[1];
	strcpy(m_id, g1.m_id);
	m_length = g1.m_length + g2.m_length;
	m_heterozygous_num = g1.m_heterozygous_num + g2.m_heterozygous_num;
	m_weight = g1.m_weight;
	m_is_phased = g1.m_is_phased;
	return *this;
}

bool Genotype::isMatch(const AlleleSequence &as, int start1, int start2, int len) const
{
	int i;
	bool match;
	if (start1+len > m_length || start2+len > as.length()) {
		match = false;
	}
	else {
		match = true;
		for (i=0; i<len; i++) {
			if (!isMatch(as[start2+i], start1+i)) {
				match = false;
				break;
			}
		}
	}
	return match;
}

bool Genotype::isMatch(const Genotype &g, int i, bool reversed) const
{
	bool match;
	if (reversed) {
		if (m_haplotypes[0].isMatch(g.m_haplotypes[1][i], i) && m_haplotypes[1].isMatch(g.m_haplotypes[0][i], i)) {
			match = true;
		}
		else {
			match = false;
		}
	}
	else {
		if (m_haplotypes[0].isMatch(g.m_haplotypes[0][i], i) && m_haplotypes[1].isMatch(g.m_haplotypes[1][i], i)) {
			match = true;
		}
		else {
			match = false;
		}
	}
	return match;
}

bool Genotype::isMatch(const Genotype &g) const
{
	if ((m_haplotypes[0].isMatch(g(0)) && m_haplotypes[1].isMatch(g(1))) ||
		(m_haplotypes[0].isMatch(g(1)) && m_haplotypes[1].isMatch(g(0)))) {
		return true;
	}
	else {
		return false;
	}
}

bool Genotype::isMatchUnphased(const Genotype &g) const
{
	int i;
	bool match;
	match = true;
	for (i=0; i<m_length; i++) {
		if (!isMatch(g, i, true) && !isMatch(g, i, false)) {
			match = false;
			break;
		}
	}
	return match;
}

int Genotype::getDiffNum(const Genotype &g) const
{
	int i, diff1, diff2;
	diff1 = diff2 = 0;
	for (i=0; i<m_length; i++) {
		if (!isMatch(g, i, true)) {
			diff1++;
		}
		if (!isMatch(g, i, false)) {
			diff2++;
		}
	}
	return (diff1 < diff2 ? diff1 : diff2);
}

int Genotype::getDiffNumIgnoreMissing(const Genotype &g) const
{
	int i, diff1, diff2;
	diff1 = diff2 = 0;
	for (i=0; i<m_length; i++) {
		if (!hasMissing(i)) {
			if (!isMatch(g, i, true)) {
				diff1++;
			}
			if (!isMatch(g, i, false)) {
				diff2++;
			}
		}
	}
	return (diff1 < diff2 ? diff1 : diff2);
}

int Genotype::getSwitchDistance(const Genotype &g) const
{
	int switch_distance, start, i;
	bool reversed;
 	if (m_length != g.m_length) {
 		Logger::error("Inconsistent genotype length while calculate switch distance!");
 		exit(1);
 	}
	start = m_length;
	for (i=0; i<m_length; i++) {
		if (isMatch(g, i, true) && isMatch(g, i, false)) {
			continue;
		}
		else {
			start = i;
			break;
		}
	}
	switch_distance = 0;
	if (start < m_length) {
		if (isMatch(g, start, true)) {
			reversed = true;
		}
		else if (isMatch(g, start, false)) {
			reversed = false;
		}
		else {
			Logger::error("Inconsistent genotypes at locus %d!", start);
			exit(1);
		}
		for (i=start+1; i<m_length; i++) {
			if (isMatch(g, i, reversed)) {
				continue;
			}
			else if (isMatch(g, i, !reversed)) {
				reversed = !reversed;
				switch_distance++;
			}
			else {
				Logger::error("Inconsistent genotypes at locus %d!", i);
				exit(1);
			}
		}
	}
	return switch_distance;
}

int Genotype::getSwitchDistanceIgnoreMissing(const Genotype &g) const
{
	int switch_distance, start, i;
	bool reversed;
 	if (m_length != g.m_length) {
 		Logger::error("Inconsistent genotype length while calculate switch distance!");
 		exit(1);
 	}
	start = m_length;
	for (i=0; i<m_length; i++) {
		if (hasMissing(i) || (isMatch(g, i, true) && isMatch(g, i, false))) {
			continue;
		}
		else {
			start = i;
			break;
		}
	}
	switch_distance = 0;
	if (start < m_length) {
		if (isMatch(g, start, true)) {
			reversed = true;
		}
		else if (isMatch(g, start, false)) {
			reversed = false;
		}
		else {
			Logger::error("Inconsistent genotypes at locus %d!", start);
			exit(1);
		}
		for (i=start+1; i<m_length; i++) {
			if (hasMissing(i) || isMatch(g, i, reversed)) {
				continue;
			}
			else if (isMatch(g, i, !reversed)) {
				reversed = !reversed;
				switch_distance++;
			}
			else {
				Logger::error("Inconsistent genotypes at locus %d!", i);
				exit(1);
			}
		}
	}
	return switch_distance;
}
