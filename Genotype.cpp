
#include "Genotype.h"

#include "MemLeak.h"


////////////////////////////////
//
// class Genotype

Genotype::Genotype(const Haplotype &h1, const Haplotype &h2)
: m_heterozygous_num(0),
  m_missing_num(0),
  m_missing_allele_num(0),
  m_prior_probability(0),
  m_posterior_probability(0),
  m_genotype_probability(0),
  m_is_phased(false)
{
	if (h1.length() == h2.length()) {
		m_haplotypes[0] = h1;
		m_haplotypes[1] = h2;
		checkGenotype();
	}
	else {
		Logger::error("Attempt to combine two inconsistent haplotypes!");
		exit(1);
	}
}

void Genotype::setHaplotypes(Haplotype &h1, Haplotype &h2)
{
	m_haplotypes[0] = h1;
	m_haplotypes[1] = h2;
	m_prior_probability = 1.0;
	m_posterior_probability = 1.0;
	m_genotype_probability = 1.0;
	m_is_phased = false;
	checkGenotype();
}

void Genotype::checkGenotype()
{
	int i;
	m_heterozygous_num = 0;
	m_missing_num = 0;
	m_missing_allele_num = 0;
	for (i=0; i<length(); i++) {
		if (isHeterozygous(i)) m_heterozygous_num++;
		if (hasMissing(i)) m_missing_num++;
		if (m_haplotypes[0][i].isMissing()) m_missing_allele_num++;
		if (m_haplotypes[1][i].isMissing()) m_missing_allele_num++;
	}
}

void Genotype::randomizePhase()
{
	int i;
	for (i=0; i<length(); i++) {
		if (2.0 * rand() < RAND_MAX ) swap(m_haplotypes[0][i], m_haplotypes[1][i]);
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

bool Genotype::isMatch(const AlleleSequence &as, int start1, int start2, int len) const
{
	int i;
	bool match;
	if (start1+len > length() || start2+len > as.length()) {
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
		if (m_haplotypes[0][i].isMatch(g.m_haplotypes[1][i]) && m_haplotypes[1][i].isMatch(g.m_haplotypes[0][i])) {
			match = true;
		}
		else {
			match = false;
		}
	}
	else {
		if (m_haplotypes[0][i].isMatch(g.m_haplotypes[0][i]) && m_haplotypes[1][i].isMatch(g.m_haplotypes[1][i])) {
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
	for (i=0; i<length(); i++) {
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
	for (i=0; i<length(); i++) {
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
	for (i=0; i<length(); i++) {
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
 	if (length() != g.length()) {
 		Logger::error("Inconsistent genotype length while calculate switch distance!");
 		exit(1);
 	}
	start = length();
	for (i=0; i<length(); i++) {
		if (isMatch(g, i, true) && isMatch(g, i, false)) {
			continue;
		}
		else {
			start = i;
			break;
		}
	}
	switch_distance = 0;
	if (start < length()) {
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
		for (i=start+1; i<length(); i++) {
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
 	if (length() != g.length()) {
 		Logger::error("Inconsistent genotype length while calculate switch distance!");
 		exit(1);
 	}
	start = length();
	for (i=0; i<length(); i++) {
		if (hasMissing(i) || (isMatch(g, i, true) && isMatch(g, i, false))) {
			continue;
		}
		else {
			start = i;
			break;
		}
	}
	switch_distance = 0;
	if (start < length()) {
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
		for (i=start+1; i<length(); i++) {
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
