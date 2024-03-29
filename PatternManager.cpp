
#include <algorithm>

#include "PatternManager.h"
#include "Allele.h"
#include "Genotype.h"
#include "HaploPattern.h"
#include "HaploBuilder.h"


PatternManager::~PatternManager()
{
	DeleteAll_Clear()(m_patterns);
	DeleteAll_Clear()(m_candidates);
}

HaploPattern *PatternManager::getSingleAllelePattern(int end, int index) const
{
	return pattern_tree()->getSingleAllelePattern(end, index);
}

HaploPattern *PatternManager::getSingleAllelePattern(int end, Allele a) const
{
	return getSingleAllelePattern(end, m_builder.genos()->getAlleleIndex(end-1, a));
}

void PatternManager::findPatternByFreq(double min_freq, int min_len, int max_len)
{
	int geno_len = m_builder.genotype_len();
	max_len = max_len <= 0 ? geno_len : max_len;
	min_len = max(min_len, 1);
	max_len = max(max_len, min_len);
	m_min_len.resize(geno_len, min_len);
	m_max_len.resize(geno_len, max_len);
	DeleteAll_Clear()(m_patterns);
	generateCandidates();
	m_min_freq = min_freq;
	searchPattern();
	DeleteAll_Clear()(m_candidates);
	Logger::verbose("Found haplotype patterns: %d     ", m_patterns.size());
	initialize();
}

void PatternManager::findPatternByNum(int max_num, int min_len, int max_len)
{
	int last_size;
	int geno_len = m_builder.genotype_len();
	max_len = max_len <= 0 ? geno_len : max_len;
	min_len = max(min_len, 1);
	max_len = max(max_len, min_len);
	m_min_len.resize(geno_len, min_len);
	m_max_len.resize(geno_len, max_len);
	DeleteAll_Clear()(m_patterns);
	generateCandidates();
	m_min_freq = 1.0;
	searchPattern(true);
	max_num = max(max_num, m_patterns.size());
	while (m_patterns.size() < max_num && m_min_freq > 1e-38) {
		last_size = m_patterns.size();
		m_min_freq *= 0.9;
		searchPattern(true);
	};
	if (m_patterns.size() > max_num) {
		vector<HaploPattern*>::iterator i_hp = m_patterns.begin() + last_size;
		sort(i_hp, m_patterns.end(), HaploPattern::greater_frequency());
		i_hp = m_patterns.begin() + max_num;
		for_each(i_hp, m_patterns.end(), DeletePtr());
		m_patterns.resize(max_num);
	}
	DeleteAll_Clear()(m_candidates);
	Logger::verbose("Found haplotype patterns: %d     ", m_patterns.size());
	initialize();
}

void PatternManager::findPatternBlock(int len)
{
	int geno_len = m_builder.genotype_len();
	len = max(1, len);
	m_min_len.resize(geno_len, len);
	m_max_len.resize(geno_len, len);
	DeleteAll_Clear()(m_patterns);
	generateCandidates();
	m_min_freq = -1.0;
	searchPattern();
	DeleteAll_Clear()(m_candidates);
	Logger::verbose("Found haplotype patterns: %d     ", m_patterns.size());
	initialize();
}

void PatternManager::generateCandidates()
{
	int geno_len = m_builder.genotype_len();
	DeleteAll_Clear()(m_candidates);
	for (int i=0; i<geno_len; ++i) {
		PatternCandidate *pc = new PatternCandidate(m_builder.genos(), i);
		m_candidates.push_back(pc);
	}
}

void PatternManager::searchPattern(bool reserve_candidates)
{
	int geno_len = m_builder.genotype_len();
	const GenoData *genos = m_builder.genos();
	PatternCandidate *pc_new = new PatternCandidate(m_builder.genos());
	vector<PatternCandidate*> new_candidates;
	while (!m_candidates.empty()) {
		PatternCandidate *pc = m_candidates.back();
		m_candidates.pop_back();
		const HaploPattern *hp = pc->pattern;
		if (hp->frequency() >= m_min_freq || hp->length() < m_min_len[hp->start()]) {
			if (hp->end() < geno_len && hp->length() < m_max_len[hp->start()]) {
				for (int i=0; i<genos->allele_num(hp->end()); ++i) {
					if (genos->allele_frequency(hp->end(), i) > 0) {
						HaploPattern *hp_new = pc_new->pattern;
						hp_new->assign(*hp, genos->allele_symbol(hp->end(), i));
						checkFrequencyWithExtension(hp_new, pc_new->state, pc->state, hp->end());
						hp_new->setPrefixFreq(hp->frequency());
						if (hp_new->prefix_freq() > 0) {
							hp_new->setTransitionProb(hp_new->frequency() / hp_new->prefix_freq());
						}
						else {
							hp_new->setTransitionProb(hp_new->frequency());
						}
						m_candidates.push_back(pc_new);
						pc_new = new PatternCandidate(genos);
					}
				}
			}
		}
		if (hp->frequency() >= m_min_freq || hp->length() <= m_min_len[hp->start()]) {
			if (hp->length() > 0 && hp->length() >= m_min_len[hp->start()]) {
				m_patterns.push_back(pc->release());
			}
			delete pc;
		}
		else if (reserve_candidates) {
			new_candidates.push_back(pc);
		}
		else {
			delete pc;
		}
	}
	m_candidates.swap(new_candidates);
}

void PatternManager::checkFrequency(HaploPattern *hp, MatchingState &ms) const
{
	ms.clear();
	if (hp->length() == 0) {
		hp->setFrequency(1.0);
	}
	else
	{
		double total_freq = 0;
		if (m_builder.samples()->empty()) {
			const GenoData &genos = *m_builder.genos();
			int geno_num = genos.genotype_num();
			for (int i=0; i<geno_num; ++i) {
				const Genotype &g = genos[i];
				double freq = 0;
				if (g.isPhased()) {
					if (hp->isMatch(g(0))) {
						freq += 1.0;
						total_freq += 0.5;
					}
					if (hp->isMatch(g(1))) {
						freq += 2.0;
						total_freq += 0.5;
					}
					if (freq > 0.5) ms.push_back(make_pair(i, freq));
				}
				else if (hp->isMatch(g)) {
					freq = getMatchingFrequency(g, &(*hp)[0], hp->start(), hp->length());
					total_freq += freq;
					ms.push_back(make_pair(i, freq));
				}
			}
			hp->setFrequency(total_freq / geno_num);
		}
		else {
			const HaploData &haplos = *m_builder.samples();
			int haplo_num = haplos.haplotype_num();
			for (int i=0; i<haplo_num; ++i) {
				const Haplotype &h = haplos[i];
				if (hp->isMatch(h)) {
					total_freq += h.weight();
					ms.push_back(make_pair(i, 0));
				}
			}
			hp->setFrequency(total_freq / haplos.total_weight());
		}
	}
}

void PatternManager::checkFrequencyWithExtension(HaploPattern *hp, MatchingState &ms, const MatchingState &oms, int start, int len) const
{
	if (hp->length() == 0 || oms.empty()) {
		checkFrequency(hp, ms);
	}
	else
	{
		ms.clear();
		double total_freq = 0;
		MatchingState::const_iterator i_ms = oms.begin();
		if (m_builder.samples()->empty()) {
			const GenoData &genos = *m_builder.genos();
			while (i_ms != oms.end()) {
				const Genotype &g = genos[i_ms->first];
				double freq = i_ms->second;
				if (g.isPhased()) {
					if (freq > 2.5) {
						if (hp->isMatch(g(0), start, len)) {
							total_freq += 0.5;
						}
						else {
							freq -= 1.0;
						}
						if (hp->isMatch(g(1), start, len)) {
							total_freq += 0.5;
						}
						else {
							freq -= 2.0;
						}
					}
					else if (freq > 1.5) {
						if (hp->isMatch(g(1), start, len)) {
							total_freq += 0.5;
						}
						else {
							freq -= 2.0;
						}
					}
					else if (freq > 0.5) {
						if (hp->isMatch(g(0), start, len)) {
							total_freq += 0.5;
						}
						else {
							freq -= 1.0;
						}
					}
					if (freq > 0.5) ms.push_back(make_pair(i_ms->first, freq));
				}
				else if (hp->isMatch(g, start, len)) {
					freq *= getMatchingFrequency(g, &(*hp)[start-hp->start()], start, len);
					total_freq += freq;
					ms.push_back(make_pair(i_ms->first, freq));
				}
				++i_ms;
			}
			hp->setFrequency(total_freq / m_builder.genotype_num());
		}
		else {
			const HaploData &haplos = *m_builder.samples();
			while (i_ms != oms.end()) {
				const Haplotype &h = haplos[i_ms->first];
				if (hp->isMatch(h, start, len)) {
					total_freq += h.weight();
					ms.push_back(make_pair(i_ms->first, 0));
				}
				++i_ms;
			}
			hp->setFrequency(total_freq / haplos.total_weight());
		}
	}
}

double PatternManager::getMatchingFrequency(const Genotype &g, const Allele *pa, int start, int len) const
{
	double total_freq, freq;
	int i, j;
	Allele b;
	total_freq = 1.0;
	for (i=0; i<len; ++i) {
		if (!pa[i].isMissing()) {			// pa[i] is not missing
			freq = 0;
			for (j=0; j<2; ++j) {
				b = g(j)[start+i];
				if (b.isMissing()) {		// b is missing
					freq += m_builder.genos()->allele_frequency(start+i, pa[i]);
// 					freq += m_builder.genos()->allele_frequency(start+i, pa[i]) > 0 ?
// 						(1.0/m_builder.genos()->allele_num(start+i)) : 0;
				}
				else if (b == pa[i]) {
					freq += 1.0;
				}
			}
			total_freq *= 0.5 * freq;
		}
	}
	return total_freq;
}

void PatternManager::initialize()
{
	int i, j, n;
	AlleleSequence temp;
	m_pattern_tree.reset(new BackwardPatternTree(*m_builder.genos()));
	m_head_list.clear();
	n = m_patterns.size();
	for (i=0; i<n; ++i) {
		HaploPattern *hp = m_patterns[i];
		hp->setID(i);
		m_pattern_tree->addPattern(hp);
		if (hp->start() == 0 && hp->length() == head_len()) {
			m_head_list.push_back(hp);
		}
	}
	for (i=0; i<n; ++i) {
		HaploPattern *hp = m_patterns[i];
		if (hp->end() < m_builder.genotype_len()) {
			temp.assign(*hp, Allele());			// append empty allele
			for (j=0; j<m_builder.genos()->allele_num(hp->end()); ++j) {
				temp[temp.length()-1] = m_builder.genos()->allele_symbol(hp->end(), j);
				hp->setSuccessor(j, m_pattern_tree->findLongestMatchPattern(hp->end()+1, &temp, hp->start()));
			}
		}
	}
}

void PatternManager::adjustFrequency()
{
	int geno_len = m_builder.genotype_len();
	int i, j, n;
	vector<vector<double> > total_freq;
	total_freq.resize(geno_len);
	for (i=0; i<geno_len; ++i) {
		total_freq[i].assign(geno_len-i+1, 0);
		n = min(m_min_len[i], total_freq[i].size());
		for (j=0; j<n; ++j) {
			total_freq[i][j] = 1.0;
		}
	}
	n = m_patterns.size();
	for (i=0; i<n; ++i) {
		HaploPattern *hp = m_patterns[i];
		total_freq[hp->start()][hp->length()] += hp->frequency();
	}
	for (i=0; i<n; ++i) {
		HaploPattern *hp = m_patterns[i];
		if (hp->frequency() > 1.0 || hp->transition_prob() > 1.0 || hp->prefix_freq() > 1.0 ||
			hp->frequency() < 0 || hp->transition_prob() < 0 || hp->prefix_freq() < 0) {
			Logger::error("================ %d, %d, %f, %f, %f", hp->start(), hp->length(), hp->frequency(), hp->prefix_freq(), hp->transition_prob());
		}
	}
}

void PatternManager::estimateFrequency()
{
	int i, n;
	vector<HaploPattern*> patterns;
	n = m_patterns.size();
	for (i=0; i<n; ++i) {
		patterns.push_back(new HaploPattern(*m_patterns[i]));
	}
	m_builder.estimateFrequency(patterns);
	for (i=0; i<n; ++i) {
		m_patterns[i]->setFrequency(patterns[i]->frequency());
		m_patterns[i]->setPrefixFreq(patterns[i]->prefix_freq());
		m_patterns[i]->setTransitionProb(patterns[i]->transition_prob());
	}
	DeleteAll_Clear()(patterns);
}

void PatternManager::estimatePatterns()
{
	int geno_len = m_builder.genotype_len();
	int i, j, n;
	vector<HaploPattern*> patterns, seeds, candidates;
	if (m_min_freq < 0) {
		estimateFrequency();
	}
	else {
		n = m_patterns.size();
		for (i=0; i<n; ++i) {
			HaploPattern *hp = m_patterns[i];
			patterns.push_back(new HaploPattern(*hp));
			if (hp->end() < geno_len && hp->length() < m_max_len[hp->start()]) {
				int an = m_builder.genos()->allele_num(hp->end());
				for (j=0; j<an; ++j) {
					const HaploPattern *succ = hp->successors(j);
					if (!succ || succ->start() != hp->start()) {
						HaploPattern *hp_new = new HaploPattern(*m_builder.genos());
						hp_new->assign(*hp, m_builder.genos()->allele_symbol(hp->end(), j));
						patterns.push_back(hp_new);
						seeds.push_back(hp_new);
					}
				}
			}
		}
		while (!patterns.empty()) {
			m_builder.estimateFrequency(patterns);
			candidates.insert(candidates.end(), patterns.begin(), patterns.end());
			patterns.clear();
			extendPatterns(patterns, seeds);
		}
		DeleteAll_Clear()(m_patterns);
		n = candidates.size();
		for (i=0; i<n; ++i) {
			HaploPattern *hp = candidates[i];
			if (hp->frequency() >= m_min_freq || hp->length() <= m_min_len[hp->start()]) {
				m_patterns.push_back(hp);
  			}
  			else {
  				delete hp;
  			}
		}
		Logger::verbose("Adjust haplotype patterns: %d     ", m_patterns.size());
		initialize();
	}
}

void PatternManager::extendPatterns(vector<HaploPattern*> &patterns, vector<HaploPattern*> &seeds)
{
	int geno_len = m_builder.genotype_len();
	int i, j, n;
	vector<HaploPattern*> new_seeds;
	int level = 4;
	while (level-- > 0) {
		n = seeds.size();
		for (i=0; i<n; ++i) {
			HaploPattern *hp = seeds[i];
			if (hp->end() < geno_len && hp->length() < m_max_len[hp->start()]) {
				if (hp->frequency() >= m_min_freq) {
					int an = m_builder.genos()->allele_num(hp->end());
					for (j=0; j<an; ++j) {
						HaploPattern *hp_new = new HaploPattern(*m_builder.genos());
						hp_new->assign(*hp, m_builder.genos()->allele_symbol(hp->end(), j));
						hp_new->setFrequency(hp->frequency());
						patterns.push_back(hp_new);
						new_seeds.push_back(hp_new);
					}
				}
			}
		}
		seeds.clear();
		seeds.swap(new_seeds);
	}
}
