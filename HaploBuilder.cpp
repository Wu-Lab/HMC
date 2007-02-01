
#include <deque>

#include "HaploBuilder.h"
#include "HaploPair.h"
#include "HaploData.h"

#include "MemLeak.h"


HaploBuilder::HaploBuilder()
{
}

HaploBuilder::~HaploBuilder()
{
	DeleteAll_Clear()(m_haplo_pattern);
	for_each(m_haplopairs.begin(), m_haplopairs.end(), DeleteAll_Clear());
}

void HaploBuilder::setHaploData(HaploData &hd)
{
	m_haplo_data = &hd;
	m_genotype_len = hd.genotype_len(); 
}

void HaploBuilder::initialize()
{
	int i, j;
	HaploPattern *hp;
	AlleleSequence temp;
	vector<HaploPattern*>::iterator i_hp;
	list<HaploPattern*>::iterator head;
	m_pattern_num = m_haplo_pattern.size();
	m_pattern_tree.reset(new BackwardPatternTree(m_haplo_data));
	m_head_len = m_genotype_len;
	m_head_list.clear();
	i = 0;
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		hp->setID(i++);
		m_pattern_tree->addPattern(hp);
		if (hp->start() == 0) {
			m_head_list.push_back(hp);
			if (hp->length() < m_head_len) m_head_len = hp->length();
		}
	}
	head = m_head_list.begin();
	while (head != m_head_list.end()) {
		if ((*head)->length() > m_head_len) {
			head = m_head_list.erase(head);
		}
		else {
			head++;
		}
	}
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		if (hp->end() < m_genotype_len) {
			temp.assign(*hp, Allele());			// append empty allele
			for (j=0; j<m_haplo_data->allele_num(hp->end()); j++) {
				temp[temp.length()-1] = m_haplo_data->allele_symbol(hp->end(), j);
				hp->setSuccessor(j, m_pattern_tree->findLongestMatchPattern(hp->end()+1, &temp, hp->start()));
			}
		}
	}
	for_each(m_haplopairs.begin(), m_haplopairs.end(), DeleteAll_Clear());
	m_haplopairs.resize(m_genotype_len+1);
	m_best_pair.resize(m_pattern_num);
}

void HaploBuilder::adjust(double min_freq)
{
	int i;
	vector<double> freq;
	HaploPattern *hp;
	vector<HaploPattern*>::iterator i_hp;
	Genotype res;
	vector<HaploPair*> res_list;
	freq.resize(m_pattern_num, 0);
	for (i=0; i<m_haplo_data->genotype_num(); ++i) {
		resolve((*m_haplo_data)[i], res, res_list);
		calcBackwardLikelihood();
		for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
			hp = *i_hp;
			if (hp->isMatch((*m_haplo_data)[i])) {
				freq[hp->id()] += HaploPair::evaluatePattern(hp, m_haplopairs[hp->start()+1]) / (*m_haplo_data)[i].likelihood();
			}
		}
		Logger::status("Adjust haplotype patterns: %d     ", i);
	}
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		hp->setFrequency(freq[hp->id()] / m_haplo_data->genotype_num());
	}
}

void HaploBuilder::resolve(const Genotype &genotype, Genotype &resolution, vector<HaploPair*> &res_list, HaploPattern *target_pattern)
{
	int i, j, k;
	Allele a, b;
	double total_likelihood;
	HaploPair *best_hp;
	vector<HaploPair*>::iterator i_hp;
	for_each(m_haplopairs.begin(), m_haplopairs.end(), DeleteAll_Clear());
	for (i=0; i<m_pattern_num; i++) {
		m_best_pair[i].clear();
	}
	initHeadList(genotype);
	for (i=m_head_len; i<m_genotype_len; i++) {
		if (genotype.isMissing(i)) {
			for (j=0; j<m_haplo_data->allele_num(i); j++) {
				if (m_haplo_data->allele_frequency(i, j) > 0) {
					for (k=j; k<m_haplo_data->allele_num(i); k++) {
						if (m_haplo_data->allele_frequency(i, k) > 0) {
							a = m_haplo_data->allele_symbol(i, j);
							b = m_haplo_data->allele_symbol(i, k);
							extendAll(i, a, b);
						}
					}
				}
			}
		}
		else if (genotype(0)[i].isMissing()) {
			for (j=0; j<m_haplo_data->allele_num(i); j++) {
				if (m_haplo_data->allele_frequency(i, j) > 0) {
					a = m_haplo_data->allele_symbol(i, j);
					extendAll(i, a, genotype(1)[i]);
				}
			}
		}
		else if (genotype(1)[i].isMissing()) {
			for (j=0; j<m_haplo_data->allele_num(i); j++) {
				if (m_haplo_data->allele_frequency(i, j) > 0) {
					a = m_haplo_data->allele_symbol(i, j);
					extendAll(i, a, genotype(0)[i]);
				}
			}
		}
		else {
			extendAll(i, genotype(0)[i], genotype(1)[i]);
		}
		if (m_haplopairs[i+1].size() <= 0) {
			break;
		}
	}
	if (m_haplopairs[m_genotype_len].size() > 0) {
		total_likelihood = 0;
		res_list.clear();
		for (i_hp = m_haplopairs[m_genotype_len].begin(); i_hp != m_haplopairs[m_genotype_len].end(); i_hp++) {
			total_likelihood += (*i_hp)->m_forward_likelihood;
			res_list.push_back(*i_hp);
		}
		sort(res_list.begin(), res_list.end(), HaploPair::greater_likelihood());
		best_hp = res_list.front();
		resolution = best_hp->getGenotype();
		resolution.setLikelihood(total_likelihood);
		resolution.setWeight(best_hp->best_likelihood() / total_likelihood);
	}
	else {
		res_list.clear();
		resolution = genotype;
		resolution.setLikelihood(0);
		resolution.setWeight(0);
	}
}

double HaploBuilder::getLikelihood(const Haplotype &haplotype)
{
	int i;
	HaploPattern *hp;
	double likelihood = 1.0;
	for (i=m_head_len; i<=m_genotype_len; i++) {
		hp = m_pattern_tree->findLongestMatchPattern(i, &haplotype);
		if (hp) {
			likelihood *= hp->transition_prob();
		}
		else {
			likelihood = 0;
			break;
		}
	}
	return likelihood;
}

double HaploBuilder::getLikelihood(const Genotype &genotype)
{
	return getLikelihood(genotype(0)) * getLikelihood(genotype(1));
}

void HaploBuilder::initHeadList(const Genotype &genotype)
{
	int j, k;
	list<HaploPattern*>::iterator head;
	deque<AlleleSequence*> as_list;
	deque<AlleleSequence*>::iterator i_as;
	AlleleSequence *as, *new_as;
	HaploPattern *hp;
	for (head = m_head_list.begin(); head != m_head_list.end(); head++) {
		if ((*head)->isMatch(genotype)) {
			as_list.push_front(new AlleleSequence);
			for (j=0; j<m_head_len; j++) {
				i_as = as_list.begin();
				if (genotype.isMissing(j) || (genotype.hasMissing(j) && genotype.hasAllele(j, (**head)[j]))) {
					while (i_as != as_list.end()) {
						as = *i_as;
						for (k=0; k<m_haplo_data->allele_num(j); k++) {
							if (m_haplo_data->allele_frequency(j, k) > 0) {
								new_as = new AlleleSequence;
								new_as->assign(*as, m_haplo_data->allele_symbol(j, k));
								as_list.push_front(new_as);
							}
						}
						delete *i_as;
						i_as = as_list.erase(i_as);
					}
				}
				else if (genotype.isHeterozygous(j)) {
					while (i_as != as_list.end()) {
						as = *i_as;
						new_as = new AlleleSequence;
						if ((**head)[j] == genotype(0)[j]) {
							new_as->assign(*as, genotype(1)[j]);
						}
						else {
							new_as->assign(*as, genotype(0)[j]);
						}
						as_list.push_front(new_as);
						delete *i_as;
						i_as = as_list.erase(i_as);
					}
				}
				else {
					while (i_as != as_list.end()) {
						as = *i_as;
						new_as = new AlleleSequence;
						new_as->assign(*as, genotype(0)[j]);
						as_list.push_front(new_as);
						delete *i_as;
						i_as = as_list.erase(i_as);
					}
				}

			}
			i_as = as_list.begin();
			while (i_as != as_list.end()) {
				hp = m_pattern_tree->findLongestMatchPattern(m_head_len, *i_as);
				if (hp && hp->start() == 0) {
					HaploPair *new_hp = new HaploPair(*head, hp);
					m_haplopairs[m_head_len].push_back(new_hp);
				 	m_best_pair[new_hp->id_a()].insert(make_pair(new_hp->id_b(), m_haplopairs[m_head_len].size()));
				}
				else {
					Logger::error("Can not find matching pattern!");
					exit(1);
				}
				delete *i_as;
				i_as = as_list.erase(i_as);
			}
		}
	}
}

void HaploBuilder::extendAll(int i, Allele a1, Allele a2)
{
	vector<HaploPair*>::iterator i_hp;
	for (i_hp = m_haplopairs[i].begin(); i_hp != m_haplopairs[i].end(); ++i_hp) {
		extend(*i_hp, a1, a2);
		if (a1 != a2) extend(*i_hp, a2, a1);
	}
}

void HaploBuilder::extend(HaploPair *hp, Allele a1, Allele a2)
{
	const HaploPattern *hpa, *hpb;
	hpa = hp->successor_a(a1);
	hpb = hp->successor_b(a2);
	if (hpa && hpb) {
		map<int, int>::iterator i = m_best_pair[hpa->id()].lower_bound(hpb->id());
		if ((*i).first != hpb->id()) {
			m_haplopairs[hp->end()+1].push_back(new HaploPair(hp, hpa, hpb));
		 	m_best_pair[hpa->id()].insert(i, make_pair(hpb->id(), m_haplopairs[hp->end()+1].size()));
		}
		else {
			m_haplopairs[hp->end()+1][(*i).second-1]->add(hp, hpa, hpb);
		}
	}
}

void HaploBuilder::calcBackwardLikelihood()
{
	int i;
	vector<HaploPair*>::iterator i_hp;

	for (i=m_genotype_len-1; i>=m_head_len; --i) {
		for (i_hp=m_haplopairs[i].begin(); i_hp!=m_haplopairs[i].end(); ++i_hp) {
			(*i_hp)->calcBackwardLikelihood();
		}
	}
}

void HaploBuilder::findHaploPatternByFreq(double min_freq, int min_len, int max_len)
{
	vector<int> min_len_vector, max_len_vector;
	vector<HaploPattern*> candidates;
	list<HaploPattern*>::iterator i_hp;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector.resize(m_genotype_len, min_len);
	max_len_vector.resize(m_genotype_len, max_len);
	DeleteAll_Clear()(m_haplo_pattern);
	generateHaploPatternCandidate(candidates);
	searchHaploPattern(candidates, min_freq, -1, min_len_vector, max_len_vector);
	DeleteAll_Clear()(candidates);
	Logger::verbose("Found haplotype patterns: %d     ", m_haplo_pattern.size());
}

void HaploBuilder::findHaploPatternByNum(int max_num, int min_len, int max_len)
{
	int last_size;
	vector<int> min_len_vector, max_len_vector;
	double min_freq = 10;
	vector<HaploPattern*> candidates;
	vector<HaploPattern*>::iterator i_hp;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector.resize(m_genotype_len, min_len);
	max_len_vector.resize(m_genotype_len, max_len);
	DeleteAll_Clear()(m_haplo_pattern);
	generateHaploPatternCandidate(candidates);
	do {
		last_size = m_haplo_pattern.size();
		min_freq *= 0.9;
		searchHaploPattern(candidates, min_freq, max_num, min_len_vector, max_len_vector);
	} while (m_haplo_pattern.size() < max_num && min_freq > 1e-38);
	if (m_haplo_pattern.size() > max_num) {
		i_hp = m_haplo_pattern.begin() + last_size;
		sort(i_hp, m_haplo_pattern.end(), HaploPattern::greater_frequency());
		i_hp = m_haplo_pattern.begin() + max_num;
		for_each(i_hp, m_haplo_pattern.end(), DeletePtr());
		m_haplo_pattern.resize(max_num);
	}
	DeleteAll_Clear()(candidates);
	Logger::verbose("Found haplotype patterns: %d     ", m_haplo_pattern.size());
}

void HaploBuilder::findHaploPatternBlockByFreq(double min_freq, int min_len, int max_len)
{
	findHaploPatternByFreq(min_freq);
	searchHaploPatternBlock(min_len, max_len);
	Logger::verbose("Found haplotype patterns by block: %d     ", m_haplo_pattern.size());
}

void HaploBuilder::findHaploPatternBlockByNum(int max_num, int min_len, int max_len)
{
	findHaploPatternByNum(max_num);
	searchHaploPatternBlock(min_len, max_len);
	Logger::verbose("Found haplotype patterns by block: %d     ", m_haplo_pattern.size());
}

void HaploBuilder::generateHaploPatternCandidate(vector<HaploPattern*> &candidates)
{
	int i;
	HaploPattern *hp;
	DeleteAll_Clear()(candidates);
	for (i=0; i<m_genotype_len; i++) {
		hp = new HaploPattern(*m_haplo_data, i);
		candidates.push_back(hp);
	}
}

void HaploBuilder::searchHaploPattern(vector<HaploPattern*> &candidates, double min_freq, int min_num, vector<int> min_len, vector<int> max_len)
{
	int i;
	vector<HaploPattern*> new_candidates;
	HaploPattern *hp, *hp_candidate;
	while (!candidates.empty()) {
		hp = candidates.back();
		candidates.pop_back();
		if (hp->frequency() >= min_freq || hp->length() <= 2) {
			if (hp->end() < m_genotype_len && hp->length() < max_len[hp->end()]) {
				for (i=0; i<m_haplo_data->allele_num(hp->end()); i++) {
					if (m_haplo_data->allele_frequency(hp->end(), i) > 0) {
						hp_candidate = new HaploPattern(*m_haplo_data, hp->end());
						hp_candidate->assign(*hp, m_haplo_data->allele_symbol(hp->end(), i));
						if (hp_candidate->frequency() >= min_freq || hp->length() <= 2)
						{
							candidates.push_back(hp_candidate);
						}
						else if (hp_candidate->frequency() > 0)
						{
							new_candidates.push_back(hp_candidate);
						}
						else {
							delete hp_candidate;
						}
					}
				}
			}
			if (hp->length() > 0 && (hp->length() >= min_len[hp->end()-1] || (min_freq >= 0 && hp->length() <= 2))) {
				hp->releaseMatchGenotype();
				m_haplo_pattern.push_back(hp);
			}
			else {
				delete hp;
			}
		}
		else if (min_num < 0) {
			delete hp;
		}
		else {
			new_candidates.push_back(hp);
		}
	}
	candidates.swap(new_candidates);
}

void HaploBuilder::searchHaploPatternBlock(int min_len, int max_len)
{
	int i;
	vector<int> min_len_vector, max_len_vector;
	HaploPattern *hp;
	vector<HaploPattern*> candidates;
	vector<HaploPattern*>::iterator i_hp;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector.resize(m_genotype_len, min_len);
	max_len_vector.resize(m_genotype_len, min_len+1);
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		if (hp->length() >= min_len && hp->length() <= max_len 
			&& hp->length() > max_len_vector[hp->end()-1]) max_len_vector[hp->end()-1] = hp->length();
		delete hp;
	}
	m_haplo_pattern.clear();
	for (i=0; i<m_genotype_len; i++) {
		min_len_vector[i] = max_len_vector[i]-1;
		Logger::debug("#%05d %5d %5d", i, min_len_vector[i], max_len_vector[i]);
	}
	generateHaploPatternCandidate(candidates);
	searchHaploPattern(candidates, 0, -1, min_len_vector, max_len_vector);
	DeleteAll_Clear()(candidates);
}
