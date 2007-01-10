
#include "HaploBuilder.h"

#include <boost/pool/object_pool.hpp>


boost::pool<> HaploPair_pool(sizeof(HaploPair));

void del_HaploPair(HaploPair *p)
{
	p->~HaploPair();
	HaploPair_pool.free(p);
}


PatternTree::PatternTree(const HaploData *haplo)
{
	int i;
	m_haplo_data = haplo;
	m_genotype_len = haplo->genotype_len();
	m_trees = new Tree<HaploPattern, double> [m_genotype_len+1];
	for (i=0; i<=m_genotype_len; i++) {
		m_trees[i].setChildrenNum(Constant::max_allele_num());
	}
}

PatternTree::~PatternTree()
{
	int i;
	for (i=0; i<=m_genotype_len; i++) {
		m_trees[i].release();
	}
	delete[] m_trees;
}

void PatternTree::addPattern(int end, HaploPattern *hp)
{
	if (hp->m_end != end) {
		Logger::error("Attempt to add HaploPattern to incorrect PatternTree!");
		exit(1);
	}
	else {
		addPattern(&m_trees[end].root(), hp, hp->m_length);
	}
}

void PatternTree::addPattern(PatternNode *node, HaploPattern *hp, int len)
{
	int i, n;
	i = hp->getAlleleIndex(len-1);
	if (i < 0) {										// allele is missing
		n = m_haplo_data->allele_num(hp->getGlobalLocus(len-1));
		for (i=0; i<n; i++) {
			if (node->getChildNode(i) == NULL) {
				node->addChild(i);
			}
		}
		if (len == 1) {
			for (i=0; i<n; i++) {
				node->getChildNode(i)->setObject(hp, hp->m_frequency);
			}
		}
		else {
			for (i=0; i<n; i++) {
				addPattern(node->getChildNode(i), hp, len-1);
			}
		}
	}
	else {
		if (node->getChildNode(i) == NULL) {
			node->addChild(i);
		}
		if (len == 1) {
			node->getChildNode(i)->setObject(hp, hp->m_frequency);
		}
		else {
			addPattern(node->getChildNode(i), hp, len-1);
		}
	}
}

HaploPattern *PatternTree::findLongestMatchPattern(int end, const HaploPattern *hp, int len)
{
	if (hp->m_end < end || hp->m_start >= end) {
		return NULL;
	}
	else {
		if (len <= 0) len = end - hp->m_start;
		return findLongestMatchPattern(&m_trees[end].root(), hp, end, len, end-hp->m_start-len);
	}
}

HaploPattern *PatternTree::findLongestMatchPattern(int end, const AlleleSequence *as, int len, int offset)
{
	if (len <= 0) len = as->length();
	return findLongestMatchPattern(&m_trees[end].root(), as, end, len, offset);
}

HaploPattern *PatternTree::findLikelyMatchPattern(int end, const HaploPattern *hp, int len)
{
	if (hp->m_end < end || hp->m_start >= end) {
		return NULL;
	}
	else {
		if (len <= 0) len = end - hp->m_start;
		return findLikelyMatchPattern(&m_trees[end].root(), hp, end, len, end-hp->m_start-len);
	}
}

HaploPattern *PatternTree::findLikelyMatchPattern(int end, const AlleleSequence *as, int len, int offset)
{
	if (len <= 0) len = as->length();
	return findLikelyMatchPattern(&m_trees[end].root(), as, end, len, offset);
}

HaploPattern *PatternTree::findLongestMatchPattern(PatternNode *node, const AlleleSequence *as, int end, int len, int offset)
{
	int i, n;
	HaploPattern *result, *temp;
	result = node->object();
	if (as->isMissing(offset+len-1)) {						// allele is missing
		n = m_haplo_data->allele_num(end-1);
		for (i=0; i<n; i++) {
			if (node->getChildNode(i) != NULL) {			// previous locus is matching
				if (len > 1) {
					temp = findLongestMatchPattern(node->getChildNode(i), as, end-1, len-1, offset);
				}
				else {
					temp = node->getChild(i);
				}
				if (result == NULL || (temp != NULL && temp->m_length > result->m_length)
					|| (temp != NULL && temp->m_length >= result->m_length 
					&& temp->m_transition_prob > result->m_transition_prob)) {
					result = temp;
				}
			}
		}
	}
	else {
		i = m_haplo_data->getAlleleIndex(end-1, (*as)[offset+len-1]);
		if (node->getChildNode(i) != NULL) {				// previous locus is matching
			if (len > 1) {
				temp = findLongestMatchPattern(node->getChildNode(i), as, end-1, len-1, offset);
			}
			else {
				temp = node->getChild(i);
			}
			if (result == NULL || (temp != NULL && temp->m_length > result->m_length)
				|| (temp != NULL && temp->m_length >= result->m_length 
				&& temp->m_transition_prob > result->m_transition_prob)) {
				result = temp;
			}
		}
	}
	return result;
}

HaploPattern *PatternTree::findLikelyMatchPattern(PatternNode *node, const AlleleSequence *as, int end, int len, int offset)
{
	int i, n;
	HaploPattern *result, *temp;
	result = node->object();
	if (as->isMissing(offset+len-1)) {						// allele is missing
		n = m_haplo_data->allele_num(end-1);
		for (i=0; i<n; i++) {
			if (node->getChildNode(i) != NULL) {			// previous locus is matching
				if (len > 1) {
					temp = findLikelyMatchPattern(node->getChildNode(i), as, end-1, len-1, offset);
				}
				else {
					temp = node->getChild(i);
				}
				if (result == NULL || (temp != NULL && temp->m_transition_prob >= result->m_transition_prob)) {
					result = temp;
				}
			}
		}
	}
	else {
		i = m_haplo_data->getAlleleIndex(end-1, (*as)[offset+len-1]);
		if (node->getChildNode(i) != NULL) {				// previous locus is matching
			if (len > 1) {
				temp = findLikelyMatchPattern(node->getChildNode(i), as, end-1, len-1, offset);
			}
			else {
				temp = node->getChild(i);
			}
			if (result == NULL || (temp != NULL && temp->m_transition_prob >= result->m_transition_prob)) {
				result = temp;
			}
		}
	}
	return result;
}

HaploBuilder::HaploBuilder()
{
	m_pattern_tree = NULL;
	m_best_pair = NULL;
	m_target_pattern = NULL;
}

HaploBuilder::~HaploBuilder()
{
	delete m_pattern_tree;
	delete[] m_best_pair;
	DeleteAllObjects(m_haplo_pattern);
}

void HaploBuilder::clear()
{
	delete m_pattern_tree;
	delete[] m_best_pair;
	m_pattern_tree = NULL;
	m_best_pair = NULL;
	m_target_pattern = NULL;
}

void HaploBuilder::setHaploData(HaploData &hd)
{
	m_haplo_data = &hd;
	m_genotype_len = hd.m_genotype_len; 
}

void HaploBuilder::initialize()
{
	int i, j;
	AlleleSequence temp;
	HaploPattern *hp;
	list<HaploPattern*>::iterator i_hp;
	list<HaploPattern*>::iterator head;
	clear();
	m_pattern_num = m_haplo_pattern.size();
	m_pattern_tree = new PatternTree(m_haplo_data);
	m_head_len = m_genotype_len;
	m_head_list.clear();
	i = 0;
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		hp->m_id = i++;
		m_pattern_tree->addPattern(hp->m_end, hp);
		if (hp->m_start == 0) {
			m_head_list.push_back(hp);
			if (hp->m_length < m_head_len) m_head_len = hp->m_length;
		}
	}
	head = m_head_list.begin();
	while (head != m_head_list.end()) {
		if ((*head)->m_length > m_head_len) {
			head = m_head_list.erase(head);
		}
		else {
			head++;
		}
	}
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		hp->m_prefix = m_pattern_tree->findLongestMatchPattern(hp->m_end-1, hp);
		if (hp->m_prefix != NULL) {
			hp->m_transition_prob = hp->m_frequency / hp->m_prefix->m_frequency;
		}
		else {
			hp->m_transition_prob = hp->m_frequency / m_haplo_data->m_genotype_num;
		}
		hp->m_weight = 0;
	}
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		if (hp->m_end < m_genotype_len) {
			temp.concatenate(*hp, -1);			// append empty allele
			for (j=0; j<m_haplo_data->allele_num(hp->m_end); j++) {
				temp.allele(temp.length()-1) = m_haplo_data->allele_symbol(hp->m_end, j);
				hp->m_successors[j] = m_pattern_tree->findLongestMatchPattern(hp->m_end+1, &temp);
			}
		}
	}
	m_best_pair = new vector<pair<int, int> > [m_pattern_num];
}

void HaploBuilder::adjust(double min_freq)
{
	int i;
	double *new_frequency;
	HaploPattern *hp;
	list<HaploPattern*>::iterator i_hp;
	Genotype res;
	vector<HaploPair*> res_list;
	new_frequency = new double [m_pattern_num];
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		new_frequency[hp->id()] = 0;
		for (i=0; i<m_haplo_data->m_genotype_num; i++) {
			if (hp->isMatch(m_haplo_data->genotype(i))) {
				resolve(m_haplo_data->genotype(i), res, res_list, hp);
				new_frequency[hp->id()] += res.likelihood() / m_haplo_data->genotype(i).likelihood();
			}
		}
		Logger::status("Adjust haplotype patterns: %d     ", hp->id());
	}
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		hp->m_frequency = new_frequency[hp->id()];
	}
	delete[] new_frequency;
	i_hp = m_haplo_pattern.begin();
	while (i_hp != m_haplo_pattern.end()) {
		hp = *i_hp;
		if (hp->m_prefix != NULL) {
			if (hp->m_frequency > hp->m_prefix->m_frequency * 1.0001) {
				Logger::error("Mistake in adjusting the pattern frequencies!");
				exit(1);
			}
		}
		if (hp->m_frequency < min_freq && hp->m_length > 1) {
			delete hp;
			i_hp = m_haplo_pattern.erase(i_hp);
		}
		else {
			i_hp++;
		}
	}
	initialize();
}

void HaploBuilder::resolve(const Genotype &genotype, Genotype &resolution, vector<HaploPair*> &res_list, HaploPattern *target_pattern)
{
	int i, j, k;
	Allele a, b;
	double weight, total_likelihood;
	HaploPair *best_hp;
	vector<HaploPair*>::iterator i_hp;
	m_target_pattern = target_pattern;
	m_last_list.clear();
	m_new_list.clear();
	for (i=0; i<m_pattern_num; i++) {
		m_best_pair[i].clear();
	}
	initHeadList(genotype);
	for (i=m_head_len; i<m_genotype_len; i++) {
		if (genotype.isMissing(i)) {
			for (j=0; j<m_haplo_data->allele_num(i); j++) {
				if (m_haplo_data->m_allele_frequency[i][j] > 0) {
					for (k=j; k<m_haplo_data->allele_num(i); k++) {
						if (m_haplo_data->m_allele_frequency[i][k] > 0) {
							a = m_haplo_data->allele_symbol(i, j);
							b = m_haplo_data->allele_symbol(i, k);
							extendAll(a, b);
						}
					}
				}
			}
		}
		else if (genotype(0).isMissing(i)) {
			for (j=0; j<m_haplo_data->allele_num(i); j++) {
				if (m_haplo_data->m_allele_frequency[i][j] > 0) {
					a = m_haplo_data->allele_symbol(i, j);
					extendAll(a, genotype(1)[i]);
				}
			}
		}
		else if (genotype(1).isMissing(i)) {
			for (j=0; j<m_haplo_data->allele_num(i); j++) {
				if (m_haplo_data->m_allele_frequency[i][j] > 0) {
					a = m_haplo_data->allele_symbol(i, j);
					extendAll(a, genotype(0)[i]);
				}
			}
		}
		else {
			extendAll(genotype(0)[i], genotype(1)[i]);
		}
		for_each(m_last_list.begin(), m_last_list.end(), del_HaploPair);
		m_last_list.swap(m_new_list);
		m_new_list.clear();
		if (m_last_list.size() <= 0) {
			break;
		}
	}
	if (m_last_list.size() > 0) {
		weight = m_haplo_data->m_genotype_num * m_haplo_data->m_genotype_num;
		total_likelihood = 0;
		res_list.clear();
		for (i_hp = m_last_list.begin(); i_hp != m_last_list.end(); i_hp++) {
			(*i_hp)->m_likelihood /= weight;
			(*i_hp)->m_total_likelihood /= weight;
			total_likelihood += (*i_hp)->m_total_likelihood;
			res_list.push_back(*i_hp);
		}
		sort(res_list.begin(), res_list.end(), greater_likelihood);
		best_hp = res_list.front();
		resolution = best_hp->getGenotype();
		resolution.setLikelihood(total_likelihood);
		resolution.setWeight(best_hp->likelihood() / total_likelihood);
	}
	else {
		res_list.clear();
		resolution = genotype;
		resolution.setLikelihood(0);
		resolution.setWeight(0);
	}
	for_each(m_last_list.begin(), m_last_list.end(), del_HaploPair);
	m_last_list.clear();
}

double HaploBuilder::getLikelihood(const Haplotype &haplotype)
{
	int i;
	HaploPattern *hp;
	double likelihood = 1.0;
	for (i=m_head_len; i<=m_genotype_len; i++) {
		hp = m_pattern_tree->findLongestMatchPattern(i, &haplotype, i);
		if (hp != NULL) {
			likelihood *= hp->m_transition_prob;
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
							if (m_haplo_data->m_allele_frequency[j][k] > 0) {
								new_as = new AlleleSequence;
								new_as->concatenate(*as, m_haplo_data->allele_symbol(j, k));
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
							new_as->concatenate(*as, genotype(1)[j]);
						}
						else {
							new_as->concatenate(*as, genotype(0)[j]);
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
						new_as->concatenate(*as, genotype(0)[j]);
						as_list.push_front(new_as);
						delete *i_as;
						i_as = as_list.erase(i_as);
					}
				}

			}
			i_as = as_list.begin();
			while (i_as != as_list.end()) {
				hp = m_pattern_tree->findLongestMatchPattern(m_head_len, *i_as);
				if (hp != NULL && hp->m_id >= (*head)->m_id) addHaploPair(m_last_list, new (HaploPair_pool.malloc()) HaploPair(*head, hp, m_target_pattern));
				delete *i_as;
				i_as = as_list.erase(i_as);
			}
		}
	}
}

void HaploBuilder::extendAll(Allele a1, Allele a2)
{
	vector<HaploPair*>::iterator i_hp;
	for (i_hp = m_last_list.begin(); i_hp != m_last_list.end(); i_hp++) {
		extend(*i_hp, a1, a2);
		if (a1 != a2) extend(*i_hp, a2, a1);
	}
}

void HaploBuilder::extend(HaploPair *hp, Allele a1, Allele a2)
{
	HaploPair *new_hp;
	double likelihood, total_likelihood;
	int b1, b2, id[2], hp_index;
	vector<pair<int, int> >::iterator i_p;
	b1 = m_haplo_data->getAlleleIndex(hp->end(), a1);
	b2 = m_haplo_data->getAlleleIndex(hp->end(), a2);
	if (hp->extendable(b1, b2, m_target_pattern)) {
		hp->extend_trial(b1, b2, id, likelihood, total_likelihood);
		hp_index = -1;
		for (i_p = m_best_pair[id[0]].begin(); i_p != m_best_pair[id[0]].end(); i_p++) {
			if (i_p->first == id[1]) {
				hp_index = i_p->second;
				break;
			}
		}
		if (hp_index < 0) {
			new_hp = new (HaploPair_pool.malloc()) HaploPair(hp, hp->successor_a(b1), hp->successor_b(b2));
			addHaploPair(m_new_list, new_hp);
		}
		else if (likelihood > m_new_list[hp_index]->m_likelihood) {
			new_hp = new (HaploPair_pool.malloc()) HaploPair(hp, hp->successor_a(b1), hp->successor_b(b2));
			new_hp->m_total_likelihood += m_new_list[hp_index]->m_total_likelihood;
			del_HaploPair(m_new_list[hp_index]);
			m_new_list[hp_index] = new_hp;
		}
		else {
			m_new_list[hp_index]->m_total_likelihood += total_likelihood;
		}
	}
}

void HaploBuilder::addHaploPair(vector<HaploPair*> &hp_list, HaploPair *hp)
{
	hp_list.push_back(hp);
	m_best_pair[hp->m_id[0]].push_back(make_pair(hp->m_id[1], hp_list.size()-1));
}

void HaploBuilder::findHaploPatternByFreq(double min_freq, int min_len, int max_len)
{
	int i;
	int *min_len_vector, *max_len_vector;
	vector<HaploPattern*> candidates;
	list<HaploPattern*>::iterator i_hp;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector = new int [m_genotype_len];
	max_len_vector = new int [m_genotype_len];
	for (i=0; i<m_genotype_len; i++) {
		min_len_vector[i] = min_len;
		max_len_vector[i] = max_len;
	}
	DeleteAllObjects(m_haplo_pattern);
	m_haplo_pattern.clear();
	generateHaploPatternCandidate(candidates);
	searchHaploPattern(candidates, min_freq, -1, min_len_vector, max_len_vector);
	DeleteAllObjects(candidates);
	candidates.clear();
	Logger::verbose("Found haplotype patterns: %d     ", m_haplo_pattern.size());
	delete[] min_len_vector;
	delete[] max_len_vector;
}

void HaploBuilder::findHaploPatternByNum(int max_num, int min_len, int max_len)
{
	int i, size, last_size;
	int *min_len_vector, *max_len_vector, *list_ind;
	double min_freq = 10, *list_freq;
	vector<HaploPattern*> candidates;
	list<HaploPattern*>::iterator i_hp;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector = new int [m_genotype_len];
	max_len_vector = new int [m_genotype_len];
	for (i=0; i<m_genotype_len; i++) {
		min_len_vector[i] = min_len;
		max_len_vector[i] = max_len;
	}
	DeleteAllObjects(m_haplo_pattern);
	m_haplo_pattern.clear();
	generateHaploPatternCandidate(candidates);
	do {
		last_size = m_haplo_pattern.size();
		min_freq *= 0.9;
		searchHaploPattern(candidates, min_freq, max_num, min_len_vector, max_len_vector);
	} while ((int)m_haplo_pattern.size() < max_num && min_freq > 1e-38);
	if ((int)m_haplo_pattern.size() > max_num) {
		size = m_haplo_pattern.size() - last_size;
		list_ind = new int [size];
		list_freq = new double [size];
		i_hp = m_haplo_pattern.begin();
		advance(i_hp, last_size);
		i = 0;
		while (i_hp != m_haplo_pattern.end()) {
			list_ind[i] = last_size+i;
			list_freq[i] = (*i_hp)->m_frequency;
			i_hp++;
			i++;
		}
		Utils::quick_sort_max(list_ind, list_freq, size);
		min_freq = list_freq[max_num-last_size];
		i_hp = m_haplo_pattern.begin();
		advance(i_hp, last_size);
		while (i_hp != m_haplo_pattern.end()) {
			if ((*i_hp)->m_frequency <= min_freq && (*i_hp)->m_length > 1) {
				delete *i_hp;
				i_hp = m_haplo_pattern.erase(i_hp);
			}
			else {
				i_hp++;
			}
		}
		delete[] list_ind;
		delete[] list_freq;
	}
	DeleteAllObjects(candidates);
	candidates.clear();
	Logger::verbose("Found haplotype patterns: %d     ", m_haplo_pattern.size());
	delete[] min_len_vector;
	delete[] max_len_vector;
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
	DeleteAllObjects(candidates);
	candidates.clear();
	for (i=0; i<m_genotype_len; i++) {
		hp = new HaploPattern(m_haplo_data, i);
		candidates.push_back(hp);
	}
}

void HaploBuilder::searchHaploPattern(vector<HaploPattern*> &candidates, double min_freq, int min_num, int *min_len, int *max_len)
{
	int i;
	vector<HaploPattern*> new_candidates;
	HaploPattern *hp, *hp_candidate;
	while (!candidates.empty()) {
		hp = candidates.back();
		candidates.pop_back();
		if (hp->m_frequency >= min_freq || hp->m_length <= 1) {
			if (hp->m_end < m_genotype_len && hp->m_length < max_len[hp->m_end]) {
				for (i=0; i<m_haplo_data->m_allele_num[hp->m_end]; i++) {
					if (m_haplo_data->m_allele_frequency[hp->m_end][i] > 0) {
						hp_candidate = new HaploPattern;
						hp_candidate->concatenate(*hp, m_haplo_data->m_allele_symbol[hp->m_end][i]);
						if (hp_candidate->m_frequency >= min_freq || hp->m_length <= 1)
						{
							candidates.push_back(hp_candidate);
						}
						else if (hp_candidate->m_frequency > 0)
						{
							if ((int)candidates.size() > Constant::ram_limit()) hp_candidate->releaseMatchGenotype();
							new_candidates.push_back(hp_candidate);
						}
						else {
							delete hp_candidate;
						}
					}
				}
			}
			if (hp->m_end > 0 && (hp->m_length >= min_len[hp->m_end-1] || hp->m_length == 1)) {
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
	int *min_len_vector, *max_len_vector;
	HaploPattern *hp;
	vector<HaploPattern*> candidates;
	list<HaploPattern*>::iterator i_hp;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector = new int [m_genotype_len];
	max_len_vector = new int [m_genotype_len];
	for (i=0; i<m_genotype_len; i++) {
		max_len_vector[i] = min_len+1;
	}
	for (i_hp = m_haplo_pattern.begin(); i_hp != m_haplo_pattern.end(); i_hp++) {
		hp = *i_hp;
		if (hp->m_length >= min_len && hp->m_length <= max_len && hp->m_length > max_len_vector[hp->m_end-1]) max_len_vector[hp->m_end-1] = hp->m_length;
		delete hp;
	}
	m_haplo_pattern.clear();
	for (i=0; i<m_genotype_len; i++) {
		min_len_vector[i] = max_len_vector[i]-1;
		Logger::debug("#%05d %5d %5d", i, min_len_vector[i], max_len_vector[i]);
	}
	generateHaploPatternCandidate(candidates);
	searchHaploPattern(candidates, 0, -1, min_len_vector, max_len_vector);
	DeleteAllObjects(candidates);
	candidates.clear();
	delete[] min_len_vector;
	delete[] max_len_vector;
}
