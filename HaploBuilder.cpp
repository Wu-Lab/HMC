
#include "HaploBuilder.h"


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
	m_heads = NULL;
	m_best_pair = NULL;
	m_target_pattern = NULL;
}

HaploBuilder::~HaploBuilder()
{
	delete m_pattern_tree;
	delete[] m_heads;
	delete[] m_best_pair;
}

void HaploBuilder::clear()
{
	delete m_pattern_tree;
	delete[] m_heads;
	delete[] m_best_pair;
	m_pattern_tree = NULL;
	m_heads = NULL;
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
	clear();
	m_pattern_num = m_haplo_pattern.size();
	m_pattern_tree = new PatternTree(m_haplo_data);
	List<HaploPattern> head_list;
	m_head_len = m_genotype_len;
	hp = m_haplo_pattern.getFirst();
	i = 0;
	while (hp != NULL) {
		hp->m_id = i++;
		m_pattern_tree->addPattern(hp->m_end, hp);
		if (hp->m_start == 0) {
			head_list.addLast(hp);
			if (hp->m_length < m_head_len) m_head_len = hp->m_length;
		}
		hp = m_haplo_pattern.getNext();
	}
	hp = head_list.getFirst();
	while (hp != NULL) {
		if (hp->m_length > m_head_len) {
			head_list.releaseCurrent();
		}
		else {
			head_list.moveNext();
		}
		hp = head_list.getCurrent();
	}
	m_head_num = head_list.toArray(m_heads);
	head_list.releaseAll();
	hp = m_haplo_pattern.getFirst();
	while (hp != NULL) {
		hp->m_prefix = m_pattern_tree->findLongestMatchPattern(hp->m_end-1, hp);
		if (hp->m_prefix != NULL) {
			hp->m_transition_prob = hp->m_frequency / hp->m_prefix->m_frequency;
		}
		else {
			hp->m_transition_prob = hp->m_frequency / m_haplo_data->m_genotype_num;
		}
		hp->m_weight = 0;
		hp = m_haplo_pattern.getNext();
	}
	hp = m_haplo_pattern.getFirst();
	while (hp != NULL) {
		if (hp->m_end < m_genotype_len) {
			temp.concatenate(*hp, -1);			// append empty allele
			for (j=0; j<m_haplo_data->allele_num(hp->m_end); j++) {
				temp.allele(temp.length()-1) = m_haplo_data->allele_symbol(hp->m_end, j);
				hp->m_successors[j] = m_pattern_tree->findLongestMatchPattern(hp->m_end+1, &temp);
			}
		}
		hp = m_haplo_pattern.getNext();
	}
	m_best_pair = new List<ListItem<HaploPair>, int> [m_pattern_num];
}

void HaploBuilder::adjust()
{
	int i;
	double *new_frequency;
	HaploPattern *hp;
	Genotype res;
	List<HaploPair, double> res_list;
	new_frequency = new double [m_pattern_num];
	hp = m_haplo_pattern.getFirst();
	while (hp != NULL) {
		new_frequency[hp->id()] = 0;
		for (i=0; i<m_haplo_data->m_genotype_num; i++) {
			if (hp->isMatch(m_haplo_data->genotype(i))) {
				resolve(m_haplo_data->genotype(i), res, res_list, hp);
				new_frequency[hp->id()] += res.likelihood() / m_haplo_data->genotype(i).likelihood();
			}
		}
		Logger::status("Adjust haplotype patterns: %d     ", hp->id());
		hp = m_haplo_pattern.getNext();
	}
	hp = m_haplo_pattern.getFirst();
	while (hp != NULL) {
		hp->m_frequency = new_frequency[hp->id()];
		hp = m_haplo_pattern.getNext();
	}
	hp = m_haplo_pattern.getFirst();
	while (hp != NULL) {
		if (hp->m_prefix != NULL) {
			if (hp->m_frequency > hp->m_prefix->m_frequency) {
				Logger::error("Mistake in adjusting the pattern frequencies!");
				exit(1);
			}
			if (hp->m_prefix->m_frequency > 1e-6) {
				hp->m_transition_prob = hp->m_frequency / hp->m_prefix->m_frequency;
			}
			else {
				hp->m_transition_prob = 0;
			}
		}
		else {
			hp->m_transition_prob = hp->m_frequency / m_haplo_data->m_genotype_num;
		}
		hp = m_haplo_pattern.getNext();
	}
	delete[] new_frequency;
}

void HaploBuilder::resolve(const Genotype &genotype, Genotype &resolution, List<HaploPair, double> &res_list, HaploPattern *target_pattern)
{
	int i, j, k;
	Allele a, b;
	double weight, total_likelihood;
	HaploPair *hp, *best;
	m_target_pattern = target_pattern;
	m_last_list = new List<HaploPair>;
	m_new_list = new List<HaploPair>;
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
		delete m_last_list;
		m_last_list = m_new_list;
		m_new_list = new List<HaploPair>;
		if (m_last_list->size() <= 0) {
			break;
		}
	}
	for (i=0; i<m_pattern_num; i++) {
		m_best_pair[i].releaseAll();
	}
	if (m_last_list->size() > 0) {
		weight = m_haplo_data->m_genotype_num * m_haplo_data->m_genotype_num;
		total_likelihood = 0;
		res_list.removeAll();
		hp = m_last_list->getFirst();
		while (hp != NULL) {
			hp->m_likelihood /= weight;
			hp->m_total_likelihood /= weight;
			total_likelihood += hp->m_total_likelihood;
			res_list.addDescent(hp, hp->m_likelihood);
			hp = m_last_list->getNext();
		}
		m_last_list->releaseAll();
		best = res_list.getFirst();
		resolution = best->getGenotype();
		resolution.setLikelihood(total_likelihood);
		resolution.setWeight(best->likelihood() / total_likelihood);
	}
	else {
		res_list.removeAll();
		resolution = genotype;
		resolution.setLikelihood(0);
		resolution.setWeight(0);
	}
	delete m_last_list;
	delete m_new_list;
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
	int i, j, k;
	AlleleSequence *as, *new_as;
	HaploPattern *hp;
	List<AlleleSequence> as_list;
	for (i=0; i<m_head_num; i++) {
		if (m_heads[i]->isMatch(genotype)) {
			as_list.addLast(new AlleleSequence);
			for (j=0; j<m_head_len; j++) {
				as = as_list.getFirst();
				if (genotype.isMissing(j) || (genotype.hasMissing(j) && genotype.hasAllele(j, (*m_heads[i])[j]))) {
					while (as != NULL) {
						for (k=0; k<m_haplo_data->allele_num(j); k++) {
							if (m_haplo_data->m_allele_frequency[j][k] > 0) {
								new_as = new AlleleSequence;
								new_as->concatenate(*as, m_haplo_data->allele_symbol(j, k));
								as_list.addFirst(new_as);
							}
						}
						as_list.removeCurrent();
						as = as_list.getCurrent();
					}
				}
				else if (genotype.isHeterozygous(j)) {
					while (as != NULL) {
						new_as = new AlleleSequence;
						if ((*m_heads[i])[j] == genotype(0)[j]) {
							new_as->concatenate(*as, genotype(1)[j]);
						}
						else {
							new_as->concatenate(*as, genotype(0)[j]);
						}
						as_list.addFirst(new_as);
						as_list.removeCurrent();
						as = as_list.getCurrent();
					}
				}
				else {
					while (as != NULL) {
						new_as = new AlleleSequence;
						new_as->concatenate(*as, genotype(0)[j]);
						as_list.addFirst(new_as);
						as_list.removeCurrent();
						as = as_list.getCurrent();
					}
				}

			}
			as = as_list.getFirst();
			while (as != NULL) {
				hp = m_pattern_tree->findLongestMatchPattern(m_head_len, as);
				if (hp != NULL && hp->m_id >= m_heads[i]->m_id) addHaploPair(m_last_list, new HaploPair(m_heads[i], hp, m_target_pattern));
				as = as_list.getNext();
			}
			as_list.removeAll();
		}
	}
}

void HaploBuilder::extendAll(Allele a1, Allele a2)
{
	HaploPair *hp;
	hp = m_last_list->getFirst();
	while (hp != NULL) {
		extend(hp, a1, a2);
		if (a1 != a2) extend(hp, a2, a1);
		hp = m_last_list->getNext();
	}
}

void HaploBuilder::extend(HaploPair *hp, Allele a1, Allele a2)
{
	HaploPair *new_hp;
	int b1, b2;
	b1 = m_haplo_data->getAlleleIndex(hp->end(), a1);
	b2 = m_haplo_data->getAlleleIndex(hp->end(), a2);
	if (hp->extendable(b1, b2, m_target_pattern)) {
		new_hp = new HaploPair(*hp);
		new_hp->extend(b1, b2);
		addHaploPair(m_new_list, new_hp);
	}
}

void HaploBuilder::addHaploPair(List<HaploPair> *hp_list, HaploPair *hp)
{
	int id[2];
	ListItem<HaploPair> *hp_item;
	hp->getID(id);
	hp_item = m_best_pair[id[0]].get(id[1]);
	if (hp_item == NULL) {
		hp_list->addLast(hp);
		m_best_pair[id[0]].addLast(hp_list->last(), id[1]);
	}
	else if (hp->m_likelihood > hp_item->object()->m_likelihood) {
		hp->m_total_likelihood += hp_item->object()->m_total_likelihood;
		hp_item->removeObject();
		hp_item->setObject(hp);
	}
	else {
		hp_item->object()->m_total_likelihood += hp->m_total_likelihood;
		delete hp;
	}
}

void HaploBuilder::findHaploPatternByFreq(double min_freq, int min_len, int max_len)
{
	int i;
	int *min_len_vector, *max_len_vector;
	List<HaploPattern> candidates;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector = new int [m_genotype_len];
	max_len_vector = new int [m_genotype_len];
	for (i=0; i<m_genotype_len; i++) {
		min_len_vector[i] = min_len;
		max_len_vector[i] = max_len;
	}
	m_haplo_pattern.removeAll();
	generateHaploPatternCandidate(candidates);
	searchHaploPattern(candidates, min_freq, -1, min_len_vector, max_len_vector);
	Logger::verbose("Found haplotype patterns: %d     ", m_haplo_pattern.size());
	delete[] min_len_vector;
	delete[] max_len_vector;
}

void HaploBuilder::findHaploPatternByNum(int max_num, int min_len, int max_len)
{
	int i, size, top_size;
	int *min_len_vector, *max_len_vector, *list_ind;
	double min_freq = 10, *list_freq;
	HaploPattern *hp;
	List<HaploPattern> candidates;
	ListItem<HaploPattern> *top;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector = new int [m_genotype_len];
	max_len_vector = new int [m_genotype_len];
	for (i=0; i<m_genotype_len; i++) {
		min_len_vector[i] = min_len;
		max_len_vector[i] = max_len;
	}
	m_haplo_pattern.removeAll();
	generateHaploPatternCandidate(candidates);
	do {
		top = m_haplo_pattern.last();
		top_size = m_haplo_pattern.size();
		min_freq *= 0.9;
		searchHaploPattern(candidates, min_freq, max_num, min_len_vector, max_len_vector);
	} while (m_haplo_pattern.size() < max_num && min_freq > 1e-38);
	if (m_haplo_pattern.size() > max_num) {
		size = m_haplo_pattern.size() - top_size;
		list_ind = new int [size];
		list_freq = new double [size];
		m_haplo_pattern.getItem(top);
		hp = m_haplo_pattern.getNext();
		i = 0;
		while (hp != NULL) {
			list_ind[i] = i;
			list_freq[i] = hp->m_frequency;
			hp = m_haplo_pattern.getNext();
			i++;
		}
		Utils::quick_sort_max(list_ind, list_freq, size);
		min_freq = list_freq[max_num-top_size];
		m_haplo_pattern.getItem(top);
		hp = m_haplo_pattern.getNext();
		while (hp != NULL) {
			if (hp->m_frequency <= min_freq && hp->m_length > 1) {
				m_haplo_pattern.removeCurrent();
			}
			else {
				m_haplo_pattern.moveNext();
			}
			hp = m_haplo_pattern.getCurrent();
		}
		delete[] list_ind;
		delete[] list_freq;
	}
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

void HaploBuilder::generateHaploPatternCandidate(List<HaploPattern> &candidates)
{
	int i;
	HaploPattern *hp;
	candidates.removeAll();
	for (i=0; i<m_genotype_len; i++) {
		hp = new HaploPattern(m_haplo_data, i);
		candidates.addLast(hp);
	}
}

void HaploBuilder::searchHaploPattern(List<HaploPattern> &candidates, double min_freq, int min_num, int *min_len, int *max_len)
{
	int i;
	HaploPattern *hp, *hp_candidate;
	hp = candidates.getFirst();
	while (hp != NULL) {
		if (hp->m_frequency >= min_freq || hp->m_length <= 1) {
			if (hp->m_end < m_genotype_len && hp->m_length < max_len[hp->m_end]) {
				for (i=0; i<m_haplo_data->m_allele_num[hp->m_end]; i++) {
					if (m_haplo_data->m_allele_frequency[hp->m_end][i] > 0) {
						hp_candidate = new HaploPattern;
						hp_candidate->concatenate(*hp, m_haplo_data->m_allele_symbol[hp->m_end][i]);
						if (hp_candidate->m_frequency >= min_freq || hp->m_length <= 1)
						{
							candidates.addNext(hp_candidate);
						}
						else if (hp_candidate->m_frequency > 0)
						{
							if (candidates.size() > Constant::ram_limit()) hp_candidate->releaseMatchGenotype();
							candidates.addFirst(hp_candidate);
						}
						else {
							delete hp_candidate;
						}
					}
				}
			}
			if (hp->m_end > 0 && hp->m_length >= min_len[hp->m_end-1]) {
				hp->releaseMatchGenotype();
				candidates.releaseCurrent();
				m_haplo_pattern.addLast(hp);
				Logger::status("Found haplotype patterns: %d     ", m_haplo_pattern.size());
			}
			else {
				candidates.removeCurrent();
			}
		}
		else if (min_num < 0) {
			candidates.removeCurrent();
		}
		else {
			candidates.moveNext();
		}
		hp = candidates.getCurrent();
	}
}

void HaploBuilder::searchHaploPatternBlock(int min_len, int max_len)
{
	int i;
	int *min_len_vector, *max_len_vector;
	HaploPattern *hp;
	List<HaploPattern> candidates;
	if (min_len < 0) min_len = 1;
	if (max_len < 0) max_len = m_genotype_len;
	min_len_vector = new int [m_genotype_len];
	max_len_vector = new int [m_genotype_len];
	for (i=0; i<m_genotype_len; i++) {
		max_len_vector[i] = min_len+1;
	}
	hp = m_haplo_pattern.getFirst();
	while (hp != NULL) {
		if (hp->m_length >= min_len && hp->m_length <= max_len && hp->m_length > max_len_vector[hp->m_end-1]) max_len_vector[hp->m_end-1] = hp->m_length;
		hp = m_haplo_pattern.getNext();
	}
	for (i=0; i<m_genotype_len; i++) {
		min_len_vector[i] = max_len_vector[i]-1;
		Logger::debug("#%05d %5d %5d", i, min_len_vector[i], max_len_vector[i]);
	}
	m_haplo_pattern.removeAll();
	generateHaploPatternCandidate(candidates);
	searchHaploPattern(candidates, 0, -1, min_len_vector, max_len_vector);
	delete[] min_len_vector;
	delete[] max_len_vector;
}
