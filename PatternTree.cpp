
#include "PatternTree.h"
#include "HaploPattern.h"
#include "HaploData.h"

#include "MemLeak.h"


PatternTree::PatternTree(const HaploData *haplo)
: m_haplo_data(haplo),
  m_genotype_len(haplo->genotype_len()),
  m_trees(haplo->genotype_len()+1)
{
	for (int i=0; i<=m_genotype_len; ++i) {
		m_trees[i].resize(haplo->max_allele_num());
	}
}

void PatternTree::addPattern(int end, HaploPattern *hp)
{
	if (hp->end() != end) {
		Logger::error("Attempt to add HaploPattern to incorrect PatternTree!");
		exit(1);
	}
	else {
		addPattern(&m_trees[end], hp, hp->length());
	}
}

void PatternTree::addPattern(PatternNode *node, HaploPattern *hp, int len)
{
	int i = hp->getAlleleIndex(len-1);
	if (i < 0) {										// allele is missing
		int n = m_haplo_data->allele_num(hp->getGlobalLocus(len-1));
		if (len == 1) {
			for (int j=0; j<n; ++j) {
				node->setChild(j, hp);
			}
		}
		else {
			for (int j=0; j<n; ++j) {
				addPattern(node->addChild(j), hp, len-1);
			}
		}
	}
	else {
		if (len == 1) {
			node->setChild(i, hp);
		}
		else {
			addPattern(node->addChild(i), hp, len-1);
		}
	}
}

HaploPattern *PatternTree::findLongestMatchPattern(int end, const HaploPattern *hp, int len)
{
	if (hp->end() < end || hp->start() >= end) {
		return NULL;
	}
	else {
		if (len <= 0) len = end - hp->start();
		return findLongestMatchPattern(&m_trees[end], hp, end, len, end-hp->start()-len);
	}
}

HaploPattern *PatternTree::findLongestMatchPattern(int end, const AlleleSequence *as, int len, int offset)
{
	if (len <= 0) len = as->length();
	return findLongestMatchPattern(&m_trees[end], as, end, len, offset);
}

HaploPattern *PatternTree::findLikelyMatchPattern(int end, const HaploPattern *hp, int len)
{
	if (hp->end() < end || hp->start() >= end) {
		return NULL;
	}
	else {
		if (len <= 0) len = end - hp->start();
		return findLikelyMatchPattern(&m_trees[end], hp, end, len, end-hp->start()-len);
	}
}

HaploPattern *PatternTree::findLikelyMatchPattern(int end, const AlleleSequence *as, int len, int offset)
{
	if (len <= 0) len = as->length();
	return findLikelyMatchPattern(&m_trees[end], as, end, len, offset);
}

HaploPattern *PatternTree::findLongestMatchPattern(PatternNode *node, const AlleleSequence *as, int end, int len, int offset)
{
	int i, n;
	HaploPattern *result, *temp;
	result = node->data();
	if ((*as)[offset+len-1].isMissing()) {						// allele is missing
		n = m_haplo_data->allele_num(end-1);
		for (i=0; i<n; i++) {
			if (node->getChild(i)) {							// previous locus is matching
				if (len > 1) {
					temp = findLongestMatchPattern(node->getChild(i), as, end-1, len-1, offset);
				}
				else {
					temp = node->getChild(i)->data();
				}
				if (result == NULL || (temp != NULL && temp->length() > result->length())
					|| (temp != NULL && temp->length() >= result->length() 
					&& temp->transition_prob() > result->transition_prob())) {
					result = temp;
				}
			}
		}
	}
	else {
		i = m_haplo_data->getAlleleIndex(end-1, (*as)[offset+len-1]);
		if (node->getChild(i)) {								// previous locus is matching
			if (len > 1) {
				temp = findLongestMatchPattern(node->getChild(i), as, end-1, len-1, offset);
			}
			else {
				temp = node->getChild(i)->data();
			}
			if (result == NULL || (temp != NULL && temp->length() > result->length())
				|| (temp != NULL && temp->length() >= result->length() 
				&& temp->transition_prob() > result->transition_prob())) {
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
	result = node->data();
	if ((*as)[offset+len-1].isMissing()) {						// allele is missing
		n = m_haplo_data->allele_num(end-1);
		for (i=0; i<n; i++) {
			if (node->getChild(i)) {							// previous locus is matching
				if (len > 1) {
					temp = findLikelyMatchPattern(node->getChild(i), as, end-1, len-1, offset);
				}
				else {
					temp = node->getChild(i)->data();
				}
				if (result == NULL || (temp != NULL && temp->transition_prob() >= result->transition_prob())) {
					result = temp;
				}
			}
		}
	}
	else {
		i = m_haplo_data->getAlleleIndex(end-1, (*as)[offset+len-1]);
		if (node->getChild(i)) {								// previous locus is matching
			if (len > 1) {
				temp = findLikelyMatchPattern(node->getChild(i), as, end-1, len-1, offset);
			}
			else {
				temp = node->getChild(i)->data();
			}
			if (result == NULL || (temp != NULL && temp->transition_prob() >= result->transition_prob())) {
				result = temp;
			}
		}
	}
	return result;
}
