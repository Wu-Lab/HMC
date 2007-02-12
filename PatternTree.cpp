
#include "PatternTree.h"
#include "HaploPattern.h"
#include "HaploData.h"

#include "MemLeak.h"


boost::pool<> PatternNode::m_pool(sizeof(PatternNode));

////////////////////////////////
//
// class BackwardPatternTree

BackwardPatternTree::BackwardPatternTree(const HaploData &hd)
: m_haplodata(hd),
  m_trees(hd.genotype_len()+1)
{
	for (int i=0; i<=hd.genotype_len(); ++i) {
		m_trees[i].resize(hd.max_allele_num());
	}
}

void BackwardPatternTree::addPattern(PatternNode *node, HaploPattern *hp, int len)
{
	int i = hp->getAlleleIndex(len-1);
	if (i < 0) {										// allele is missing
		int n = m_haplodata.allele_num(hp->getGlobalLocus(len-1));
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

HaploPattern *BackwardPatternTree::findLongestMatchPattern(int end, const HaploPattern *hp, int len) const
{
	if (end <= hp->start() || end > hp->end()) {
		return 0;
	}
	else {
		int max_len = end - hp->start();
		if (len <= 0 || len > max_len) len = max_len;
		return findLongestMatchPattern(&m_trees[end], hp, end-1, max_len-1, len);
	}
}

HaploPattern *BackwardPatternTree::findLongestMatchPattern(int end, const AlleleSequence *as, int as_start, int len) const
{
	if (end <= as_start || end > (as_start+as->length())) {
		return 0;
	}
	else {
		int max_len = end - as_start;
		if (len <= 0 || len > max_len) len = max_len;
		return findLongestMatchPattern(&m_trees[end], as, end-1, max_len-1, len);
	}
}

HaploPattern *BackwardPatternTree::findLikelyMatchPattern(int end, const HaploPattern *hp, int len) const
{
	if (end <= hp->start() || end > hp->end()) {
		return 0;
	}
	else {
		int max_len = end - hp->start();
		if (len <= 0 || len > max_len) len = max_len;
		return findLikelyMatchPattern(&m_trees[end], hp, end-1, max_len-1, len);
	}
}

HaploPattern *BackwardPatternTree::findLikelyMatchPattern(int end, const AlleleSequence *as, int as_start, int len) const
{
	if (end <= as_start || end > (as_start+as->length())) {
		return 0;
	}
	else {
		int max_len = end - as_start;
		if (len <= 0 || len > max_len) len = max_len;
		return findLikelyMatchPattern(&m_trees[end], as, end-1, max_len-1, len);
	}
}

HaploPattern *BackwardPatternTree::findLongestMatchPattern(const PatternNode *node, const AlleleSequence *as, int lg, int ll, int len) const
{
	int i, n;
	HaploPattern *result, *temp;
	result = node->data();
	if ((*as)[ll].isMissing()) {						// allele is missing
		n = m_haplodata.allele_num(lg);
		for (i=0; i<n; i++) {
			if (node->getChild(i)) {					// previous locus is matching
				if (len > 1) {
					temp = findLongestMatchPattern(node->getChild(i), as, lg-1, ll-1, len-1);
				}
				else {
					temp = node->getChild(i)->data();
				}
				if (result == 0 || (temp && temp->length() > result->length())) {
					result = temp;
				}
			}
		}
	}
	else {
		i = m_haplodata.getAlleleIndex(lg, (*as)[ll]);
		if (node->getChild(i)) {								// previous locus is matching
			if (len > 1) {
				temp = findLongestMatchPattern(node->getChild(i), as, lg-1, ll-1, len-1);
			}
			else {
				temp = node->getChild(i)->data();
			}
			if (result == 0 || (temp && temp->length() > result->length())) {
				result = temp;
			}
		}
	}
	return result;
}

HaploPattern *BackwardPatternTree::findLikelyMatchPattern(const PatternNode *node, const AlleleSequence *as, int lg, int ll, int len) const
{
	int i, n;
	HaploPattern *result, *temp;
	result = node->data();
	if ((*as)[ll].isMissing()) {						// allele is missing
		n = m_haplodata.allele_num(lg);
		for (i=0; i<n; i++) {
			if (node->getChild(i)) {					// previous locus is matching
				if (len > 1) {
					temp = findLikelyMatchPattern(node->getChild(i), as, lg-1, ll-1, len-1);
				}
				else {
					temp = node->getChild(i)->data();
				}
				if (result == 0 || (temp && temp->transition_prob() >= result->transition_prob())) {
					result = temp;
				}
			}
		}
	}
	else {
		i = m_haplodata.getAlleleIndex(lg, (*as)[ll]);
		if (node->getChild(i)) {						// previous locus is matching
			if (len > 1) {
				temp = findLikelyMatchPattern(node->getChild(i), as, lg-1, ll-1, len-1);
			}
			else {
				temp = node->getChild(i)->data();
			}
			if (result == 0 || (temp && temp->transition_prob() >= result->transition_prob())) {
				result = temp;
			}
		}
	}
	return result;
}


////////////////////////////////
//
// class ForwardPatternTree

ForwardPatternTree::ForwardPatternTree(const HaploData &hd)
: m_haplodata(hd),
  m_trees(hd.genotype_len()+1)
{
	for (int i=0; i<=hd.genotype_len(); ++i) {
		m_trees[i].resize(hd.max_allele_num());
	}
}

void ForwardPatternTree::addPattern(PatternNode *node, HaploPattern *hp, int len)
{
	int i = hp->getAlleleIndex(hp->length()-len);
	if (i < 0) {										// allele is missing
		int n = m_haplodata.allele_num(hp->getGlobalLocus(hp->length()-len));
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
