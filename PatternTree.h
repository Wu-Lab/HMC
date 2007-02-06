
#ifndef __PATTERNTREE_H
#define __PATTERNTREE_H


#include <vector>

#include "Tree.h"
#include "HaploPattern.h"


class AlleleSequence;
class HaploData;


typedef TreeNode<HaploPattern*> PatternNode;


class BackwardPatternTree {
	const HaploData &m_haplodata;
	vector<PatternNode> m_trees;

public:
	explicit BackwardPatternTree(const HaploData &hd);

	void addPattern(HaploPattern *hp) {	addPattern(&m_trees[hp->end()], hp, hp->length()); }

	HaploPattern *findLongestMatchPattern(int end, const HaploPattern *hp, int len = 0);
	HaploPattern *findLongestMatchPattern(int end, const AlleleSequence *as, int as_start = 0, int len = 0);
	HaploPattern *findLikelyMatchPattern(int end, const HaploPattern *hp, int len = 0);
	HaploPattern *findLikelyMatchPattern(int end, const AlleleSequence *as, int as_start = 0, int len = 0);

	HaploPattern *getSingleAllelePattern(int end, int index);

protected:
	void addPattern(PatternNode *node, HaploPattern *hp, int len);
	HaploPattern *findLongestMatchPattern(PatternNode *node, const AlleleSequence *as, int lg, int ll, int len);
	HaploPattern *findLikelyMatchPattern(PatternNode *node, const AlleleSequence *as, int lg, int ll, int len);
};

inline HaploPattern *BackwardPatternTree::getSingleAllelePattern(int end, int index)
{
	return m_trees[end].getChild(index)->data();
}


class ForwardPatternTree {
	const HaploData &m_haplodata;
	vector<PatternNode> m_trees;

public:
	explicit ForwardPatternTree(const HaploData &hd);

	void addPattern(HaploPattern *hp) {	addPattern(&m_trees[hp->start()], hp, hp->length()); }

protected:
	void addPattern(PatternNode *node, HaploPattern *hp, int len);
};


#endif // __PATTERNTREE_H
