
#ifndef __PATTERNTREE_H
#define __PATTERNTREE_H


#include <vector>

#include "Tree.h"
#include "HaploPattern.h"


class AlleleSequence;
class GenoData;


typedef TreeNode<HaploPattern*> PatternNode;


class BackwardPatternTree {
	const GenoData &m_genos;
	vector<PatternNode> m_trees;

public:
	explicit BackwardPatternTree(const GenoData &genos);

	void addPattern(HaploPattern *hp) {	addPattern(&m_trees[hp->end()], hp, hp->length()); }

	HaploPattern *findLongestMatchPattern(int end, const HaploPattern *hp, int len = 0) const;
	HaploPattern *findLongestMatchPattern(int end, const AlleleSequence *as, int as_start = 0, int len = 0) const;
	HaploPattern *findLikelyMatchPattern(int end, const HaploPattern *hp, int len = 0) const;
	HaploPattern *findLikelyMatchPattern(int end, const AlleleSequence *as, int as_start = 0, int len = 0) const;

	HaploPattern *getSingleAllelePattern(int end, int index) const;

protected:
	void addPattern(PatternNode *node, HaploPattern *hp, int len);
	HaploPattern *findLongestMatchPattern(const PatternNode *node, const AlleleSequence *as, int lg, int ll, int len) const;
	HaploPattern *findLikelyMatchPattern(const PatternNode *node, const AlleleSequence *as, int lg, int ll, int len) const;
};

inline HaploPattern *BackwardPatternTree::getSingleAllelePattern(int end, int index) const
{
	return m_trees[end].getChild(index)->data();
}


class ForwardPatternTree {
	const GenoData &m_genos;
	vector<PatternNode> m_trees;

public:
	explicit ForwardPatternTree(const GenoData &genos);

	PatternNode *root(int i) { return &m_trees[i]; }

	void addPattern(HaploPattern *hp) {	addPattern(&m_trees[hp->start()], hp, hp->length()); }

protected:
	void addPattern(PatternNode *node, HaploPattern *hp, int len);
};


#endif // __PATTERNTREE_H
