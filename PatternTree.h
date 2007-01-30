
#ifndef __PATTERNTREE_H
#define __PATTERNTREE_H


#include <vector>

#include "Tree.h"


class AlleleSequence;
class HaploPattern;
class HaploData;


typedef TreeNode<HaploPattern*> PatternNode;


class PatternTree {
	const HaploData *m_haplo_data;
	int m_genotype_len;
	vector<PatternNode> m_trees;

public:
	explicit PatternTree(const HaploData *haplo);

	void addPattern(int end, HaploPattern *hp);
	HaploPattern *findLongestMatchPattern(int end, const HaploPattern *hp, int len = 0);
	HaploPattern *findLongestMatchPattern(int end, const AlleleSequence *as, int len = 0, int offset = 0);
	HaploPattern *findLikelyMatchPattern(int end, const HaploPattern *hp, int len = 0);
	HaploPattern *findLikelyMatchPattern(int end, const AlleleSequence *as, int len = 0, int offset = 0);

	HaploPattern *getSingleAllelePattern(int end, int index);

protected:
	void addPattern(PatternNode *node, HaploPattern *hp, int len);
	HaploPattern *findLongestMatchPattern(PatternNode *node, const AlleleSequence *as, int end, int len, int offset = 0);
	HaploPattern *findLikelyMatchPattern(PatternNode *node, const AlleleSequence *as, int end, int len, int offset = 0);

	friend class HaploBuilder;
};

inline HaploPattern *PatternTree::getSingleAllelePattern(int end, int index)
{
	return m_trees[end].getChild(index)->data();
}


#endif // __PATTERNTREE_H
