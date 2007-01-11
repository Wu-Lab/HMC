
#ifndef __HAPLOBUILDER_H
#define __HAPLOBUILDER_H


#include <list>
#include <vector>

#include "Utils.h"
#include "Constant.h"
#include "Allele.h"


class AlleleSequence;
class Haplotype;
class Genotype;
class HaploPattern;
class HaploPair;
class HaploData;
typedef TreeNode<HaploPattern, double> PatternNode;


class PatternTree {
	const HaploData *m_haplo_data;
	int m_genotype_len;
	Tree<HaploPattern, double> *m_trees;

public:
	PatternTree(const HaploData *haplo);
	~PatternTree();

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
	return m_trees[end].root().getChild(index);
}


class HaploBuilder {
protected:
	HaploData *m_haplo_data;
	list<HaploPattern*> m_haplo_pattern;
	int m_genotype_len;
	int m_pattern_num;

	PatternTree *m_pattern_tree;
	list<HaploPattern*> m_head_list;
	int m_head_len;

	Haplotype *m_genotype;
	vector<pair<int, int> > *m_best_pair;
	vector<HaploPair*> m_last_list;
	vector<HaploPair*> m_new_list;
	HaploPattern *m_target_pattern;

public:
	HaploBuilder();
	~HaploBuilder();

	HaploData *haplo_data() { return m_haplo_data; }
	list<HaploPattern*> &haplo_pattern() { return m_haplo_pattern; }

	void setHaploData(HaploData &hd);

	void initialize();
	void adjust(double min_freq = 0);

	void resolve(const Genotype &genotype, Genotype &resolution, vector<HaploPair*> &res_list, HaploPattern *target_pattern = NULL);

	double getLikelihood(const Haplotype &haplotype);
	double getLikelihood(const Genotype &genotype);

protected:
	void clear();
	void initHeadList(const Genotype &genotype);

	void extendAll(Allele a1, Allele a2);
	void extend(HaploPair *hp, Allele a1, Allele a2);
	void addHaploPair(vector<HaploPair*> &hp_list, HaploPair *hp);

	void findHaploPatternByFreq(double min_freq, int min_len = 1, int max_len = -1);
	void findHaploPatternByNum(int max_num, int min_len = 1, int max_len = -1);
	void findHaploPatternBlockByFreq(double min_freq, int min_len = 1, int max_len = -1);
	void findHaploPatternBlockByNum(int max_num, int min_len = 1, int max_len = -1);

	void generateHaploPatternCandidate(vector<HaploPattern*> &candidates);
	void searchHaploPattern(vector<HaploPattern*> &candidates, double min_freq, int min_num, int *min_len, int *max_len);
	void searchHaploPatternBlock(int min_len, int max_len);

	friend class HaploData;
};


#endif // __HAPLOBUILDER_H
