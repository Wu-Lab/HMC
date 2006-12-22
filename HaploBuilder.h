
#ifndef __HAPLOBUILDER_H
#define __HAPLOBUILDER_H


#include "Utils.h"
#include "Constant.h"
#include "Allele.h"
#include "Haplotype.h"
#include "Genotype.h"
#include "HaploData.h"
#include "HaploPattern.h"
#include "HaploPair.h"


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
	List<HaploPattern> m_haplo_pattern;
	int m_genotype_len;
	int m_pattern_num;

	PatternTree *m_pattern_tree;
	HaploPattern **m_heads;
	int m_head_num;
	int m_head_len;

	Haplotype *m_genotype;
	List<ListItem<HaploPair>, int> *m_best_pair;
	List<HaploPair> *m_last_list;
	List<HaploPair> *m_new_list;

public:
	HaploBuilder();
	~HaploBuilder();

	HaploData *haplo_data() { return m_haplo_data; }
	List<HaploPattern> &haplo_pattern() { return m_haplo_pattern; }

	void setHaploData(HaploData &hd);

	void initialize();
	void adjust(const Haplotype &h_old, const Haplotype &h_new);
	void adjust(const Genotype &g_old, const Genotype &g_new);

	void resolve(const Genotype &genotype, Genotype &resolution, List<HaploPair, double> &res_list);

	double getLikelihood(const Haplotype &haplotype);
	double getLikelihood(const Genotype &genotype);

protected:
	void clear();
	void initHeadList(const Genotype &genotype);

	void extendAll(Allele a1, Allele a2);
	void extend(HaploPair *hp, Allele a1, Allele a2);
	void addHaploPair(List<HaploPair> *hp_list, HaploPair *hp);

	void findHaploPatternByFreq(double min_freq, int min_len = 1, int max_len = -1);
	void findHaploPatternByNum(int max_num, int min_len = 1, int max_len = -1);
	void findHaploPatternBlockByFreq(double min_freq, int min_len = 1, int max_len = -1);
	void findHaploPatternBlockByNum(int max_num, int min_len = 1, int max_len = -1);

	void generateHaploPatternCandidate(List<HaploPattern> &candidates);
	void searchHaploPattern(List<HaploPattern> &candidates, double min_freq, int min_num, int *min_len, int *max_len);
	void searchHaploPatternBlock(int min_len, int max_len);

	friend class HaploData;
};


#endif // __HAPLOBUILDER_H
