
#ifndef __HAPLOBUILDER_H
#define __HAPLOBUILDER_H


#include <list>
#include <vector>
#include <map>

#include "Utils.h"
#include "Constant.h"
#include "Allele.h"
#include "HaploPattern.h"
#include "HaploPair.h"
#include "PatternTree.h"


class HaploBuilder {
protected:
	HaploData *m_haplo_data;
	vector<HaploPattern*> m_haplo_pattern;
	int m_genotype_len;
	int m_pattern_num;

	tr1::shared_ptr<PatternTree> m_pattern_tree;
	list<HaploPattern*> m_head_list;
	int m_head_len;

	Haplotype *m_genotype;
	vector<vector<HaploPair*> > m_haplopairs;
	vector<map<int, int> > m_best_pair;

public:
	HaploBuilder();
	~HaploBuilder();

	HaploData *haplo_data() { return m_haplo_data; }
	const vector<HaploPattern*> &haplo_pattern() { return m_haplo_pattern; }

	void setHaploData(HaploData &hd);

	void initialize();
	void adjust(double min_freq = 0);

	void resolve(const Genotype &genotype, Genotype &resolution, vector<HaploPair*> &res_list, HaploPattern *target_pattern = NULL);

	double getLikelihood(const Haplotype &haplotype);
	double getLikelihood(const Genotype &genotype);

protected:
	void clear();
	void clearHaploPairs();
	void initHeadList(const Genotype &genotype);

	void extendAll(int i, Allele a1, Allele a2);
	void extend(HaploPair *hp, Allele a1, Allele a2);
	void calcBackwardLikelihood();

	void findHaploPatternByFreq(double min_freq, int min_len = 1, int max_len = -1);
	void findHaploPatternByNum(int max_num, int min_len = 1, int max_len = -1);
	void findHaploPatternBlockByFreq(double min_freq, int min_len = 1, int max_len = -1);
	void findHaploPatternBlockByNum(int max_num, int min_len = 1, int max_len = -1);

	void generateHaploPatternCandidate(vector<HaploPattern*> &candidates);
	void searchHaploPattern(vector<HaploPattern*> &candidates, double min_freq, int min_num, vector<int> min_len, vector<int> max_len);
	void searchHaploPatternBlock(int min_len, int max_len);

	friend class HaploData;
};


#endif // __HAPLOBUILDER_H
