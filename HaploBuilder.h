
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
#include "PatternManager.h"


class HaploBuilder {
protected:
	HaploData *m_haplodata;
	PatternManager m_patterns;

	vector<vector<HaploPair*> > m_haplopairs;
	vector<map<int, int> > m_best_pair;

	double m_current_genotype_probability;

public:
	HaploBuilder();
	~HaploBuilder();

	HaploData *haplodata() { return m_haplodata; }
	const HaploData *haplodata() const { return m_haplodata; }
	const HaploPattern *patterns(int i) const { return m_patterns[i]; }

	int pattern_num() const { return m_patterns.size(); }
	int genotype_num() const { return m_haplodata->genotype_num(); }
	int genotype_len() const { return m_haplodata->genotype_len(); }

	void setHaploData(HaploData &hd);

	void initialize();

	void resolve(const Genotype &genotype, Genotype &resolution, vector<HaploPair*> &res_list, HaploPattern *target_pattern = 0);

	double getLikelihood(const Haplotype &haplotype);
	double getLikelihood(const Genotype &genotype);

	void estimateFrequency(vector<HaploPattern*> &patterns);

protected:
	void clear();
	void clearHaploPairs();
	void initHeadList(const Genotype &genotype);

	void extendAll(int i, Allele a1, Allele a2);
	void extend(HaploPair *hp, Allele a1, Allele a2);

	void calcBackwardLikelihood();
	double estimateFrequency(PatternNode *node, int locus, const Allele &a, double last_freq, const map<HaploPair*, double> last_match[2]);
};


#endif // __HAPLOBUILDER_H
