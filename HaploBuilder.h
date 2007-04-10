
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
#include "HaploData.h"
#include "PatternTree.h"
#include "PatternManager.h"


class HaploBuilder {
protected:
	GenoData *m_genos;
	PatternManager m_patterns;

	HaploData m_samples;
	int m_sample_size;

	vector<vector<HaploPair*> > m_haplopairs;
	vector<map<int, int> > m_best_pair;

	double m_current_genotype_probability;

public:
	HaploBuilder();
	~HaploBuilder();

	GenoData *genos() { return m_genos; }
	const GenoData *genos() const { return m_genos; }
	const HaploPattern *patterns(int i) const { return m_patterns[i]; }
	HaploData *samples() { return &m_samples; }
	const HaploData *samples() const { return &m_samples; }

	int pattern_num() const { return m_patterns.size(); }
	int genotype_num() const { return m_genos->genotype_num(); }
	int genotype_len() const { return m_genos->genotype_len(); }

	void setGenoData(GenoData &genos);

	void initialize();

	void resolve(const Genotype &genotype, Genotype &resolution, vector<Genotype> &res_list, int sample_size = 1);

	double getLikelihood(const Haplotype &haplotype);
	double getLikelihood(const Genotype &genotype);

	void estimateFrequency(vector<HaploPattern*> &patterns);

protected:
	void clear();
	void clearHaploPairs();
	void initHeadList(const Genotype &genotype);

	void extendAll(int i, Allele a1, Allele a2);
	void extend(HaploPair *hp, Allele a1, Allele a2);
	void addHaploPair(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb, double prob = 1.0);

	void calcBackwardLikelihood();
	double estimateFrequency(PatternNode *node, int locus, const Allele &a, double last_freq, const map<HaploPair*, double> last_match[3]);
};


#endif // __HAPLOBUILDER_H
