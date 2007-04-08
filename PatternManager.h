
#ifndef __PATTERNMANAGER_H
#define __PATTERNMANAGER_H


#include <vector>
#include <list>

#include "HaploPattern.h"
#include "PatternTree.h"


class Allele;
class Genotype;
class HaploBuilder;


class PatternManager {
	typedef list<pair<int, double> > MatchingState;

	struct PatternCandidate {
		HaploPattern *pattern;
		MatchingState state;

		PatternCandidate(const GenoData *genos, int start = 0)
			: pattern(new HaploPattern(*genos, start)) { }
		~PatternCandidate() { delete pattern; }
		HaploPattern *release() {
			HaploPattern *hp = pattern;
			pattern = 0;
			state.clear();
			return hp; 
		}
	};

	HaploBuilder &m_builder;
	vector<HaploPattern*> m_patterns;
	vector<PatternCandidate*> m_candidates;
	double m_min_freq;
	vector<int> m_min_len;
	vector<int> m_max_len;

	tr1::shared_ptr<BackwardPatternTree> m_pattern_tree;
	vector<HaploPattern*> m_head_list;

public:
	PatternManager(HaploBuilder &hb) : m_builder(hb) { }
	~PatternManager();

	HaploPattern *operator [](int i) { return m_patterns[i]; }
	const HaploPattern *operator [](int i) const { return m_patterns[i]; }
	const BackwardPatternTree *pattern_tree() const { return m_pattern_tree.get(); }
	const vector<HaploPattern*> &head_list() const { return m_head_list; }
	int head_len() const { return m_min_len[0]; }

	int size() const { return m_patterns.size(); }

	HaploPattern *getSingleAllelePattern(int end, int index) const;
	HaploPattern *getSingleAllelePattern(int end, Allele a) const;

	void findPatternByFreq(double min_freq, int min_len = 2, int max_len = 0);
	void findPatternByNum(int max_num, int min_len = 2, int max_len = 0);
	void findPatternBlock(int len);

	void adjustFrequency();
	void estimateFrequency();
	void estimatePatterns();

protected:
	void generateCandidates();
	void searchPattern(bool reserve_candidates = false);

	void checkFrequency(HaploPattern *hp, MatchingState &ms) const;
	void checkFrequencyWithExtension(HaploPattern *hp, MatchingState &ms, const MatchingState &old_ms, int start, int len = 1) const;
	double getMatchingFrequency(const Genotype &g, const Allele *pa, int start, int len) const;

	void initialize();
	void extendPatterns(vector<HaploPattern*> &patterns, vector<HaploPattern*> &seeds);
};


#endif // __PATTERNMANAGER_H
