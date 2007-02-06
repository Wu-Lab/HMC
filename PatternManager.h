
#ifndef __PATTERNMANAGER_H
#define __PATTERNMANAGER_H


#include <vector>
#include <list>

#include "HaploPattern.h"


class Allele;
class Genotype;
class HaploBuilder;


class PatternManager {
	typedef list<pair<int, double> > MatchingState;

	struct PatternCandidate {
		HaploPattern *pattern;
		MatchingState state;

		PatternCandidate(const HaploData *hd, int start = 0)
			: pattern(new HaploPattern(*hd, start)) { }
		~PatternCandidate() { delete pattern; }
		HaploPattern *release() {
			HaploPattern *hp = pattern;
			pattern = 0;
			state.clear();
			return hp; 
		}
	};

	const HaploBuilder &m_builder;
	vector<HaploPattern*> m_patterns;
	vector<PatternCandidate*> m_candidates;
	vector<int> m_min_len;
	vector<int> m_max_len;

public:
	PatternManager(const HaploBuilder &hb) : m_builder(hb) { }
	~PatternManager();

	HaploPattern *operator [](int i) { return m_patterns[i]; }
	const HaploPattern *operator [](int i) const { return m_patterns[i]; }

	int size() const { return m_patterns.size(); }

	void findPatternByFreq(double min_freq, int min_len = 2, int max_len = 0);
	void findPatternByNum(int max_num, int min_len = 2, int max_len = 0);
	void findPatternBlock(int len);

protected:
	void generateCandidates();
	void searchPattern(double min_freq, bool reserve_candidates = false);

	void checkFrequency(HaploPattern *hp, MatchingState &ms) const;
	void checkFrequencyWithExtension(HaploPattern *hp, MatchingState &ms, const MatchingState &old_ms, int start, int len = 1) const;
	double getMatchingFrequency(const Genotype &g, const Allele *pa, int start, int len) const;
};


#endif // __PATTERNMANAGER_H
