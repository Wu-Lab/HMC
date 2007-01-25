
#include "HaploPair.h"
#include "HaploPattern.h"
#include "HaploData.h"

#include "MemLeak.h"


boost::pool<> HaploPair::m_pool(sizeof(HaploPair));


HaploPair::HaploPair(const HaploPattern *hpa, const HaploPattern *hpb)
: m_pattern_a(*hpa), m_pattern_b(*hpb),
  m_backward_link(0), m_backward_likelihood(1.0)
{
	if (hpa->end() != hpb->end()) {
		Logger::error("Construct HaplPair from inconsistent HaploPatterns!");
		exit(1);
	}
	if (hpa->start() != 0 || hpb->start() != 0) {
		Logger::error("Construct HaplPair from middle!");
		exit(1);
	}
	m_forward_likelihood = m_best_likelihood
		= m_transition_prob = hpa->frequency() * hpb->frequency();
}

HaploPair::HaploPair(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb)
: m_pattern_a(*hpa), m_pattern_b(*hpb),
  m_backward_link(hp), m_backward_likelihood(1.0)
{
	m_transition_prob = hpa->transition_prob() * hpb->transition_prob();
	m_best_likelihood = hp->m_best_likelihood * m_transition_prob;
	m_forward_likelihood = hp->m_forward_likelihood * m_transition_prob;
	hp->m_forward_links.push_back(make_pair(this, m_transition_prob));
}

void HaploPair::add(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb)
{
	double likelihood = hp->m_best_likelihood * m_transition_prob;
	if (m_best_likelihood < likelihood) {
		m_best_likelihood = likelihood;
		m_backward_link = hp;
	}
	m_forward_likelihood += hp->m_forward_likelihood * m_transition_prob;
	hp->m_forward_links.push_back(make_pair(this, m_transition_prob));
}

Genotype HaploPair::getGenotype()
{
	int i;
	Genotype g(m_pattern_a.haplodata().genotype_len());
	HaploPair *hp;
	i = end() - 1;
	hp = this;
	while (hp) {
		if (hp->m_backward_link) {
			g(0)[i] = hp->m_pattern_a[i-hp->m_pattern_a.start()];
			g(1)[i] = hp->m_pattern_b[i-hp->m_pattern_b.start()];
			hp = hp->m_backward_link;
			--i;
		}
		else {
			copy(&hp->m_pattern_a[0], &hp->m_pattern_a[i+1-hp->m_pattern_b.start()], &g(0)[hp->m_pattern_a.start()]);
			copy(&hp->m_pattern_b[0], &hp->m_pattern_b[i+1-hp->m_pattern_b.start()], &g(1)[hp->m_pattern_b.start()]);
			break;
		}
	}
	g.checkGenotype();
	return g;
}

struct MatchPair {
	const HaploPair *hp;
	char match;
	double likelihood;
	MatchPair(const HaploPair *h, char m, double l) : hp(h), match(m), likelihood(l) { }
};

double HaploPair::evaluatePattern(const HaploPattern *pattern, vector<HaploPair*> &hp_list)
{
	int i;
	char match;
	double weight;
	const HaploPair *hp, *last_hp;
	vector<MatchPair> last_list, new_list;
	vector<HaploPair*>::iterator i_hp;
	vector<MatchPair>::iterator i_mp;
	vector<pair<HaploPair*, double> >::const_iterator i_link;

	for (i_hp=hp_list.begin(); i_hp!=hp_list.end(); ++i_hp) {
		hp = *i_hp;
		match = 0;
		if (hp->allele_a() == (*pattern)[0]) {
			match += 1;
		}
		if (hp->allele_b() == (*pattern)[0]) {
			match += 2;
		}
		if (match > 0) {
			if (match == 3) {
				weight = hp->forward_likelihood();
			}
			else {
				weight = hp->forward_likelihood() / 2;
			}
			last_list.push_back(MatchPair(hp, match, weight));
		}
	}

	Logger::resumeTimer(3);
	for (i=1; i<pattern->length(); ++i) {
		for (i_mp=last_list.begin(); i_mp!=last_list.end(); ++i_mp) {
			last_hp = (*i_mp).hp;
			for (i_link=last_hp->m_forward_links.begin(); i_link!=last_hp->m_forward_links.end(); ++i_link) {
				hp = (*i_link).first;
				match = (*i_mp).match;
				weight = 1.0;
				if ((match & 1) && hp->allele_a() != (*pattern)[i]) {
					match -= 1;
					weight = 0.5;
				}
				if ((match & 2) && hp->allele_b() != (*pattern)[i]) {
					match -= 2;
					weight = 0.5;
				}
				if (match > 0) {
					weight *= (*i_mp).likelihood * (*i_link).second;
					new_list.push_back(MatchPair(hp, match, weight));
				}
			}
		}
		new_list.swap(last_list);
		new_list.clear();
	}
	Logger::pauseTimer(3);

	double freq = 0;
	for (i_mp=last_list.begin(); i_mp!=last_list.end(); ++i_mp) {
		hp = (*i_mp).hp;
		freq += (*i_mp).likelihood * hp->backward_likelihood();
	}

	return freq;
}
