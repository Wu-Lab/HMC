
#include <map>

#include "HaploPair.h"
#include "HaploPattern.h"
#include "GenoData.h"

#include "MemLeak.h"


boost::pool<> HaploPair::m_pool(sizeof(HaploPair));


HaploPair::HaploPair(const HaploPattern *hpa, const HaploPattern *hpb)
: m_pattern_a(*hpa), m_pattern_b(*hpb),
  m_allele_a((*hpa)[hpa->length()-1]), m_allele_b((*hpb)[hpb->length()-1]),
  m_backward_likelihood(1.0)
{
	if (hpa->end() != hpb->end()) {
		Logger::error("Construct HaplPair from inconsistent HaploPatterns (end %d, %d) !", hpa->end(), hpb->end());
		exit(1);
	}
	if (hpa->start() != 0 || hpb->start() != 0) {
		Logger::error("Construct HaplPair from middle!");
		exit(1);
	}
	m_forward_likelihood = m_transition_prob = hpa->frequency() * hpb->frequency();
	if (hpa->id() != hpb->id()) {
		m_forward_likelihood *= 2.0;
	}
	m_best_links.resize(1);
	m_best_links[0].set(0, 0, false, m_transition_prob);
}

HaploPair::HaploPair(const HaploPattern *hpa, const HaploPattern *hpb, HaploPair *hp, bool reversed, double prob)
: m_pattern_a(*hpa), m_pattern_b(*hpb),
  m_allele_a((*hpa)[hpa->length()-1]), m_allele_b((*hpb)[hpb->length()-1]),
  m_backward_likelihood(1.0)
{
	int i, n;
	m_transition_prob = hpa->transition_prob() * hpb->transition_prob() * prob;
	m_forward_likelihood = hp->m_forward_likelihood * m_transition_prob;
	hp->m_forward_links[reversed ? 1 : 0].push_back(this);
	n = hp->m_best_links.size();
	m_best_links.resize(n);
	for (i=0; i<n; ++i) {
		m_best_links[i].set(hp, i, reversed, hp->m_best_links[i].likelihood * m_transition_prob);
	}
}

void HaploPair::add(HaploPair *hp, bool reversed, int best_num)
{
	int i, m, n;
	double likelihood;
	vector<HaploPairLink>::iterator i_link;
	m_forward_likelihood += hp->m_forward_likelihood * m_transition_prob;
	hp->m_forward_links[reversed ? 1 : 0].push_back(this);
	m = max(best_num-m_best_links.size(), 0);
	n = hp->m_best_links.size();
	m = min(m, n);
	for (i=0; i<m; ++i) {
		m_best_links.push_back(HaploPairLink(hp, i, reversed, hp->m_best_links[i].likelihood * m_transition_prob));
	}
	for (i=m; i<n; ++i) {
		i_link = min_element(m_best_links.begin(), m_best_links.end());
		likelihood = hp->m_best_links[i].likelihood * m_transition_prob;
		if (likelihood > i_link->likelihood) {
			i_link->set(hp, i, reversed, likelihood);
		}
		else {
			break;
		}
	}
}

Genotype HaploPair::getGenotype(int index) const
{
	int a = 0, b = 1;
	int i = end() - 1;
	Genotype g(m_pattern_a.genos().genotype_len());
	g.setPriorProbability(getLikelihood(index));
	const HaploPair *hp = this;
	while (hp) {
		if ((index < hp->m_best_links.size()) && hp->m_best_links[index].link) {
			g(a)[i] = hp->m_pattern_a[i-hp->m_pattern_a.start()];
			g(b)[i] = hp->m_pattern_b[i-hp->m_pattern_b.start()];
			if (hp->m_best_links[index].reversed) swap(a, b);
			hp = hp->m_best_links[index].link;
			index = hp->m_best_links[index].index;
			--i;
		}
		else {
			const Allele *first = &hp->m_pattern_a[0];
			copy(first, first+i+1-hp->m_pattern_a.start(), &g(a)[hp->m_pattern_a.start()]);
			first = &hp->m_pattern_b[0];
			copy(first, first+i+1-hp->m_pattern_b.start(), &g(b)[hp->m_pattern_b.start()]);
			break;
		}
	}
	g.checkGenotype();
	return g;
}

void HaploPair::calcBackwardLikelihood()
{
	vector<HaploPair*>::const_iterator i_link;
	m_backward_likelihood = 0;
	for (int i=0; i<2; ++i) {
		for (i_link=m_forward_links[i].begin(); i_link!=m_forward_links[i].end(); ++i_link) {
			HaploPair *next_hp = *i_link;
			m_backward_likelihood += next_hp->backward_likelihood() * next_hp->transition_prob();
		}
	}
}
