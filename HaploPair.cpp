
#include <map>

#include "HaploPair.h"
#include "HaploPattern.h"
#include "HaploData.h"

#include "MemLeak.h"


boost::pool<> HaploPair::m_pool(sizeof(HaploPair));


HaploPair::HaploPair(const HaploPattern *hpa, const HaploPattern *hpb)
: m_pattern_a(*hpa), m_pattern_b(*hpb),
  m_allele_a((*hpa)[hpa->length()-1]), m_allele_b((*hpb)[hpb->length()-1]),
  m_backward_link(0), m_backward_likelihood(1.0)
{
	if (hpa->end() != hpb->end()) {
		Logger::error("Construct HaplPair from inconsistent HaploPatterns (end %d, %d) !", hpa->end(), hpb->end());
		exit(1);
	}
	if (hpa->start() != 0 || hpb->start() != 0) {
		Logger::error("Construct HaplPair from middle!");
		exit(1);
	}
	m_forward_likelihood = m_best_likelihood
		= m_transition_prob = hpa->frequency() * hpb->frequency();
}

HaploPair::HaploPair(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb, double prob)
: m_pattern_a(*hpa), m_pattern_b(*hpb),
  m_allele_a((*hpa)[hpa->length()-1]), m_allele_b((*hpb)[hpb->length()-1]),
  m_backward_link(hp), m_backward_likelihood(1.0)
{
	m_transition_prob = hpa->transition_prob() * hpb->transition_prob() * prob;
// 	double prob_recomb = 0.1;
// 	if (hpa->length() == 1) m_transition_prob * prob_recomb;
// 	if (hpb->length() == 1) m_transition_prob * prob_recomb;
	m_best_likelihood = hp->m_best_likelihood * m_transition_prob;
	m_forward_likelihood = hp->m_forward_likelihood * m_transition_prob;
	hp->m_forward_links.push_back(this);
}

void HaploPair::add(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb)
{
	double likelihood = hp->m_best_likelihood * m_transition_prob;
	if (m_best_likelihood < likelihood) {
		m_best_likelihood = likelihood;
		m_backward_link = hp;
	}
	m_forward_likelihood += hp->m_forward_likelihood * m_transition_prob;
	hp->m_forward_links.push_back(this);
}

Genotype HaploPair::getGenotype()
{
	int i = end() - 1;
	Genotype g(m_pattern_a.haplodata().genotype_len());
	HaploPair *hp = this;
	while (hp) {
		if (hp->m_backward_link) {
			g(0)[i] = hp->m_pattern_a[i-hp->m_pattern_a.start()];
			g(1)[i] = hp->m_pattern_b[i-hp->m_pattern_b.start()];
			hp = hp->m_backward_link;
			--i;
		}
		else {
			const Allele *first = &hp->m_pattern_a[0];
			copy(first, first+i+1-hp->m_pattern_a.start(), &g(0)[hp->m_pattern_a.start()]);
			first = &hp->m_pattern_b[0];
			copy(first, first+i+1-hp->m_pattern_b.start(), &g(1)[hp->m_pattern_b.start()]);
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
	for (i_link=m_forward_links.begin(); i_link!=m_forward_links.end(); ++i_link) {
		HaploPair *next_hp = *i_link;
		m_backward_likelihood += next_hp->backward_likelihood() * next_hp->transition_prob();
	}
}
