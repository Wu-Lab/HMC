
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
	bool homo = (hpa->id() == hpb->id());
	if (!homo) {
		m_forward_likelihood *= 2.0;
	}
	m_best_links.push_back(HaploPairLink(0, 0, false, homo, m_transition_prob));
}

HaploPair::HaploPair(const HaploPattern *hpa, const HaploPattern *hpb, HaploPair *hp, bool reversed)
: m_pattern_a(*hpa), m_pattern_b(*hpb),
  m_allele_a((*hpa)[hpa->length()-1]), m_allele_b((*hpb)[hpb->length()-1]),
  m_backward_likelihood(1.0)
{
	int i, n;
	double likelihood;
	m_transition_prob = hpa->transition_prob() * hpb->transition_prob();
	m_forward_likelihood = hp->m_forward_likelihood * m_transition_prob;
	hp->m_forward_links[reversed ? 1 : 0].push_back(this);
	n = hp->m_best_links.size();
	m_best_links = hp->m_best_links;
	for (i=0; i<n; ++i) {
		m_best_links[i].link = hp;
		m_best_links[i].index = i;
		m_best_links[i].reversed = reversed;
		m_best_links[i].likelihood *= m_transition_prob;
	}
	if (m_allele_a != m_allele_b) {
		for (i=0; i<n; ++i) {
			if (m_best_links[i].homozygous) {
				if (reversed) m_best_links[i].likelihood = 0;
				m_best_links[i].homozygous = false;
			}
		}
	}
}

void HaploPair::add(HaploPair *hp, bool reversed, int best_num)
{
	int i, k, n;
	m_forward_likelihood += hp->m_forward_likelihood * m_transition_prob;
	hp->m_forward_links[reversed ? 1 : 0].push_back(this);
	k = m_best_links.size();
	n = hp->m_best_links.size();
	m_best_links.insert(m_best_links.end(), hp->m_best_links.begin(), hp->m_best_links.end());
	for (i=k; i<k+n; ++i) {
		m_best_links[i].link = hp;
		m_best_links[i].index = i-k;
		m_best_links[i].reversed = reversed;
		m_best_links[i].likelihood *= m_transition_prob;
	}
	if (m_allele_a != m_allele_b) {
		for (i=k; i<k+n; ++i) {
			if (m_best_links[i].homozygous) {
				if (reversed) m_best_links[i].likelihood = 0;
				m_best_links[i].homozygous = false;
			}
		}
	}
	if (m_best_links.size() > best_num) {
		nth_element(m_best_links.begin(), m_best_links.begin()+best_num-1, m_best_links.end(), greater<HaploPairLink>());
		m_best_links.resize(best_num);
	}
}

Genotype HaploPair::getGenotype(int index) const
{
	const HaploPair *hp, *next_hp;
	int a = 0, b = 1;
	int i = end() - 1;
	Genotype g(m_pattern_a.genos().genotype_len());
	if (m_best_links[index].homozygous) {
		g.setPriorProbability(getLikelihood(index));
	}
	else {
		g.setPriorProbability(getLikelihood(index) * 2.0);
	}
	hp = this;
	while (hp) {
		if ((index < hp->m_best_links.size()) && hp->m_best_links[index].link) {
			g(a)[i] = hp->m_pattern_a[i-hp->m_pattern_a.start()];
			g(b)[i] = hp->m_pattern_b[i-hp->m_pattern_b.start()];
			if (hp->m_best_links[index].reversed) swap(a, b);
			next_hp = hp->m_best_links[index].link;
			index = hp->m_best_links[index].index;
			hp = next_hp;
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
