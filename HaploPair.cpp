
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

double HaploPair::evaluatePattern(const HaploPattern *pattern, vector<HaploPair*> &hp_list)
{
	int i;
	char match;
	HaploPair *hp, *last_hp;
	vector<pair<HaploPair*, double> > last_list[3], new_list[3];
	vector<HaploPair*>::iterator i_hp;
	vector<pair<HaploPair*, double> >::iterator i_mp;
	vector<pair<HaploPair*, double> >::const_iterator i_link;

	for (i_hp=hp_list.begin(); i_hp!=hp_list.end(); ++i_hp) {
		hp = *i_hp;
		match = 0;
		if (hp->allele_a() == (*pattern)[0]) {
			match |= 1;
		}
		if (hp->allele_b() == (*pattern)[0]) {
			match |= 2;
		}
		switch (match) {
		case 1:
			last_list[1].push_back(make_pair(hp, hp->forward_likelihood() / 2.0));
			break;
		case 2:
			last_list[2].push_back(make_pair(hp, hp->forward_likelihood() / 2.0));
			break;
		case 3:
			last_list[0].push_back(make_pair(hp, hp->forward_likelihood()));
			break;
		}
	}

	Logger::resumeTimer(3);
	for (i=1; i<pattern->length(); ++i) {
		for (i_mp=last_list[0].begin(); i_mp!=last_list[0].end(); ++i_mp) {
			last_hp = (*i_mp).first;
			for (i_link=last_hp->m_forward_links.begin(); i_link!=last_hp->m_forward_links.end(); ++i_link) {
				hp = (*i_link).first;
				match = 0;
				if (hp->allele_a() == (*pattern)[i]) {
					match |= 1;
				}
				if (hp->allele_b() == (*pattern)[i]) {
					match |= 2;
				}
				switch (match) {
				case 1:
					new_list[1].push_back(make_pair(hp, (*i_mp).second * (*i_link).second / 2.0));
					break;
				case 2:
					new_list[2].push_back(make_pair(hp, (*i_mp).second * (*i_link).second / 2.0));
					break;
				case 3:
					new_list[0].push_back(make_pair(hp, (*i_mp).second * (*i_link).second));
					break;
				}
			}
		}

		for (i_mp=last_list[1].begin(); i_mp!=last_list[1].end(); ++i_mp) {
			last_hp = (*i_mp).first;
			for (i_link=last_hp->m_forward_links.begin(); i_link!=last_hp->m_forward_links.end(); ++i_link) {
				hp = (*i_link).first;
				if (hp->allele_a() == (*pattern)[i]) {
					new_list[1].push_back(make_pair(hp, (*i_mp).second * (*i_link).second));
				}
			}
		}

		for (i_mp=last_list[2].begin(); i_mp!=last_list[2].end(); ++i_mp) {
			last_hp = (*i_mp).first;
			for (i_link=last_hp->m_forward_links.begin(); i_link!=last_hp->m_forward_links.end(); ++i_link) {
				hp = (*i_link).first;
				if (hp->allele_b() == (*pattern)[i]) {
					new_list[2].push_back(make_pair(hp, (*i_mp).second * (*i_link).second));
				}
			}
		}

		new_list[0].swap(last_list[0]);
		new_list[1].swap(last_list[1]);
		new_list[2].swap(last_list[2]);
		new_list[0].clear();
		new_list[1].clear();
		new_list[2].clear();
	}
	Logger::pauseTimer(3);

	double freq = 0;
	for (i_mp=last_list[0].begin(); i_mp!=last_list[0].end(); ++i_mp) {
		hp = (*i_mp).first;
		freq += (*i_mp).second * hp->backward_likelihood();
	}
	for (i_mp=last_list[1].begin(); i_mp!=last_list[1].end(); ++i_mp) {
		hp = (*i_mp).first;
		freq += (*i_mp).second * hp->backward_likelihood();
	}
	for (i_mp=last_list[2].begin(); i_mp!=last_list[2].end(); ++i_mp) {
		hp = (*i_mp).first;
		freq += (*i_mp).second * hp->backward_likelihood();
	}

	return freq;
}
