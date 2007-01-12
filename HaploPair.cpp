
#include "HaploPair.h"
#include "HaploPattern.h"
#include "HaploData.h"

#include "MemLeak.h"


boost::pool<> PatternPair::m_pool(sizeof(PatternPair));
boost::pool<> HaploPair::m_pool(sizeof(HaploPair));


HaploPair::HaploPair(const HaploPattern *hpa, const HaploPattern *hpb, HaploPattern *target_pattern)
: m_pair(new PatternPair(hpa, hpb)),
  m_end(hpa->end())
{
	if (hpa->end() != hpb->end()) {
		Logger::error("Inconsistent HaploPattern!");
		exit(1);
	}
	getID(m_id, hpa, hpb);
	m_likelihood = hpa->frequency() * hpb->frequency();
	m_total_likelihood = m_likelihood;
	m_homogenous = false;
	if (hpa->id() == hpb->id()) {
		m_homogenous = true;
	}
	else if (hpa->id() > hpb->id()) {
		Logger::warning("HaploPair maybe repeatedly counted!");
	}
	m_half = false;
	m_match_a = m_match_b = true;
	m_match_next_a = m_match_next_b = true;
	if (target_pattern != NULL && end() > target_pattern->start()) {
		m_match_a = hpa->isMatch(*target_pattern);
		m_match_b = hpb->isMatch(*target_pattern);
		if (end() >= target_pattern->end()) {
			if (!m_match_a || !m_match_b){
				m_likelihood *= 0.5;
				m_total_likelihood *= 0.5;
			}
		}
	}
}

HaploPair::HaploPair(HaploPair *hp, const HaploPattern *hpa, const HaploPattern *hpb)
: m_pair(new PatternPair(hpa, hpb, hp->m_pair)),
  m_end(hp->end()+1),
  m_likelihood(hp->m_likelihood * hpa->transition_prob() * hpb->transition_prob()),
  m_total_likelihood(hp->m_total_likelihood * hpa->transition_prob() * hpb->transition_prob()),
  m_homogenous(hp->m_homogenous && hpa->id() == hpb->id()),
  m_half(hp->m_half),
  m_match_a(hp->m_match_next_a),
  m_match_b(hp->m_match_next_b)
{
	if (m_half) {
		m_likelihood *= 0.5;
		m_total_likelihood *= 0.5;
		m_half = false;
	}
	getID(m_id, hpa, hpb);
}

Genotype HaploPair::getGenotype()
{
	int i;
	Genotype g(m_pair->m_pattern_a.haplodata().genotype_len());
	SP_PatternPair pp;
	i = end() - 1;
	pp = m_pair;
	while (pp) {
		if (pp->m_prev) {
			g(0)[i] = pp->m_pattern_a[i-pp->m_pattern_a.start()];
			g(1)[i] = pp->m_pattern_b[i-pp->m_pattern_b.start()];
		}
		else {
			copy(&pp->m_pattern_a[0], &pp->m_pattern_a[i+1-pp->m_pattern_b.start()], &g(0)[pp->m_pattern_a.start()]);
			copy(&pp->m_pattern_b[0], &pp->m_pattern_b[i+1-pp->m_pattern_b.start()], &g(1)[pp->m_pattern_b.start()]);
			break;
		}
		pp = pp->m_prev;
		--i;
	}
	g.checkGenotype();
	return g;
}

bool HaploPair::extendable(int a, int b, HaploPattern *target_pattern)
{
	const HaploPattern *hpa, *hpb;
	if (!m_match_a && !m_match_b) {
		return false;
	}
	m_half = false;
	m_match_next_a = m_match_a;
	m_match_next_b = m_match_b;
	hpa = successor_a(a);
	hpb = successor_b(b);
	if (hpa != NULL && hpb != NULL) {
		if (m_homogenous && hpa->id() > hpb->id()) {
			return false;
		}
		if (target_pattern != NULL && end() >= target_pattern->start() && end() < target_pattern->end()) {
			m_match_next_a = m_match_a && hpa->isMatch(*target_pattern, end());
			m_match_next_b = m_match_b && hpb->isMatch(*target_pattern, end());
			if (!m_match_next_a && !m_match_next_b) {
				return false;
			}
			if (end()+1 == target_pattern->end()) {
				if (!m_match_next_a || !m_match_next_b){
					m_half = true;
				}
			}
		}
		return true;
	}
	return false;
}

void HaploPair::extend_trial(int a, int b, int id[2], double &likelihood, double &total_likelihood)
{
	const HaploPattern *hpa, *hpb;
	hpa = successor_a(a);
	hpb = successor_b(b);
	getID(id, hpa, hpb);
	likelihood = m_likelihood * hpa->transition_prob() * hpb->transition_prob();
	total_likelihood = m_total_likelihood * hpa->transition_prob() * hpb->transition_prob();
	if (m_half) {
		likelihood *= 0.5;
		total_likelihood *= 0.5;
	}
}

void HaploPair::getID(int id[2], const HaploPattern *hpa, const HaploPattern *hpb)
{
	if (hpa == NULL || hpb == NULL) {
		id[0] = -1;
		id[1] = -1;
		Logger::warning("Invalid HaploPair ID!");
	}
	else if (hpa->start() < hpb->start()) {
		id[0] = hpa->id();
		id[1] = hpb->id();
	}
	else if (hpa->start() > hpb->start()) {
		id[0] = hpb->id();
		id[1] = hpa->id();
	}
	else if (hpa->id() < hpb->id()) {
		id[0] = hpa->id();
		id[1] = hpb->id();
	}
	else {
		id[0] = hpb->id();
		id[1] = hpa->id();
	}
}
