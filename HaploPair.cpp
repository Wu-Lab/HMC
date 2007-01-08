
#include "HaploPair.h"


HaploPair::HaploPair(HaploPattern *hpa, HaploPattern *hpb, HaploPattern *target_pattern)
: m_pair(new PatternPair(hpa, hpb)),
  m_end(hpa->m_end)
{
	if (hpa->m_end != hpb->m_end) {
		Logger::error("Inconsistent HaploPattern!");
		exit(1);
	}
	getID(m_id, hpa, hpb);
	m_likelihood = hpa->m_frequency * hpb->m_frequency;
	m_total_likelihood = m_likelihood;
	m_homogenous = false;
	if (hpa->m_id == hpb->m_id) {
		m_homogenous = true;
	}
	else if (hpa->m_id > hpb->m_id) {
		Logger::warning("HaploPair maybe repeatedly counted!");
	}
	m_half = false;
	m_match_a = m_match_b = true;
	m_match_next_a = m_match_next_b = true;
	if (target_pattern != NULL && m_end > target_pattern->m_start) {
		m_match_a = hpa->isMatch(*target_pattern);
		m_match_b = hpb->isMatch(*target_pattern);
		if (m_end >= target_pattern->m_end) {
			if (!m_match_a || !m_match_b){
				m_likelihood *= 0.5;
				m_total_likelihood *= 0.5;
			}
		}
	}
}

Genotype HaploPair::getGenotype()
{
	int i;
	Genotype g(m_pair->pattern_a.m_haplo_data->m_genotype_len);
	tr1::shared_ptr<PatternPair> pp;
	i = m_end - 1;
	pp = m_pair;
	while (pp) {
		if (pp->prev) {
			g.haplotypes(0).allele(i) = pp->pattern_a[i-pp->pattern_a.start()];
			g.haplotypes(1).allele(i) = pp->pattern_b[i-pp->pattern_b.start()];
		}
		else {
			while (i >= 0) {
				g.haplotypes(0).allele(i) = pp->pattern_a[i-pp->pattern_a.start()];
				g.haplotypes(1).allele(i) = pp->pattern_b[i-pp->pattern_b.start()];
				i--;
			}
			break;
		}
		pp = pp->prev;
		i--;
	}
	g.checkGenotype();
	return g;
}

bool HaploPair::extendable(int a, int b, HaploPattern *target_pattern)
{
	HaploPattern *hpa, *hpb;
	if (!m_match_a && !m_match_b) {
		return false;
	}
	m_half = false;
	m_match_next_a = m_match_a;
	m_match_next_b = m_match_b;
	hpa = m_pair->pattern_a.m_successors[a];
	hpb = m_pair->pattern_b.m_successors[b];
	if (hpa != NULL && hpb != NULL) {
		if (m_homogenous && hpa->m_id > hpb->m_id) {
			return false;
		}
		if (target_pattern != NULL && m_end >= target_pattern->m_start && m_end < target_pattern->m_end) {
			m_match_next_a = m_match_a && hpa->isMatch(*target_pattern, m_end);
			m_match_next_b = m_match_b && hpb->isMatch(*target_pattern, m_end);
			if (!m_match_next_a && !m_match_next_b) {
				return false;
			}
			if (m_end+1 == target_pattern->m_end) {
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
	HaploPattern *hpa, *hpb;
	hpa = m_pair->pattern_a.m_successors[a];
	hpb = m_pair->pattern_b.m_successors[b];
	getID(id, hpa, hpb);
	likelihood = m_likelihood * hpa->m_transition_prob * hpb->m_transition_prob;
	total_likelihood = m_total_likelihood * hpa->m_transition_prob * hpb->m_transition_prob;
	if (m_half) {
		likelihood *= 0.5;
		total_likelihood *= 0.5;
	}
}

void HaploPair::extend(int a, int b)
{
	HaploPattern *hpa, *hpb;
	hpa = m_pair->pattern_a.m_successors[a];
	hpb = m_pair->pattern_b.m_successors[b];
	tr1::shared_ptr<PatternPair> pair(new PatternPair(hpa, hpb));
	m_pair.swap(pair->prev);
	m_pair.swap(pair);
	m_likelihood *= hpa->m_transition_prob;
	m_likelihood *= hpb->m_transition_prob;
	m_total_likelihood *= hpa->m_transition_prob;
	m_total_likelihood *= hpb->m_transition_prob;
	if (m_homogenous && hpa->m_id != hpb->m_id) {
		m_homogenous = false;
	}
	if (m_half) {
		m_likelihood *= 0.5;
		m_total_likelihood *= 0.5;
	}
	m_half = false;
	m_match_a = m_match_next_a;
	m_match_b = m_match_next_b;
	getID(m_id, hpa, hpb);
	m_end++;
}

void HaploPair::getID(int id[2], HaploPattern *hpa, HaploPattern *hpb)
{
	if (hpa == NULL || hpb == NULL) {
		id[0] = -1;
		id[1] = -1;
		Logger::warning("Invalid HaploPair ID!");
	}
	else if (hpa->m_start < hpb->m_start) {
		id[0] = hpa->m_id;
		id[1] = hpb->m_id;
	}
	else if (hpa->m_start > hpb->m_start) {
		id[0] = hpb->m_id;
		id[1] = hpa->m_id;
	}
	else if (hpa->m_id < hpb->m_id) {
		id[0] = hpa->m_id;
		id[1] = hpb->m_id;
	}
	else {
		id[0] = hpb->m_id;
		id[1] = hpa->m_id;
	}
}
