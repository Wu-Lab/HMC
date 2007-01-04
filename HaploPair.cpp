
#include "HaploPair.h"


HaploPair::HaploPair(HaploPattern *hpa, HaploPattern *hpb, HaploPattern *target_pattern)
{
	int i;
	if (hpa->m_end != hpb->m_end) {
		Logger::error("Inconsistent HaploPattern!");
		exit(1);
	}
	m_end = hpa->m_end;
	m_genotype_len = hpa->m_haplo_data->m_genotype_len;
	m_patterns[0] = new HaploPattern * [m_genotype_len+1];
	m_patterns[1] = new HaploPattern * [m_genotype_len+1];
	for (i=0; i<m_end; i++) {
		m_patterns[0][i] = NULL;
		m_patterns[1][i] = NULL;
	}
	m_patterns[0][m_end] = hpa;
	m_patterns[1][m_end] = hpb;
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

HaploPair::HaploPair(const HaploPair &hp)
{
	int i;
	Logger::resumeTimer(3);
	m_end = hp.m_end;
	m_genotype_len = hp.m_genotype_len;
	m_patterns[0] = new HaploPattern * [m_genotype_len+1];
	m_patterns[1] = new HaploPattern * [m_genotype_len+1];
	for (i=0; i<=m_end; i++) {
		m_patterns[0][i] = hp.m_patterns[0][i];
		m_patterns[1][i] = hp.m_patterns[1][i];
	}
	m_likelihood = hp.m_likelihood;
	m_total_likelihood = hp.m_total_likelihood;
	m_homogenous = hp.m_homogenous;
	m_half = hp.m_half;
	m_match_a = hp.m_match_a;
	m_match_b = hp.m_match_b;
	m_match_next_a = hp.m_match_next_a;
	m_match_next_b = hp.m_match_next_b;
	Logger::pauseTimer(3);
}

HaploPair::~HaploPair()
{
	delete[] m_patterns[0];
	delete[] m_patterns[1];
}

Genotype HaploPair::getGenotype()
{
	int i, j;
	Genotype g;
	for (i=0; i<2; i++) {
		j = 1;
		while (j <= m_end && m_patterns[i][j] == NULL) j++;
		if (j <= m_end) g.haplotypes(i) = (*m_patterns[i][j++]);
		while (j <= m_end && m_patterns[i][j] != NULL) {
			g.haplotypes(i) += (*m_patterns[i][j])[m_patterns[i][j]->m_length-1];
			j++;
		}
	}
	g.setLength(m_genotype_len);
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
	hpa = m_patterns[0][m_end]->m_successors[a];
	hpb = m_patterns[1][m_end]->m_successors[b];
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
	hpa = m_patterns[0][m_end]->m_successors[a];
	hpb = m_patterns[1][m_end]->m_successors[b];
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
	hpa = m_patterns[0][m_end]->m_successors[a];
	hpb = m_patterns[1][m_end]->m_successors[b];
	m_patterns[0][m_end+1] = hpa;
	m_patterns[1][m_end+1] = hpb;
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
