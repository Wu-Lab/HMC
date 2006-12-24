
#include "HaploPair.h"


HaploPair::HaploPair(HaploPattern *hpa, HaploPattern *hpb)
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
	m_likelihood = hpa->m_frequency * hpb->m_frequency;
	m_total_likelihood = m_likelihood;
	m_homogenous = false;
	if (hpa->m_id == hpb->m_id) {
		m_homogenous = true;
	}
	else if (hpa->m_id > hpb->m_id) {
		Logger::warning("HaploPair maybe repeatedly counted!");
	}
}

HaploPair::HaploPair(const HaploPair &hp)
{
	int i;
	m_end = hp.m_end;
	m_genotype_len = hp.m_genotype_len;
	m_patterns[0] = new HaploPattern * [m_genotype_len+1];
	m_patterns[1] = new HaploPattern * [m_genotype_len+1];
	for (i=0; i<=m_end; i++) {
		m_patterns[0][i] = hp.m_patterns[0][i];
		m_patterns[1][i] = hp.m_patterns[1][i];
	}
	m_weight = hp.m_weight;
	m_likelihood = hp.m_likelihood;
	m_total_likelihood = hp.m_total_likelihood;
	m_homogenous = hp.m_homogenous;
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
	g.setWeight(m_weight);
	g.setLength(m_genotype_len);
	g.checkGenotype();
	return g;
}

void HaploPair::getID(int id[2])
{
	if (m_patterns[0][m_end] == NULL || m_patterns[1][m_end] == NULL) {
		id[0] = -1;
		id[1] = -1;
		Logger::warning("Invalid HaploPair!");
	}
	else if (m_patterns[0][m_end]->m_start < m_patterns[1][m_end]->m_start) {
		id[0] = m_patterns[0][m_end]->m_id;
		id[1] = m_patterns[1][m_end]->m_id;
	}
	else if (m_patterns[0][m_end]->m_start > m_patterns[1][m_end]->m_start) {
		id[0] = m_patterns[1][m_end]->m_id;
		id[1] = m_patterns[0][m_end]->m_id;
	}
	else if (m_patterns[0][m_end]->m_id < m_patterns[1][m_end]->m_id) {
		id[0] = m_patterns[0][m_end]->m_id;
		id[1] = m_patterns[1][m_end]->m_id;
	}
	else {
		id[0] = m_patterns[1][m_end]->m_id;
		id[1] = m_patterns[0][m_end]->m_id;
	}
}

bool HaploPair::extendable(int a, int b)
{
	if (m_patterns[0][m_end]->m_successors[a] != NULL && m_patterns[1][m_end]->m_successors[b] != NULL) {
		if (m_homogenous && m_patterns[0][m_end]->m_successors[a]->m_id > m_patterns[1][m_end]->m_successors[b]->m_id) {
			return false;
		}
		else {
			return true;
		}
	}
	else {
		return false;
	}
}

void HaploPair::extend(int a, int b, HaploPattern *hpa, HaploPattern *hpb)
{
	if (hpa != NULL) {
		m_patterns[0][m_end+1] = hpa;
		m_likelihood *= 0.05;
		m_total_likelihood *= 0.05;
	}
	else {
		m_patterns[0][m_end+1] = m_patterns[0][m_end]->m_successors[a];
	}
	if (hpb != NULL) {
		m_patterns[1][m_end+1] = hpb;
		m_likelihood *= 0.05;
		m_total_likelihood *= 0.05;
	}
	else {
		m_patterns[1][m_end+1] = m_patterns[1][m_end]->m_successors[b];
	}
	m_end++;
	m_likelihood *= m_patterns[0][m_end]->m_transition_prob;
	m_likelihood *= m_patterns[1][m_end]->m_transition_prob;
	m_total_likelihood *= m_patterns[0][m_end]->m_transition_prob;
	m_total_likelihood *= m_patterns[1][m_end]->m_transition_prob;
	if (m_homogenous && m_patterns[0][m_end]->m_id != m_patterns[1][m_end]->m_id) {
		m_homogenous = false;
	}
}
