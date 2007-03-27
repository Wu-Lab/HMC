
#include "HaploBuilder.h"
#include "HaploPair.h"
#include "HaploData.h"

#include "MemLeak.h"


HaploBuilder::HaploBuilder()
: m_patterns(*this)
{
}

HaploBuilder::~HaploBuilder()
{
	for_each(m_haplopairs.begin(), m_haplopairs.end(), DeleteAll_Clear());
}

void HaploBuilder::setHaploData(HaploData &hd)
{
	m_haplodata = &hd;
}

void HaploBuilder::initialize()
{
	for_each(m_haplopairs.begin(), m_haplopairs.end(), DeleteAll_Clear());
	m_haplopairs.resize(genotype_len()+1);
	m_best_pair.resize(pattern_num());
}

void HaploBuilder::resolve(const Genotype &genotype, Genotype &resolution, vector<HaploPair*> &res_list, HaploPattern *target_pattern)
{
	int pn = pattern_num();
	int head_len = m_patterns.head_len();
	int i, j, k;
	Allele a, b;
	double total_likelihood;
	HaploPair *best_hp;
	vector<HaploPair*>::iterator i_hp;
	for_each(m_haplopairs.begin(), m_haplopairs.end(), DeleteAll_Clear());
	m_haplopairs.resize(genotype_len()+1);
	m_best_pair.resize(pn);
	for (i=0; i<pn; ++i) {
		m_best_pair[i].clear();
	}
	initHeadList(genotype);
	for (i=head_len; i<genotype_len(); ++i) {
		if (genotype.isMissing(i)) {
			for (j=0; j<m_haplodata->allele_num(i); ++j) {
				if (m_haplodata->allele_frequency(i, j) > 0) {
					for (k=j; k<m_haplodata->allele_num(i); ++k) {
						if (m_haplodata->allele_frequency(i, k) > 0) {
							a = m_haplodata->allele_symbol(i, j);
							b = m_haplodata->allele_symbol(i, k);
							extendAll(i, a, b);
						}
					}
				}
			}
		}
		else if (genotype(0)[i].isMissing()) {
			for (j=0; j<m_haplodata->allele_num(i); ++j) {
				if (m_haplodata->allele_frequency(i, j) > 0) {
					a = m_haplodata->allele_symbol(i, j);
					extendAll(i, a, genotype(1)[i]);
				}
			}
		}
		else if (genotype(1)[i].isMissing()) {
			for (j=0; j<m_haplodata->allele_num(i); ++j) {
				if (m_haplodata->allele_frequency(i, j) > 0) {
					a = m_haplodata->allele_symbol(i, j);
					extendAll(i, a, genotype(0)[i]);
				}
			}
		}
		else {
			extendAll(i, genotype(0)[i], genotype(1)[i]);
		}
		if (m_haplopairs[i+1].size() <= 0) {
			break;
		}
	}
	if (m_haplopairs[genotype_len()].size() > 0) {
		total_likelihood = 0;
		res_list.clear();
		for (i_hp = m_haplopairs[genotype_len()].begin(); i_hp != m_haplopairs[genotype_len()].end(); ++i_hp) {
			total_likelihood += (*i_hp)->forward_likelihood();
			res_list.push_back(*i_hp);
		}
		sort(res_list.begin(), res_list.end(), HaploPair::greater_likelihood());
		best_hp = res_list.front();
		resolution = best_hp->getGenotype();
		resolution.setLikelihood(total_likelihood);
		resolution.setWeight(best_hp->best_likelihood() / total_likelihood);
	}
	else {
		res_list.clear();
		resolution = genotype;
		resolution.setLikelihood(0);
		resolution.setWeight(0);
	}
}

double HaploBuilder::getLikelihood(const Haplotype &haplotype)
{
	int head_len = m_patterns.head_len();
	const BackwardPatternTree *pattern_tree = m_patterns.pattern_tree();
	int i;
	HaploPattern *hp;
	double likelihood = 1.0;
	for (i=head_len; i<=genotype_len(); ++i) {
		hp = pattern_tree->findLongestMatchPattern(i, &haplotype);
		if (hp) {
			likelihood *= hp->transition_prob();
		}
		else {
			likelihood = 0;
			break;
		}
	}
	return likelihood;
}

double HaploBuilder::getLikelihood(const Genotype &genotype)
{
	return getLikelihood(genotype(0)) * getLikelihood(genotype(1));
}

void HaploBuilder::initHeadList(const Genotype &genotype)
{
	int head_len = m_patterns.head_len();
	const vector<HaploPattern*> &head_list = m_patterns.head_list();
	const BackwardPatternTree *pattern_tree = m_patterns.pattern_tree();
	vector<AlleleSequence*> new_list, last_list;
	vector<AlleleSequence*>::iterator i_as;
	vector<HaploPattern*>::const_iterator head;
	for (head = head_list.begin(); head != head_list.end(); ++head) {
		if ((*head)->isMatch(genotype)) {
			last_list.push_back(new AlleleSequence);
			for (int j=0; j<head_len; ++j) {
				i_as = last_list.begin();
				if (genotype.isMissing(j) || (genotype.hasMissing(j) && genotype.hasAllele(j, (**head)[j]))) {
					while (i_as != last_list.end()) {
						AlleleSequence *as = *i_as;
						for (int k=0; k<m_haplodata->allele_num(j); ++k) {
							if (m_haplodata->allele_frequency(j, k) > 0) {
								AlleleSequence *new_as = new AlleleSequence;
								new_as->assign(*as, m_haplodata->allele_symbol(j, k));
								new_list.push_back(new_as);
							}
						}
						++i_as;
					}
				}
				else if (genotype.isHeterozygous(j)) {
					while (i_as != last_list.end()) {
						AlleleSequence *as = *i_as;
						AlleleSequence *new_as = new AlleleSequence;
						if ((**head)[j] == genotype(0)[j]) {
							new_as->assign(*as, genotype(1)[j]);
						}
						else {
							new_as->assign(*as, genotype(0)[j]);
						}
						new_list.push_back(new_as);
						++i_as;
					}
				}
				else {
					while (i_as != last_list.end()) {
						AlleleSequence *as = *i_as;
						AlleleSequence *new_as = new AlleleSequence;
						new_as->assign(*as, genotype(0)[j]);
						new_list.push_back(new_as);
						++i_as;
					}
				}
				DeleteAll_Clear()(last_list);
				last_list.swap(new_list);
			}
			i_as = last_list.begin();
			while (i_as != last_list.end()) {
				HaploPattern *hp = pattern_tree->findLongestMatchPattern(head_len, *i_as);
				if (hp && hp->start() == 0) {
					HaploPair *new_hp = new HaploPair(*head, hp);
					m_haplopairs[head_len].push_back(new_hp);
				 	m_best_pair[new_hp->id_a()].insert(make_pair(new_hp->id_b(), m_haplopairs[head_len].size()));
				}
				else {
					Logger::error("Can not find matching pattern!");
					exit(1);
				}
				++i_as;
			}
			DeleteAll_Clear()(last_list);
		}
	}
}

void HaploBuilder::extendAll(int i, Allele a1, Allele a2)
{
	vector<HaploPair*>::iterator i_hp;
	for (i_hp = m_haplopairs[i].begin(); i_hp != m_haplopairs[i].end(); ++i_hp) {
		extend(*i_hp, a1, a2);
		if (a1 != a2) extend(*i_hp, a2, a1);
	}
}

void HaploBuilder::extend(HaploPair *hp, Allele a1, Allele a2)
{
	const HaploPattern *hpa, *hpb;
	hpa = hp->successor_a(a1);
	hpb = hp->successor_b(a2);
	if (hpa && hpb) {
		map<int, int>::iterator i = m_best_pair[hpa->id()].lower_bound(hpb->id());
		if (i == m_best_pair[hpa->id()].end() || (*i).first != hpb->id()) {
			m_haplopairs[hp->end()+1].push_back(new HaploPair(hp, hpa, hpb));
		 	m_best_pair[hpa->id()].insert(i, make_pair(hpb->id(), m_haplopairs[hp->end()+1].size()));
		}
		else {
			m_haplopairs[hp->end()+1][(*i).second-1]->add(hp, hpa, hpb);
		}
	}
}

void HaploBuilder::calcBackwardLikelihood()
{
	int head_len = m_patterns.head_len();
	vector<HaploPair*>::iterator i_hp;
	for (int i=genotype_len()-1; i>=head_len; --i) {
		for (i_hp=m_haplopairs[i].begin(); i_hp!=m_haplopairs[i].end(); ++i_hp) {
			(*i_hp)->calcBackwardLikelihood();
		}
	}
}

void HaploBuilder::estimateFrequency(vector<HaploPattern*> &patterns)
{
	int i, n;
	int start, geno;
	Genotype res;
	vector<HaploPair*> res_list;
	ForwardPatternTree tree(*m_haplodata);
	map<HaploPair*, double> match_list[2];

	n = patterns.size();
	for (i=0; i<n; ++i) {
		HaploPattern *hp = patterns[i];
		tree.addPattern(hp);
		hp->setFrequency(0);
		hp->setPrefixFreq(0);
	}

	for (geno=0; geno<genotype_num(); ++geno) {
		resolve((*m_haplodata)[geno], res, res_list);
		calcBackwardLikelihood();
		m_current_genotype_frequency = (*m_haplodata)[geno].likelihood();

		for (start=0; start<genotype_len(); ++start) {
			match_list[0].clear();
			match_list[1].clear();

			int end = max(start, m_patterns.head_len());
			n = m_haplopairs[end].size();
			for (i=0; i<n; ++i) {
				HaploPair *hp = m_haplopairs[end][i];
				match_list[0][hp] = hp->forward_likelihood();
			}

			PatternNode *node = tree.root(start);
			n = node->size();
			for (i=0; i<node->size(); ++i) {
				if (node->getChild(i)) {
					estimateFrequency(node->getChild(i), start, m_haplodata->allele_symbol(start, i), 1.0, match_list);
				}
			}
		}
	}

	n = patterns.size();
	for (i=0; i<n; ++i) {
		HaploPattern *hp = patterns[i];
		if (hp->prefix_freq() > 0) {
			hp->setTransitionProb(hp->frequency() / hp->prefix_freq());
		}
		else {
			hp->setTransitionProb(hp->frequency() / genotype_num());
		}
		hp->setFrequency(hp->frequency() / genotype_num());
		hp->setPrefixFreq(hp->prefix_freq() / genotype_num());
	}
}

void HaploBuilder::estimateFrequency(PatternNode *node, int locus, const Allele &a, double last_freq, const map<HaploPair*, double> last_match[2])
{
	map<HaploPair*, double> match_list[2];
	map<HaploPair*, double>::const_iterator i_mp;
	vector<HaploPair*>::const_iterator i_link;

	if (locus < m_patterns.head_len()) {
		for (i_mp=last_match[0].begin(); i_mp!=last_match[0].end(); ++i_mp) {
			HaploPair *hp = (*i_mp).first;
			double weight = (*i_mp).second;
			if (hp->pattern_a()[locus] == a) {
				if (hp->pattern_b()[locus] == a) {
					match_list[0][hp] += weight;
				}
				else {
					match_list[1][hp] += weight;
				}
			}
		}
		for (i_mp=last_match[1].begin(); i_mp!=last_match[1].end(); ++i_mp) {
			HaploPair *hp = (*i_mp).first;
			double weight = (*i_mp).second;
			if (hp->pattern_a()[locus] == a) {
				match_list[1][hp] += weight;
			}
		}
	}
	else {
		for (i_mp=last_match[0].begin(); i_mp!=last_match[0].end(); ++i_mp) {
			HaploPair *last_hp = (*i_mp).first;
			double weight = (*i_mp).second;
			for (i_link=last_hp->forward_links().begin(); i_link!=last_hp->forward_links().end(); ++i_link) {
				HaploPair *hp = *i_link;
				if (hp->allele_a() == a) {
					if (hp->allele_b() == a) {
						match_list[0][hp] += weight * hp->transition_prob();
					}
					else {
						match_list[1][hp] += weight * hp->transition_prob();
					}
				}
			}
		}
		for (i_mp=last_match[1].begin(); i_mp!=last_match[1].end(); ++i_mp) {
			HaploPair *last_hp = (*i_mp).first;
			double weight = (*i_mp).second;
			for (i_link=last_hp->forward_links().begin(); i_link!=last_hp->forward_links().end(); ++i_link) {
				HaploPair *hp = *i_link;
				if (hp->allele_a() == a) {
					match_list[1][hp] += weight * hp->transition_prob();
				}
			}
		}
	}

	double freq = 0;
	for (i_mp=match_list[0].begin(); i_mp!=match_list[0].end(); ++i_mp) {
		freq += (*i_mp).second * (*i_mp).first->backward_likelihood();
	}
	for (i_mp=match_list[1].begin(); i_mp!=match_list[1].end(); ++i_mp) {
		freq += (*i_mp).second * (*i_mp).first->backward_likelihood();
	}
	freq /= m_current_genotype_frequency;

	if (node->data()) {
		HaploPattern *hp = node->data();
		hp->setFrequency(hp->frequency() + freq);
		hp->setPrefixFreq(hp->prefix_freq() + last_freq);
	}

	for (int i=0; i<node->size(); ++i) {
		if (node->getChild(i)) {
			estimateFrequency(node->getChild(i), locus+1, m_haplodata->allele_symbol(locus+1, i), freq, match_list);
		}
	}
}
