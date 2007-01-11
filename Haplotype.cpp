
#include "Haplotype.h"

#include "MemLeak.h"


////////////////////////////////
//
// class Haplotype

Haplotype::Haplotype()
{
	m_id[0] = 0;
	m_weight = 1.0;
}

Haplotype::Haplotype(int len)
	: AlleleSequence(len)
{
	m_id[0] = 0;
	m_weight = 1.0;
}

Haplotype::Haplotype(const AlleleSequence &as)
	: AlleleSequence(as)
{
	m_id[0] = 0;
	m_weight = 1.0;
}

Haplotype::Haplotype(const Haplotype &h)
	: AlleleSequence(h)
{
	strcpy(m_id, h.m_id);
	m_weight = h.m_weight;
}

Haplotype &Haplotype::assign(const Haplotype &h)
{
	strcpy(m_id, h.m_id);
	m_weight = h.m_weight;
	AlleleSequence::assign(h);
	return *this;
}

Haplotype &Haplotype::concatenate(const Haplotype &h)
{
	AlleleSequence::concatenate(h);
	return *this;
}

Haplotype &Haplotype::concatenate(const Allele &a)
{
	AlleleSequence::concatenate(a);
	return *this;
}

Haplotype &Haplotype::concatenate(const Haplotype &h1, const Haplotype &h2)
{
	strcpy(m_id, h1.m_id);
	AlleleSequence::concatenate(h1, h2);
	return *this;
}

Haplotype &Haplotype::concatenate(const Haplotype &h, const Allele &a)
{
	strcpy(m_id, h.m_id);
	AlleleSequence::concatenate(h, a);
	return *this;
}

Haplotype &Haplotype::concatenate(const Allele &a, const Haplotype &h)
{
	strcpy(m_id, h.m_id);
	AlleleSequence::concatenate(a, h);
	return *this;
}
