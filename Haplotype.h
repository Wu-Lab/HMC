
#ifndef __HAPLOTYPE_H
#define __HAPLOTYPE_H


#include "Utils.h"
#include "Constant.h"
#include "Allele.h"


class Haplotype : public AlleleSequence {
protected:
	char m_id[STR_LEN_GENOTYPE_ID];
	double m_weight;

public:
	Haplotype();
	Haplotype(const Haplotype &h);
	explicit Haplotype(int len);
	explicit Haplotype(const AlleleSequence &as);

	char *id() const { return (char *) m_id; }
	double weight() const { return m_weight; }

	void setID(const char *id) { strcpy(m_id, id); }
	void setWeight(double weight) { m_weight = weight; }

	Haplotype &assign(const Haplotype &h);
	Haplotype &assign(const Haplotype &h1, const Haplotype &h2);
	Haplotype &assign(const Haplotype &h, const Allele &a);
	Haplotype &assign(const Allele &a, const Haplotype &h);

	Haplotype &operator +=(const Haplotype &h);
	Haplotype &operator +=(const Allele &a);

	friend class HaploPattern;
	friend class HaploData;
};

inline Haplotype &Haplotype::assign(const Haplotype &h)
{
	strcpy(m_id, h.m_id);
	m_weight = h.m_weight;
	AlleleSequence::assign(h);
	return *this;
}

inline Haplotype &Haplotype::assign(const Haplotype &h1, const Haplotype &h2)
{
	strcpy(m_id, h1.m_id);
	AlleleSequence::assign(h1, h2);
	return *this;
}

inline Haplotype &Haplotype::assign(const Haplotype &h, const Allele &a)
{
	strcpy(m_id, h.m_id);
	AlleleSequence::assign(h, a);
	return *this;
}

inline Haplotype &Haplotype::assign(const Allele &a, const Haplotype &h)
{
	strcpy(m_id, h.m_id);
	AlleleSequence::assign(a, h);
	return *this;
}

inline Haplotype &Haplotype::operator +=(const Haplotype &h)
{
	AlleleSequence::operator +=(h);
	return *this;
}

inline Haplotype &Haplotype::operator +=(const Allele &a)
{
	AlleleSequence::operator +=(a);
	return *this;
}


#endif // __HAPLOTYPE_H
