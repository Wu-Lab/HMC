
#ifndef __HAPLOTYPE_H
#define __HAPLOTYPE_H


#include <string>

#include "Utils.h"
#include "Constant.h"
#include "Allele.h"


class Haplotype : public AlleleSequence {
protected:
	string m_id;
	double m_weight;

public:
	Haplotype();
	explicit Haplotype(int len);
	explicit Haplotype(const AlleleSequence &as);

	const string &id() const { return m_id; }
	const char *id_str() const { return m_id.c_str(); }
	double weight() const { return m_weight; }

	void setID(const string &id) { m_id = id; }
	void setWeight(double weight) { m_weight = weight; }

	Haplotype &assign(const Haplotype &h1, const Haplotype &h2);
	Haplotype &assign(const Haplotype &h, const Allele &a);
	Haplotype &assign(const Allele &a, const Haplotype &h);

	friend class HaploPattern;
	friend class HaploData;
};

Haplotype::Haplotype()
: m_weight(1.0)
{
}

Haplotype::Haplotype(int len)
: AlleleSequence(len),
  m_weight(1.0)
{
}

inline Haplotype::Haplotype(const AlleleSequence &as)
: AlleleSequence(as),
  m_weight(1.0)
{
}

inline Haplotype &Haplotype::assign(const Haplotype &h1, const Haplotype &h2)
{
	m_id = h1.m_id;
	AlleleSequence::assign(h1, h2);
	return *this;
}

inline Haplotype &Haplotype::assign(const Haplotype &h, const Allele &a)
{
	m_id = h.m_id;
	AlleleSequence::assign(h, a);
	return *this;
}

inline Haplotype &Haplotype::assign(const Allele &a, const Haplotype &h)
{
	m_id = h.m_id;
	AlleleSequence::assign(a, h);
	return *this;
}


#endif // __HAPLOTYPE_H
