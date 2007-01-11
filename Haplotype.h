
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

	Haplotype &operator =(const Haplotype &h);
	Haplotype &operator *=(const Haplotype &h);
	Haplotype &operator +=(const Haplotype &h);
	Haplotype &operator +=(const Allele &a);

protected:
	Haplotype &assign(const Haplotype &h);
	Haplotype &concatenate(const Haplotype &h);
	Haplotype &concatenate(const Allele &a);

public:
	friend Haplotype operator +(const Haplotype &h1, const Haplotype &h2);
	friend Haplotype operator +(const Haplotype &h, const Allele &a);
	friend Haplotype operator +(const Allele &a, const Haplotype &h);

	Haplotype &concatenate(const Haplotype &h1, const Haplotype &h2);
	Haplotype &concatenate(const Haplotype &h, const Allele &a);
	Haplotype &concatenate(const Allele &a, const Haplotype &h);

	friend class HaploPattern;
	friend class HaploData;
};

inline Haplotype &Haplotype::operator =(const Haplotype &h)
{
	return assign(h);
}

inline Haplotype &Haplotype::operator +=(const Haplotype &h)
{
	return concatenate(h);
}

inline Haplotype &Haplotype::operator +=(const Allele &a)
{
	return concatenate(a);
}

inline Haplotype operator +(const Haplotype &h1, const Haplotype &h2)
{
	return Haplotype().concatenate(h1, h2);
}

inline Haplotype operator +(const Haplotype &h, const Allele &a)
{
	return Haplotype().concatenate(h, a);
}

inline Haplotype operator +(const Allele &a, const Haplotype &h)
{
	return Haplotype().concatenate(a, h);
}


#endif // __HAPLOTYPE_H
