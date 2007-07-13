
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
	double weight() const { return m_weight; }

	void setID(const string &id) { m_id = id; string_replace(m_id, " ", "_"); }
	void setWeight(double weight) { m_weight = weight; }
};

inline Haplotype::Haplotype()
: m_weight(1.0)
{
}

inline Haplotype::Haplotype(int len)
: AlleleSequence(len),
  m_weight(1.0)
{
}

inline Haplotype::Haplotype(const AlleleSequence &as)
: AlleleSequence(as),
  m_weight(1.0)
{
}


#endif // __HAPLOTYPE_H
