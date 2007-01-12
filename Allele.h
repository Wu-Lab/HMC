
#ifndef __ALLELE_H
#define __ALLELE_H


#include <vector>

#include "Utils.h"


class Allele {
	int m_value;

public:
	Allele() { m_value = -1; }
	explicit Allele(int a) { m_value = a; }
	Allele &operator =(int a) { m_value = a; return *this; }

	bool isMissing() const { return m_value < 0; }
	bool isMatch(Allele a) const { return (isMissing() || a.isMissing() || m_value == a.m_value); }

	friend bool operator ==(const Allele &lhs, const Allele &rhs);
	friend bool operator !=(const Allele &lhs, const Allele &rhs);

	int asInt() const { return m_value; }
	char asChar() const { return m_value; }
};

inline bool operator ==(const Allele &lhs, const Allele &rhs)
{
	return ((lhs.m_value == rhs.m_value) || (lhs.isMissing() && rhs.isMissing()));
}

inline bool operator !=(const Allele &lhs, const Allele &rhs)
{
	return !(lhs == rhs);
}


class AlleleSequence {
protected:
	vector<Allele> m_alleles;

public:
	AlleleSequence() { }
	explicit AlleleSequence(int len);
	explicit AlleleSequence(const Allele &a);
	explicit AlleleSequence(const AlleleSequence &as1, const AlleleSequence &as2);

	Allele &operator [](int i) { return m_alleles[i]; }
	const Allele &operator [](int i) const { return m_alleles[i]; }
	int length() const { return m_alleles.size(); }

	bool isMatch(const AlleleSequence &as, int start1, int start2, int len) const;
	bool isMatch(const AlleleSequence &as) const;
	int getDiffNum(const AlleleSequence &as, int start1, int start2, int len) const;
	int getDiffNum(const AlleleSequence &as) const;

	int setLength(int len);

	char *read(const char *types, char *buffer, int len = 0);
	char *write(const char *types, char *buffer) const;

	AlleleSequence &assign(const AlleleSequence &as, const Allele &a);
	AlleleSequence &assign(const Allele &a, const AlleleSequence &as);
	AlleleSequence &assign(const AlleleSequence &as1, const AlleleSequence &as2);

	AlleleSequence &operator +=(const AlleleSequence &as);
	AlleleSequence &operator +=(const Allele &a);

protected:
	static char *readAllele(char type, char *buffer, Allele &allele);
	static char *writeAllele(char type, char *buffer, const Allele &allele);
};

inline AlleleSequence::AlleleSequence(int len)
{
	if (len > 0) {
		m_alleles.resize(len);
	}
	else {
		m_alleles.clear();
	}
}

inline AlleleSequence::AlleleSequence(const Allele &a)
{
	m_alleles.push_back(a);
}

inline AlleleSequence::AlleleSequence(const AlleleSequence &as1, const AlleleSequence &as2)
: m_alleles(as1.m_alleles)
{
	(*this) += as2;
}

inline bool AlleleSequence::isMatch(const AlleleSequence &as) const
{
	return (length() == as.length() && isMatch(as, 0, 0, length()));
}

inline int AlleleSequence::getDiffNum(const AlleleSequence &as) const
{
	return getDiffNum(as, 0, 0, (length() < as.length() ? length() : as.length()));
}

inline AlleleSequence &AlleleSequence::assign(const AlleleSequence &as, const Allele &a)
{
	m_alleles = as.m_alleles;
	(*this) += a;
	return *this;
}

inline AlleleSequence &AlleleSequence::assign(const Allele &a, const AlleleSequence &as)
{
	m_alleles.clear();
	m_alleles.push_back(a);
	(*this) += as;
	return *this;
}

inline AlleleSequence &AlleleSequence::assign(const AlleleSequence &as1, const AlleleSequence &as2)
{
	m_alleles = as1.m_alleles;
	(*this) += as2;
	return *this;
}

inline AlleleSequence &AlleleSequence::operator +=(const AlleleSequence &as)
{
	m_alleles.insert(m_alleles.end(), as.m_alleles.begin(), as.m_alleles.end());
	return *this;
}

inline AlleleSequence &AlleleSequence::operator +=(const Allele &a)
{
	m_alleles.push_back(a);
	return *this;
}


#endif // __ALLELE_H
