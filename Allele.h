
#ifndef __ALLELE_H
#define __ALLELE_H


#include "Utils.h"


typedef int Allele;


class AlleleSequence {
protected:
	Allele *m_alleles;
	int m_length;

public:
	AlleleSequence();
	AlleleSequence(const AlleleSequence &as);
	explicit AlleleSequence(int len);
	~AlleleSequence();

	Allele &operator [](int i) { return m_alleles[i]; }
	const Allele &operator [](int i) const { return m_alleles[i]; }
	int length() const { return m_length; }

	bool isMissing(int i) const { return (m_alleles[i] < 0); }
	bool isMatch(const Allele &a, int locus) const;
	bool isMatch(const AlleleSequence &as, int start1, int start2, int len) const;
	bool isMatch(const AlleleSequence &as) const;
	int getDiffNum(const AlleleSequence &as, int start1, int start2, int len) const;
	int getDiffNum(const AlleleSequence &as) const;

	int setLength(int i);

	char *read(const char *types, char *buffer, int len = 0);
	char *write(const char *types, char *buffer) const;

	AlleleSequence &operator =(const AlleleSequence &as);
	AlleleSequence &operator +=(const AlleleSequence &as);
	AlleleSequence &operator +=(const Allele &a);

protected:
	static char *readAllele(const char type, char *buffer, int &allele);
	static char *writeAllele(const char type, char *buffer, int &allele);

	AlleleSequence &assign(const AlleleSequence &as);
	AlleleSequence &concatenate(const AlleleSequence &as);
	AlleleSequence &concatenate(const Allele &a);

public:
	friend AlleleSequence operator +(const AlleleSequence &as1, const AlleleSequence &as2);
	friend AlleleSequence operator +(const AlleleSequence &as, const Allele &a);
	friend AlleleSequence operator +(const Allele &a, const AlleleSequence &as);

	AlleleSequence &concatenate(const AlleleSequence &as1, const AlleleSequence &as2);
	AlleleSequence &concatenate(const AlleleSequence &as, const Allele &a);
	AlleleSequence &concatenate(const Allele &a, const AlleleSequence &as);
};

inline bool AlleleSequence::isMatch(const Allele &a, int locus) const
{
	return (isMissing(locus) || a < 0 || m_alleles[locus] == a);
}

inline bool AlleleSequence::isMatch(const AlleleSequence &as) const
{
	return (m_length == as.m_length && isMatch(as, 0, 0, m_length));
}

inline int AlleleSequence::getDiffNum(const AlleleSequence &as) const
{
	return getDiffNum(as, 0, 0, (m_length < as.m_length ? m_length : as.m_length));
}

inline AlleleSequence &AlleleSequence::operator =(const AlleleSequence &as)
{
	return assign(as);
}

inline AlleleSequence &AlleleSequence::operator +=(const AlleleSequence &as)
{
	return concatenate(as);
}

inline AlleleSequence &AlleleSequence::operator +=(const Allele &a)
{
	return concatenate(a);
}

inline AlleleSequence operator +(const AlleleSequence &as1, const AlleleSequence &as2)
{
	return AlleleSequence().concatenate(as1, as2);
}

inline AlleleSequence operator +(const AlleleSequence &as, const Allele &a)
{
	return AlleleSequence().concatenate(as, a);
}

inline AlleleSequence operator +(const Allele &a, const AlleleSequence &as)
{
	return AlleleSequence().concatenate(a, as);
}


#endif // __ALLELE_H
