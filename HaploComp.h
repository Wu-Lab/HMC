
#ifndef __HAPLOCOMP_H
#define __HAPLOCOMP_H


#include "Utils.h"
#include "GenoData.h"


class HaploComp {
	const GenoData *m_genos_real, *m_genos_infer, *m_genos_input;
	int m_genotype_num;
	int m_genotype_len;

	double m_switch_error;
	double m_incorrect_genotype_percentage;
	double m_incorrect_haplotype_percentage;
	double m_missing_error;
	double m_k2_distance;

	int m_switch_error_numerator;
	int m_switch_error_denominator;
	int m_incorrect_genotype_numerator;
	int m_incorrect_genotype_denominator;
	int m_incorrect_haplotype_numerator;
	int m_incorrect_haplotype_denominator;
	int m_missing_error_numerator;
	int m_missing_error_denominator;

public:
	HaploComp();
	explicit HaploComp(const GenoData *real, const GenoData *infer, const GenoData *input = NULL);

	double switch_error() const { return m_switch_error; }
	double incorrect_genotype_percentage() const { return m_incorrect_genotype_percentage; }
	double incorrect_haplotype_percentage() const { return m_incorrect_haplotype_percentage; }
	double missing_error() const { return m_missing_error; }
	double k2_distance() const { return m_k2_distance; }

	HaploComp &operator +=(const HaploComp &hc);

protected:
	void getMissingError(const Genotype &real, const Genotype &infer, const Genotype &input);

	void calculate();
};


#endif // __HAPLOCOMP_H
