
#include "HaploComp.h"

#include "MemLeak.h"


HaploComp::HaploComp()
{
	m_haplo_real = NULL;
	m_haplo_infer = NULL;
	m_haplo_input = NULL;
	m_genotype_num = 0;
	m_genotype_len = 0;
	m_switch_error = 0;
	m_incorrect_genotype_percentage = 0;
	m_incorrect_haplotype_percentage = 0;
	m_missing_error = 0;
	m_k2_distance = 0;
	m_switch_error_numerator = 0;
	m_switch_error_denominator = 0;
	m_incorrect_genotype_numerator = 0;
	m_incorrect_genotype_denominator = 0;
	m_incorrect_haplotype_numerator = 0;
	m_incorrect_haplotype_denominator = 0;
	m_missing_error_numerator = 0;
	m_missing_error_denominator = 0;
}

HaploComp::HaploComp(const HaploData *real, const HaploData *infer, const HaploData *input)
{	
	int i, sd, ig;
	m_haplo_real = real;
	m_haplo_infer = infer;
	if (input == NULL) m_haplo_input = real;
	if (m_haplo_real->genotype_num() != m_haplo_infer->genotype_num() ||
		m_haplo_real->genotype_len() != m_haplo_infer->genotype_len()) {
		Logger::error("Attempt to compare inconsistent haplotype data!");
		exit(1);
	}
	m_genotype_num = m_haplo_real->unphased_num();
	m_genotype_len = m_haplo_real->genotype_len();
	m_switch_error = 0;
	m_incorrect_genotype_percentage = 0;
	m_incorrect_haplotype_percentage = 0;
	m_missing_error = 0;
	m_k2_distance = 0;
	m_switch_error_numerator = 0;
	m_switch_error_denominator = 0;
	m_incorrect_genotype_numerator = 0;
	m_incorrect_genotype_denominator = 0;
	m_incorrect_haplotype_numerator = 0;
	m_incorrect_haplotype_denominator = 0;
	m_missing_error_numerator = 0;
	m_missing_error_denominator = 0;
	for (i=0; i<m_genotype_num; i++) {
		const Genotype &geno_real = (*m_haplo_real)[i];
		const Genotype &geno_infer = (*m_haplo_infer)[i];
		const Genotype &geno_input = (*m_haplo_input)[i];
		// Switch Error
		sd = geno_real.getSwitchDistanceIgnoreMissing(geno_infer);
		m_switch_error_numerator += sd;
		m_switch_error_denominator += geno_real.heterozygous_num() - 1;
		// IGP
		ig = geno_real.getDiffNumIgnoreMissing(geno_infer);
		m_incorrect_genotype_numerator += ig;
		m_incorrect_genotype_denominator += m_genotype_len - geno_real.missing_num();
		// IHP
		if (sd > 0) m_incorrect_haplotype_numerator++;
		if (geno_real.heterozygous_num() > 1) m_incorrect_haplotype_denominator++;
		// Missing Error
		if (m_haplo_input != m_haplo_real) {
			getMissingError(geno_real, geno_infer, geno_input);
		}
	}
	calculate();
}

HaploComp &HaploComp::operator +=(const HaploComp &hc)
{
	m_switch_error_numerator += hc.m_switch_error_numerator;
	m_switch_error_denominator += hc.m_switch_error_denominator;
	m_incorrect_genotype_numerator += hc.m_incorrect_genotype_numerator;
	m_incorrect_genotype_denominator += hc.m_incorrect_genotype_denominator;
	m_incorrect_haplotype_numerator += hc.m_incorrect_haplotype_numerator;
	m_incorrect_haplotype_denominator += hc.m_incorrect_haplotype_denominator;
	m_missing_error_numerator += hc.m_missing_error_numerator;
	m_missing_error_denominator += hc.m_missing_error_denominator;
	calculate();
	return *this;
}

void HaploComp::getMissingError(const Genotype &real, const Genotype &infer, const Genotype &input)
{
	int i, diff1, diff2;
	diff1 = diff2 = 0;
	for (i=0; i<m_genotype_len; i++) {
		if (!input.hasMissing(i)) {
			if (!real.isMatch(infer, i, true)) {
				diff1++;
			}
			if (!real.isMatch(infer, i, false)) {
				diff2++;
			}
		}
	}
	if (diff1 < diff2) {
		for (i=0; i<m_genotype_len; i++) {
			if (input.isMissing(i)) {
				if (!real(0).isMissing(i)) {
					m_missing_error_denominator++;
					if (!real(0).isMatch(infer(0)[i], i)) m_missing_error_numerator++;			
				};
				if (!real(1).isMissing(i)) {
					m_missing_error_denominator++;
					if (!real(1).isMatch(infer(1)[i], i)) m_missing_error_numerator++;			
				};
			}
			else if (input.hasMissing(i) && !real.hasMissing(i)) {
				m_missing_error_denominator++;
				if (!real.isMatch(infer, i, true)) m_missing_error_numerator++;
			}
		}
	}
	else {
		for (i=0; i<m_genotype_len; i++) {
			if (input.isMissing(i)) {
				if (!real(0).isMissing(i)) {
					m_missing_error_denominator++;
					if (!real(0).isMatch(infer(1)[i], i)) m_missing_error_numerator++;			
				};
				if (!real(1).isMissing(i)) {
					m_missing_error_denominator++;
					if (!real(1).isMatch(infer(0)[i], i)) m_missing_error_numerator++;			
				};
			}
			else if (input.hasMissing(i) && !real.hasMissing(i)) {
				m_missing_error_denominator++;
				if (!real.isMatch(infer, i, false)) m_missing_error_numerator++;
			}
		}
	}
}

void HaploComp::calculate()
{
	m_switch_error = (double) m_switch_error_numerator / m_switch_error_denominator;
	m_incorrect_genotype_percentage = (double) m_incorrect_genotype_numerator / m_incorrect_genotype_denominator;
	m_incorrect_haplotype_percentage = (double) m_incorrect_haplotype_numerator / m_incorrect_haplotype_denominator;
	if (m_haplo_input != m_haplo_real && m_missing_error_denominator > 0) {
		m_missing_error = (double) m_missing_error_numerator / m_missing_error_denominator;
	}
	else {
		m_missing_error = 0;
	}
}
