
#ifndef __CONSTANT_H
#define __CONSTANT_H


#include "Utils.h"


#define STR_LEN_FILE_NAME 256
#define STR_LEN_FILE_LINE 10000
#define STR_LEN_ALLELE_NAME 64
#define STR_LEN_GENOTYPE_ID 256


class Constant {
protected:
	static int m_max_allele_num;
	static int m_ram_limit;
	static int m_average_marker_distance;

public:
	static int max_allele_num() { return m_max_allele_num; }
	static int ram_limit() { return m_ram_limit; }
	static int average_marker_distance() { return m_average_marker_distance; }

	static void setRAMLimit(int i) { m_ram_limit = i; }
	static void setMaxAlleleNum(int i) { m_max_allele_num = i; }
	static void setAverageMarkerDistance(int i) { m_average_marker_distance = i; }
};


#endif // __CONSTANT_H
