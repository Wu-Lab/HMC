
#ifndef __CONSTANT_H
#define __CONSTANT_H


#include "Utils.h"


#define STR_LEN_FILE_NAME 256
#define STR_LEN_FILE_LINE 10000
#define STR_LEN_ALLELE_NAME 64
#define STR_LEN_GENOTYPE_ID 256


class Constant {
protected:
	static int m_average_marker_distance;

public:
	static int average_marker_distance() { return m_average_marker_distance; }

	static void setAverageMarkerDistance(int i) { m_average_marker_distance = i; }
};


#endif // __CONSTANT_H
