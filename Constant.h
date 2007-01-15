
#ifndef __CONSTANT_H
#define __CONSTANT_H


class Constant {
protected:
	static int m_average_marker_distance;

public:
	static int average_marker_distance() { return m_average_marker_distance; }

	static void setAverageMarkerDistance(int i) { m_average_marker_distance = i; }
};


#endif // __CONSTANT_H
