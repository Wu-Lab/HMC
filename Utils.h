
#ifndef __UTILS_H
#define __UTILS_H


#include <ctime>
#include <boost/shared_ptr.hpp>

#include "Tree.h"
#include "Matrix.h"


using namespace std;


namespace std {
	namespace tr1 = ::boost;
}


class Timer {
	static int m_instance_num;
	int m_id;
	clock_t m_start, m_finish;
	bool m_running;
	double m_duration;
	char m_description[80];

public:
	Timer() { m_id = m_instance_num++; }
	void begin(const char *description);
	void end();
	void pause();
	void resume();
	double time() { return m_duration; }
};


class Logger {
	static bool m_logging;
	static int m_log_level;
	static const int m_log_level_error;
	static const int m_log_level_warning;
	static const int m_log_level_info;
	static const int m_log_level_verbose;
	static const int m_log_level_debug;

	static Timer m_timer[10];

public:
	static int log_level() { return m_log_level; }
	static int log_level_error() { return m_log_level_error; }
	static int log_level_warning() { return m_log_level_warning; }
	static int log_level_info() { return m_log_level_info; }
	static int log_level_verbose() { return m_log_level_verbose; }
	static int log_level_debug() { return m_log_level_debug; }

	static Timer &timer(int i) { return m_timer[i]; }

	static void error(const char *format, ...);
	static void warning(const char *format, ...);
	static void info(const char *format, ...);
	static void verbose(const char *format, ...);
	static void status(const char *format, ...);
	static void debug(const char *format, ...);

	static void print(int level, FILE *fp, const char *prompt, const char *format, ...);
	static void println(int level, FILE *fp, const char *prompt, const char *format, ...);

	static void setLogLevel(int level);
	static void enableLogging();
	static void disableLogging();

	static bool isDebug() { return (m_log_level >= m_log_level_debug); }

	static void beginTimer(int i, const char *description);
	static void endTimer(int i);
	static void pauseTimer(int i);
	static void resumeTimer(int i);

private:
	static void _print(int level, FILE *fp, const char *prompt, const char *format, va_list argptr);
	static void _println(int level, FILE *fp, const char *prompt, const char *format, va_list argptr);
};


// something for using with STL

struct DeletePtr {
	template<typename T>
	void operator()(const T* ptr) const { delete ptr; }
};

template<typename T>
void DeleteAllObjects(const T& container) {
	for_each(container.begin(), container.end(), DeletePtr());
};


#endif // __UTILS_H
