
#ifndef __UTILS_H
#define __UTILS_H


#include <new>
#include <ctime>
#include <vector>
#include <iosfwd>
#include <iterator>
#include <boost/shared_ptr.hpp>


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


class StandardNewDelete {
public:
	// normal new/delete
	static void *operator new(std::size_t size) throw(std::bad_alloc)
	{ return ::operator new(size); }
	static void operator delete(void *pMemory) throw()
	{ ::operator delete(pMemory); }

	// placement new/delete
	static void *operator new(std::size_t size, void *ptr) throw()
	{ return ::operator new(size, ptr); }
	static void operator delete(void *pMemory, void *ptr) throw()
	{ ::operator delete(pMemory, ptr); }

	// nothrow new/delete
	static void *operator new(std::size_t size, const std::nothrow_t &nt) throw()
	{ return ::operator new(size, nt); }
	static void operator delete(void *pMemory, const std::nothrow_t &) throw()
	{ ::operator delete(pMemory); }

	// debug new/delete
#ifdef _DEBUG
	static void *operator new(std::size_t size, int, const char *, int) throw(std::bad_alloc)
	{ return operator new(size); }
	static void operator delete(void *pMemory, int, const char *, int) throw()
	{ operator delete(pMemory); }
#endif // _DEBUG
};


class NoThrowNewDelete : public StandardNewDelete {
public:
	using StandardNewDelete::operator new;
	using StandardNewDelete::operator delete;

	// nothrow new/delete
	static void *operator new(std::size_t size, const std::nothrow_t &) throw()
	{ return operator new(size); }
	static void operator delete(void *pMemory, const std::nothrow_t &) throw()
	{ operator delete(pMemory); }
};


struct DeletePtr {
	template<typename T>
	void operator()(const T* ptr) const { delete ptr; }
};


struct DeleteAll {
	template<typename T>
	void operator()(const T& container) {
		for_each(container.begin(), container.end(), DeletePtr());
	}
};


struct DeleteAll_Clear {
	template<typename T>
	void operator()(T& container) {
		for_each(container.begin(), container.end(), DeletePtr());
		container.clear();
	}
};


template<typename T>
ostream &operator<<(ostream &os, const vector<T> &v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " ")); 
    return os;
};


// functions for string

void string_replace(string &str, const string &src, const string &dst);
string int2str(int num);
int str2int(const string &str);


#endif // __UTILS_H
