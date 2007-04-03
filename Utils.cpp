
#include <cstdlib>
#include <cstdarg>

#include "Utils.h"


////////////////////////////////
//
// class Timer


int Timer::m_instance_num = 0;

void Timer::begin(const char *description)
{
	m_running = true;
	m_duration = 0;
	strcpy(m_description, description);
	Logger::debug("TIMER %d: %s ...", m_id, m_description);
	m_start = clock();
}

void Timer::end()
{
	pause();
	Logger::debug("TIMER %d: %s = %f seconds", m_id, m_description, m_duration);
}

void Timer::pause()
{
	if (m_running) {
		m_finish = clock();
		m_duration += (double) (m_finish - m_start) / CLOCKS_PER_SEC;
		m_running = false;
	}
}

void Timer::resume()
{
	if (!m_running) {
		m_start = clock();
		m_running = true;
	}
}


////////////////////////////////
//
// class Logger


bool Logger::m_logging = true;
int Logger::m_log_level = 3;
const int Logger::m_log_level_error = 1;
const int Logger::m_log_level_warning = 2;
const int Logger::m_log_level_info = 3;
const int Logger::m_log_level_verbose = 4;
const int Logger::m_log_level_debug = 5;
Timer Logger::m_timer[10];

void Logger::setLogLevel(int level)
{
	m_log_level = level;
}

void Logger::enableLogging()
{
	m_logging = true;
}

void Logger::disableLogging()
{
	m_logging = false;
}

void Logger::error(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_error, stderr, "ERROR: ", format, argptr);
	va_end(argptr);
}

void Logger::warning(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_warning, stderr, "WARNING: ", format, argptr);
	va_end(argptr);
}

void Logger::info(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_info, stdout, NULL, format, argptr);
	va_end(argptr);
}

void Logger::verbose(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_verbose, stderr, NULL, format, argptr);
	va_end(argptr);
}

void Logger::status(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_print(m_log_level_verbose, stderr, NULL, format, argptr);
	print(m_log_level_verbose, stderr, NULL, "               \r");
	va_end(argptr);
}

void Logger::debug(const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(m_log_level_debug, stderr, NULL, format, argptr);
	va_end(argptr);
}

void Logger::print(int level, FILE *fp, const char *prompt, const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_print(level, stderr, NULL, format, argptr);
	va_end(argptr);
}

void Logger::println(int level, FILE *fp, const char *prompt, const char *format, ...)
{
	va_list argptr;
	va_start(argptr, format);
	_println(level, stderr, NULL, format, argptr);
	va_end(argptr);
}

inline void Logger::_print(int level, FILE *fp, const char *prompt, const char *format, va_list argptr)
{
	if (m_log_level >= level) {
		if (prompt != NULL) fprintf(fp, prompt);
		vfprintf(fp, format, argptr);
	}
}

inline void Logger::_println(int level, FILE *fp, const char *prompt, const char *format, va_list argptr)
{
	if (m_log_level >= level) {
		if (prompt != NULL) fprintf(fp, prompt);
		vfprintf(fp, format, argptr);
		fprintf(fp, "\n");
	}
}

void Logger::beginTimer(int i, const char *description)
{
	m_timer[i].begin(description);
}

void Logger::endTimer(int i)
{
	m_timer[i].end();
}

void Logger::pauseTimer(int i)
{
	m_timer[i].pause();
}

void Logger::resumeTimer(int i)
{
	m_timer[i].resume();
}


////////////////////////////////
//
// functions for string

void string_replace(string &str, const string &src, const string &dst)
{
	string::size_type npos = static_cast<string::size_type>(string::npos);
	string::size_type pos = 0;
	string::size_type srclen = src.size();
	string::size_type dstlen = dst.size();
	while((pos=str.find(src, pos)) != npos) {
		str.replace(pos, srclen, dst);
		pos += dstlen;
	}
}

string int2str(int num)
{
	if (num == 0) return "0";
	string str = "";
	int num_abs = num > 0 ? num : -num;
	while (num_abs) {
		str = static_cast<char>(num_abs % 10 + '0') + str;
		num_abs /= 10;
	}
	if (num < 0) str = "-" + str;
	return str;
}

int str2int(const string &str)
{
	int len = str.size(), num = 0, i = 0;
    if (str[0] == '-') i = 1;
	while (i < len) {
		num = num * 10 + static_cast<int>(str[i] - '0');
		++i;
	}
	if (str[0] == '-') num *= -1;
	return num;
}
