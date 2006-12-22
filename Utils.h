
#ifndef __UTILS_H
#define __UTILS_H


#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#include "List.h"
#include "Tree.h"
#include "Matrix.h"


class Utils {
public:
	template <class T> static int rand_int(T max) { return (int) ((double) rand() * max / RAND_MAX); }

	template <class T> static T min (T i, T j) { return ((i < j) ? i : j); }
	template <class T> static T max (T i, T j) { return ((i > j) ? i : j); }

	template <class T> static void swap(T &i, T &j) { T temp = i; i = j; j = temp; }

	static void quick_sort_min(int *ind, double *val, int size);
	static void quick_sort_max(int *ind, double *val, int size);

	static unsigned long hash(const char *key) { return hash((const unsigned char *) key, strlen(key)); }
	static unsigned long hash(register const unsigned char *key, register unsigned long len, register unsigned long init = 0);
	static unsigned long hash_size(int n) { return (unsigned long) 1 << n; }
	static unsigned long hash_mask(int n) { return hash_size(n)-1; }
};


typedef List<char, unsigned long> _ConstStringList;
class ConstStringList: public _ConstStringList {
public:
	ConstStringList() { m_is_copy = true; }
	bool contain(const char *str) { return (get(str) != NULL); }
	const char *get(const char *str) { return _ConstStringList::get(Utils::hash(str)); }
	int addFirst(const char *str) { return _ConstStringList::addFirst((char *)str, Utils::hash(str)); }
	int addLast(const char *str) { return _ConstStringList::addLast((char *)str, Utils::hash(str)); }
	int addNext(const char *str) { return _ConstStringList::addNext((char *)str, Utils::hash(str)); }
	int addPrev(const char *str) { return _ConstStringList::addPrev((char *)str, Utils::hash(str)); }
	int addAscent(const char *str) { return _ConstStringList::addAscent((char *)str, Utils::hash(str)); }
	int addDescent(const char *str) { return _ConstStringList::addDescent((char *)str, Utils::hash(str)); }
	int remove(const char *str) { return _ConstStringList::remove(Utils::hash(str)); }
	int release(const char *str) { return _ConstStringList::release(Utils::hash(str)); }
};


class _String {
	char *m_buffer;

public:
	_String(char *buf) { m_buffer = buf; }
	~_String() { delete[] m_buffer; m_buffer = NULL; }

	static char *get(_String *s) { return (s == NULL ? NULL : s->m_buffer); }

	friend class StringList;
};

typedef List<_String, unsigned long> _StringList;
class StringList: public _StringList {
public:
	bool contain(const char *str) { return (get(str) != NULL); }
	char *get(const char *str) { return get(Utils::hash(str)); }
	int addFirst(char *str) { return _StringList::addFirst(new _String(str), Utils::hash(str)); }
	int addLast(char *str) { return _StringList::addLast(new _String(str), Utils::hash(str)); }
	int addNext(char *str) { return _StringList::addNext(new _String(str), Utils::hash(str)); }
	int addPrev(char *str) { return _StringList::addPrev(new _String(str), Utils::hash(str)); }
	int addAscent(char *str) { return _StringList::addAscent(new _String(str), Utils::hash(str)); }
	int addDescent(char *str) { return _StringList::addDescent(new _String(str), Utils::hash(str)); }
	int remove(const char *str) { return _StringList::remove(Utils::hash(str)); }
	int release(const char *str) { return _StringList::release(Utils::hash(str)); }
	char *getFirst() { return _String::get(_StringList::getFirst()); }
	char *getLast() { return _String::get(_StringList::getLast()); }
	char *getNext() { return _String::get(_StringList::getNext()); }
	char *getPrev() { return _String::get(_StringList::getPrev()); }
	char *getCurrent() const { return _String::get(_StringList::getCurrent()); }
	char *getItem(ListItem<_String, unsigned long> *item) { return _String::get(_StringList::getItem(item)); }
	char *get(unsigned long key) { return _String::get(_StringList::get(key)); }
};


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


class ArgOption;
typedef List<ArgOption, unsigned long> ArgOptionList;


class ArgOption {
	char *m_name;
	char *m_short_name;
	bool m_require_value;
	char *m_default_value;
	ConstStringList m_possible_values;
	ArgOptionList m_conflict_options;
	char *m_help_info;
	bool m_defined;
	char *m_value;

public:
	ArgOption(char *name, char *short_name = NULL, bool require_value = false, char *default_value = NULL, int possible_values_num = 0, ...);
	ArgOption(char *name, char *short_name = NULL, bool require_value = false, char *default_value = NULL, int possible_values_num = 0, va_list argptr = NULL);
	~ArgOption();

	const char *name() { return m_name; }
	const char *short_name() { return m_short_name; }
	bool require_value() { return m_require_value; }
	const char *default_value() { return m_default_value; }
	const char *help_info() { return m_help_info; }

	void setValue(char *value = NULL);
	void addConflictOption(ArgOption *option);
	void addHelpInfo(char *info);

	char *getValue() const;
	int getValueAsInt() const;
	double getValueAsDouble() const;

	bool isDefined() { return m_defined; }
	bool isConflicted();
	bool equal(const char *value) const;

	char *printDef(char *buffer);
	char *printVal(char *buffer);

	friend class ArgParser;
};


class ArgParser {
	const char *m_program_name;
	ArgOptionList m_options_with_value;
	ArgOptionList m_options_without_value;
	ArgOptionList m_options_shortname;
	ArgOptionList m_options_default;
	ArgOptionList m_args_option;
	ConstStringList m_args_non_option;
	
public:
	ArgParser();
	~ArgParser();

	int getOptionArgNum() const { return m_args_option.size(); }
	int getNonOptionArgNum() const { return m_args_non_option.size(); }

	const char *getNonOptionArgument(int i = 0);

	ArgOption *findOption(char *name);
	ArgOption *getOption(char *name);

	ArgOption *addOption(char *name, char *short_name = NULL, bool require_value = false, char *default_value = NULL, int possible_values_num = 0, ...);
	void addHelpInfo(char *name, char *info) { getOption(name)->addHelpInfo(info); }
	void addConflictOptions(char *name1, char *name2);
	void checkConflictOptions(ArgOption *option);

	void parseArguments(int argc, char *argv[]);

	bool isDefined(char *name) { return getOption(name)->isDefined(); }

	char *printDef(char *buffer);
	char *printVal(char *buffer);
};


#endif // __UTILS_H
