
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "Utils.h"


////////////////////////////////
//
// class Utils


/* quick sort function */

void Utils::quick_sort_min(int *ind, double *val, int size)
{
	int i, j;
	int t;
	double d;

	double pivot;

	if (size <= 24) {

		/* insertion sort */
		for (i=1;i<size;i++) {
			t = ind[i];
			d = val[i];
			for (j=i-1;j>=0;j--) {
				if (val[j] > d) {
					ind[j+1] = ind[j];
					val[j+1] = val[j];
				}
				else {
					break;
				}
			}
			ind[j+1] = t;
			val[j+1] = d;
		}
	}
	else {
		pivot = val[rand_int(size - 1)];

		i = -1;
		j = size;
		while (true) {
			do {
				i++;
			} while (val[i] < pivot);
			do {
				j--;
			} while (val[j] > pivot);
			if (i < j) {
				t = ind[i];
				d = val[i];
				ind[i] = ind[j];
				val[i] = val[j];
				ind[j] = t;
				val[j] = d;
			}
			else {
				if (i == 0) {
					i = 1;
				}
				break;
			}
		}

		quick_sort_min(ind, val, i);
		quick_sort_min(ind + i, val + i, size - i);
	}
}

void Utils::quick_sort_max(int *ind, double *val, int size)
{
	int i, j;
	int t;
	double d;

	double pivot;

	if (size <= 24) {

		/* insertion sort */
		for (i=1;i<size;i++) {
			t = ind[i];
			d = val[i];
			for (j=i-1;j>=0;j--) {
				if (val[j] < d) {
					ind[j+1] = ind[j];
					val[j+1] = val[j];
				}
				else {
					break;
				}
			}
			ind[j+1] = t;
			val[j+1] = d;
		}
	}
	else {
		pivot = val[rand_int(size - 1)];

		i = -1;
		j = size;
		while (true) {
			do {
				i++;
			} while (val[i] > pivot);
			do {
				j--;
			} while (val[j] < pivot);
			if (i < j) {
				t = ind[i];
				d = val[i];
				ind[i] = ind[j];
				val[i] = val[j];
				ind[j] = t;
				val[j] = d;
			}
			else {
				if (i == 0) {
					i = 1;
				}
				break;
			}
		}

		quick_sort_max(ind, val, i);
		quick_sort_max(ind + i, val + i, size - i);
	}
}

/* hash function */

#define hash_mix(a, b, c) \
{ \
	a -= b; a -= c; a ^= (c>>13); \
    b -= c; b -= a; b ^= (a<<8);  \
	c -= a; c -= b; c ^= (b>>13); \
	a -= b; a -= c; a ^= (c>>12); \
	b -= c; b -= a; b ^= (a<<16); \
	c -= a; c -= b; c ^= (b>>5);  \
	a -= b; a -= c; a ^= (c>>3);  \
	b -= c; b -= a; b ^= (a<<10); \
	c -= a; c -= b; c ^= (b>>15); \
}

unsigned long Utils::hash(register const unsigned char *key, register unsigned long len, register unsigned long init)
{
	register unsigned long a, b, c, l;
	/* set up the internal state */
	l = len;
	a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
	c = init;         /* the previous hash value */
	/* handle most of the key */
	while (l >= 12)
	{
		a += (key[0] + ((unsigned long) key[1]<<8) + ((unsigned long) key[2]<<16) + ((unsigned long) key[3]<<24));
		b += (key[4] + ((unsigned long) key[5]<<8) + ((unsigned long) key[6]<<16) + ((unsigned long) key[7]<<24));
		c += (key[8] + ((unsigned long) key[9]<<8) + ((unsigned long) key[10]<<16) + ((unsigned long) key[11]<<24));
		hash_mix(a, b, c);
		key += 12;
		l -= 12;
	}
	/* handle the last 11 bytes */
	c += len;
	switch(l)	/* all the case statements fall through */
	{
		case 11: c += ((unsigned long) key[10]<<24);
		case 10: c += ((unsigned long) key[9]<<16);
		case 9 : c += ((unsigned long) key[8]<<8);
		/* the first byte of c is reserved for the length */
		case 8 : b += ((unsigned long) key[7]<<24);
		case 7 : b += ((unsigned long) key[6]<<16);
		case 6 : b += ((unsigned long) key[5]<<8);
		case 5 : b += key[4];
		case 4 : a += ((unsigned long) key[3]<<24);
		case 3 : a += ((unsigned long) key[2]<<16);
		case 2 : a += ((unsigned long) key[1]<<8);
		case 1 : a += key[0];
		/* case 0: nothing left to add */
	}
	hash_mix(a,b,c);
	return c;
}


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
// class ArgList


ArgOption::ArgOption(char *name, char *short_name, bool require_value, char *default_value, int possible_values_num, ...)
{
	va_list possible_values;
	va_start(possible_values, possible_values_num);
	ArgOption(name, short_name, require_value, default_value, possible_values_num, possible_values);
	va_end(possible_values);
}

ArgOption::ArgOption(char *name, char *short_name, bool require_value, char *default_value, int possible_values_num, va_list possible_values)
{
	int i;
	char *value;
	m_name = name;
	m_short_name = short_name;
	m_require_value = require_value;
	m_default_value = default_value;
	for (i=0; i<possible_values_num; i++) {
		value = va_arg(possible_values, char *);
		m_possible_values.addLast(value);
	}
	m_defined = false;
	m_value = m_default_value;
	m_help_info = NULL;
}

ArgOption::~ArgOption()
{
	m_conflict_options.releaseAll();
}

void ArgOption::setValue(char *value)
{
	if (m_require_value) {
		if (m_possible_values.size() > 0) {
			if (m_possible_values.get(value) != NULL) {
				m_defined = true;
				m_value = value;
			}
			else {
				Logger::error("Invalid value for option -%s: %s!", m_name, value);
				exit(1);
			}
		}
		else {
			m_defined = true;
			m_value = value;
		}
	}
	else {
		m_defined = true;
		m_value = m_default_value;
	}
}

void ArgOption::addConflictOption(ArgOption *option)
{
	unsigned long key = Utils::hash(option->m_name);
	if (m_conflict_options.get(key) == NULL) {
		m_conflict_options.addLast(option, key);
	}
}

void ArgOption::addHelpInfo(char *info)
{
	m_help_info = info;
}

char *ArgOption::getValue() const
{
	if (m_value != NULL) {
		return m_value;
	}
	return m_default_value;
}

int ArgOption::getValueAsInt() const
{
	if (m_value != NULL) {
		return atoi(m_value);
	}
	else if (m_default_value != NULL) {
		return atoi(m_default_value);
	}
	return 0;
}

double ArgOption::getValueAsDouble() const
{
	if (m_value != NULL) {
		return atof(m_value);
	}
	else if (m_default_value != NULL) {
		return atof(m_default_value);
	}
	return 0;
}

bool ArgOption::isConflicted()
{
	ArgOption *option = m_conflict_options.getFirst();
	while (option != NULL) {
		if (option->isDefined()) {
			return true;
		}
		option = m_conflict_options.getNext();
	}
	return false;
}

bool ArgOption::equal(const char *value) const
{
	if (m_value != NULL) {
		return (strcmp(value, m_value) == 0);
	}
	else if (m_default_value != NULL) {
		return (strcmp(value, m_default_value) == 0);
	}
	return false;
}

char *ArgOption::printDef(char *buffer)
{
	const char *s;
	strcpy(buffer, "-");
	strcat(buffer, m_name);
	if (m_short_name != NULL) {
		strcat(buffer, " (-");
		strcat(buffer, m_short_name);
		strcat(buffer, ")");
	}
	if (m_possible_values.size() > 0) {
		strcat(buffer, " {");
		s = m_possible_values.getFirst();
		while (s != NULL) {
			strcat(buffer, s);
			strcat(buffer, "|");
			s = m_possible_values.getNext();
		}
		buffer[strlen(buffer)-1] = '}';
	}
	if (m_default_value != NULL) {
		strcat(buffer, ", default=");
		strcat(buffer, m_default_value);
	}
	if (m_help_info != NULL) {
		strcat(buffer, "\n\t");
		strcat(buffer, m_help_info);
	}
	return buffer;
}

char *ArgOption::printVal(char *buffer)
{
	strcpy(buffer, "-");
	strcat(buffer, m_name);
	if (m_value != NULL) {
		strcat(buffer, " = ");
		strcat(buffer, m_value);
	}
	return buffer;
}

ArgParser::ArgParser()
{
	m_program_name = NULL;
}

ArgParser::~ArgParser()
{
	m_options_shortname.releaseAll();
	m_options_default.releaseAll();
	m_args_option.releaseAll();
	m_args_non_option.releaseAll();
}

const char *ArgParser::getNonOptionArgument(int i)
{
	if (i >= 0 && i < m_args_non_option.size()) {
		m_args_non_option.getFirst();
		while (i > 0) {
			m_args_non_option.getNext();
			i--;
		}
		return m_args_non_option.getCurrent();
	}
	else {
		return NULL;
	}
}

ArgOption *ArgParser::findOption(char *name)
{
	ArgOption *arg;
	unsigned long key = Utils::hash(name);
	arg = m_options_with_value.get(key);
	if (arg != NULL) return arg;
	arg = m_options_without_value.get(key);
	if (arg != NULL) return arg;
	arg = m_options_shortname.get(key);
	if (arg != NULL) return arg;
	return NULL;
}

ArgOption *ArgParser::getOption(char *name)
{
	ArgOption *arg = findOption(name);
	if (arg == NULL) {
		Logger::error("Invalid option -%s!", name);
		exit(1);
	}
	return arg;
}

ArgOption *ArgParser::addOption(char *name, char *short_name, bool require_value, char *default_value, int possible_values_num, ...)
{
	ArgOption *arg;
	va_list possible_values;
	va_start(possible_values, possible_values_num);
	arg = new ArgOption(name, short_name, require_value, default_value, possible_values_num, possible_values);
	va_end(possible_values);
	if (findOption(name) != NULL) {
		Logger::error("[ArgList::addOption] Duplicated option name -%s!", name);
		exit(1);
	}
	if (short_name != NULL && findOption(short_name) != NULL) {
		Logger::error("[ArgList::addOption] Duplicated option short name -%s!", short_name);
		exit(1);
	}
	if (require_value) {
		m_options_with_value.addLast(arg, Utils::hash(name));
	}
	else {
		m_options_without_value.addLast(arg, Utils::hash(name));
	}
	if (short_name != NULL) {
		m_options_shortname.addLast(arg, Utils::hash(short_name));
	}
	if (default_value != NULL) {
		m_options_default.addLast(arg, Utils::hash(name));
	}
	return arg;
}

void ArgParser::addConflictOptions(char *name1, char *name2)
{
	ArgOption *option1 = getOption(name1);
	ArgOption *option2 = getOption(name2);
	if (option1 != NULL && option2 != NULL) {
		option1->addConflictOption(option2);
		option2->addConflictOption(option1);
	}
}

void ArgParser::checkConflictOptions(ArgOption *option)
{
	ArgOption *conflict = option->m_conflict_options.getFirst();
	while (conflict != NULL) {
		if (conflict->isDefined()) {
			Logger::error("The options -%s and -%s can not be used together!", option->m_name, conflict->m_name);
			exit(1);
		}
		m_options_default.release(Utils::hash(conflict->m_name));
		conflict = option->m_conflict_options.getNext();
	}
}

void ArgParser::parseArguments(int argc, char *argv[])
{
	int i;
	char *arg;
	ArgOption *option;
	m_program_name = argv[0];
	for (i=1; i<argc; i++) {
		arg = argv[i];
		if (arg[0] != '-') {
			m_args_non_option.addLast(arg);
		}
		else {
			arg++;
			option = getOption(arg);
			checkConflictOptions(option);
			if (option->require_value()) {
				i++;
				if (i < argc && argv[i][0] != '-') {
					option->setValue(argv[i]);
				}
				else {
					Logger::error("Missing value for option -%s!", option->m_name);
					exit(1);
				}
			}
			else {
				option->setValue();
			}
			m_options_default.release(Utils::hash(option->m_name));
			m_args_option.addLast(option, Utils::hash(option->m_name));
		}
	}
}

char *ArgParser::printDef(char *buffer)
{
	char buf[128];
	ArgOption *option;
	strcpy(buffer, "");
	if (m_options_with_value.size() > 0) {
		strcat(buffer, "Options with value:\n");
		option = m_options_with_value.getFirst();
		while (option != NULL) {
			strcat(buffer, "  ");
			strcat(buffer, option->printDef(buf));
			strcat(buffer, "\n");
			option = m_options_with_value.getNext();
		}
	}
	if (m_options_without_value.size() > 0) {
		if (strlen(buffer) > 0) {
			strcat(buffer, "\n");
		}
		strcat(buffer, "Options without value:\n");
		option = m_options_without_value.getFirst();
		while (option != NULL) {
			strcat(buffer, "  ");
			strcat(buffer, option->printDef(buf));
			strcat(buffer, "\n");
			option = m_options_without_value.getNext();
		}
	}
	return buffer;
}

char *ArgParser::printVal(char *buffer)
{
	char buf[10240];
	const char *s;
	ArgOption *option;
	strcpy(buffer, "");
	if (m_args_option.size() > 0) {
		strcat(buffer, "Options:\n");
		option = m_args_option.getFirst();
		while (option != NULL) {
			strcat(buffer, "  ");
			strcat(buffer, option->printVal(buf));
			strcat(buffer, "\n");
			option = m_args_option.getNext();
		}
	}
	if (m_options_default.size() > 0) {
		if (strlen(buffer) > 0) {
			strcat(buffer, "\n");
		}
		strcat(buffer, "Options with default value:\n");
		option = m_options_default.getFirst();
		while (option != NULL) {
			strcat(buffer, "  ");
			strcat(buffer, option->printVal(buf));
			strcat(buffer, "\n");
			option = m_options_default.getNext();
		}
	}
	if (m_args_non_option.size() > 0) {
		if (strlen(buffer) > 0) {
			strcat(buffer, "\n");
		}
		strcat(buffer, "Non-option arguments:\n");
		s = m_args_non_option.getFirst();
		while (s != NULL) {
			strcat(buffer, "  ");
			strcat(buffer, s);
			strcat(buffer, "\n");
			s = m_args_non_option.getNext();
		}
	}
	return buffer;
}
