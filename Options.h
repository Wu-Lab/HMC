
#ifndef __OPTIONS_H
#define __OPTIONS_H


#include <iosfwd>
#include <boost/program_options.hpp>

#include "Utils.h"


namespace po = ::boost::program_options;


struct DisplayOption {
	enum { all, defaulted, specified };

	ostream &os;
	int flag;

	DisplayOption(ostream &s, int f = 0) : os(s), flag(f) { }
	void operator()(const pair<string, po::variable_value> &opt) const;
};


void print_options(po::variables_map &vm, ostream &os, int specified = DisplayOption::all);


/* Function used to check that 'opt1' and 'opt2' are not specified
   at the same time. */
void conflicting_options(const po::variables_map& vm, 
                         const char* opt1, const char* opt2);


/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const po::variables_map& vm,
                       const char* for_what, const char* required_option);


#endif // __OPTIONS_H