
#include <iomanip>

#include "Options.h"


void DisplayOption::operator()(const pair<string, po::variable_value> &opt) const
{
	if (((flag == DisplayOption::defaulted) && !opt.second.defaulted()) ||	
		((flag == DisplayOption::specified) && opt.second.defaulted())) return;
	os << "  " << opt.first << " = ";
	if (!opt.second.empty()) {
		const std::type_info &type = opt.second.value().type();
		if (type == typeid(int)) {
			os << opt.second.as<int>();
		}
		else if (type == typeid(double)) {
			os << opt.second.as<double>();
		}
		else if (type == typeid(bool)) {
			os << (opt.second.as<bool>() ? "TRUE" : "FALSE");
		}
		else if (type == typeid(string)) {
			os << opt.second.as<string>();
		}
		else if (type == typeid(vector<string>)) {
			os << opt.second.as<vector<string> >();
		}
		else {
			os << " <unknown type> ";
		}
	}
	else {
		os << " <empty> ";
	}
	os << endl;
}


void print_options(po::variables_map &vm, ostream &os, int specified)
{
	for_each(vm.begin(), vm.end(), DisplayOption(os, specified));
}


/* Function used to check that 'opt1' and 'opt2' are not specified
   at the same time. */
void conflicting_options(const po::variables_map& vm, 
                         const char* opt1, const char* opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() 
        && vm.count(opt2) && !vm[opt2].defaulted())
        throw logic_error(string("Conflicting options '") 
                          + opt1 + "' and '" + opt2 + "'.");
}


/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const po::variables_map& vm,
                       const char* for_what, const char* required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted())
            throw logic_error(string("Option '") + for_what 
                              + "' requires option '" + required_option + "'.");
}
