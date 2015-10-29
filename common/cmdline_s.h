//////////////////////////////////////////////////////////////////////////////
//// This software module is developed by SciDM (Scientific Data Management) in 1998-2015
//// 
//// This program is free software; you can redistribute, reuse,
//// or modify it with no restriction, under the terms of the MIT License.
//// 
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//// 
//// For any questions please contact Denis Kaznadzey at dkaznadzey@yahoo.com
//////////////////////////////////////////////////////////////////////////////

#ifndef __cmdline_s__h__
#define __cmdline_s__h__

#pragma warning (disable: 4786)
#pragma warning (disable : 4503)


#include <string>
#include <vector>
#include <map>

typedef std::vector<std::string> Arglist;
typedef std::map<std::string, Arglist> Optdict;

void get_opt (int argc, char** argv, const char* optdef, Arglist& arglist, Optdict& optdict, char** longopts = NULL);
void parse_options (int argc, char** argv, std::string& progname, Arglist& arglist, Optdict& optdict, const char* optdef = "s:h", char** longopts = NULL);


#endif
