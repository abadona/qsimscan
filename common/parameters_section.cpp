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

#include <iterator>
#include "parameters_section.h"

bool ParametersSection::writeHelp (std::ostream& ostr)
{
    ostr << "    Section [" << name_.c_str () << "] : " << description_.c_str () << std::endl;
    for (parmap::iterator itr = parameters_.begin (); itr != parameters_.end (); itr ++)
        ostr << "      " << itr->second.name_.c_str () << " : type = " << itr->second.type_.c_str () << ", default = " << itr->second.def_value_.c_str () << " : " << itr->second.description_.c_str () << std::endl;
    return true;
}

bool ParametersSection::hasParameter (const char* name)
{
    parmap::iterator itr = parameters_.find (name);
    if (itr == parameters_.end ()) return false;
    else return true;
}

