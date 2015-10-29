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

#ifndef __parameters_section_h__
#define __parameters_section_h__
#pragma warning (disable:4786)

#include <string>
#include <map>
#include <ostream>

struct Parameter_descr
{
    const char* name_;
    const char* type_;
    const char* def_value_;
    const char* description_;
};


struct Parameter
{
    Parameter () {}
    Parameter (Parameter_descr& descr)
    :
    name_ (descr.name_),
    type_ (descr.type_),
    def_value_ (descr.def_value_),
    description_ (descr.description_)
    {
    }
    std::string name_;
    std::string type_;
    std::string def_value_;
    std::string description_;
    std::string value_;
};


struct ParametersSection
{
    typedef std::map<std::string, Parameter> parmap;

    std::string name_;
    std::string description_;
    parmap parameters_;
    bool writeHelp (std::ostream& ostr);
    bool hasParameter (const char* name);
};


#endif //__parameters_section_h__
