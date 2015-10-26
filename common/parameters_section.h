
//////////////////////////////////////////////////////////////////////////////
// This software module is developed by SCIDM team in 1999-2015.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// For any questions please contact SciDM team by email at team@scidm.org
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
