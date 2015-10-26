
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
// For any questions please contact SciDM team by email at scidmteam@yahoo.com
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

