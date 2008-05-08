
//////////////////////////////////////////////////////////////////////////////
// This software module is developed by SCIDM team in 1999-2008.
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
