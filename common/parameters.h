
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
// For any questions please contact SciDM team by email at team@scidm.org
//////////////////////////////////////////////////////////////////////////////

#ifndef __parameters_h__
#define __parameters_h__

#pragma warning (disable: 4786)

#include "platform.h"
#include <map>
#include <ostream>
#include "parameters_section.h"

#ifndef __parameters_cpp__
extern const char* volatile_section_name;
#endif


class Parameters
{
public:
    typedef std::map<std::string, ParametersSection> sectmap;
    sectmap sections_;

    void addSection (const char* sectname, const char* sectdescr, Parameter_descr* pars, int par_no);
    void addParameter (const char* sectname, Parameter_descr& descr);
    void writeHelp (std::ostream& o, const char* header = NULL);

    bool readFile (const char* filename);
    bool writeSection (const char* section, std::ostream& o);
    bool read (std::istream& istr);
    bool write (std::ostream& ostr);

    // parameters access
    bool hasSection (const char* sectname);
    bool hasParameter (const char* section, const char* name);
    bool hasValue (const char* section, const char* name);
    bool setParameter (const char* section, const char* name, const char* value);
    const char* getParameter (const char* section, const char* name);
    bool removeParameter (const char* section, const char* name);
    const char* getDefault (const char* section, const char* name);

    // utility
    longlong getInteger (const char* section, const char* name);
    double  getFloat (const char* section, const char* name);
    bool getBoolean (const char* section, const char* name);
    void setInteger (const char* section, const char* name	, longlong value);
    void setFloat  (const char* section, const char* name, double value);
    void setBoolean (const char* section, const char* name, bool value);
};

std::ostream& operator << (std::ostream&, Parameters&);

#endif // __parameters_h__
