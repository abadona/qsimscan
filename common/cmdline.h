
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

#ifndef __cmdline_h__
#define __cmdline_h__

//#pragma warning (disable:4786)

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

typedef std::vector <std::string> svec;

struct KeyFormat
{
    svec shortopts_;
    svec longopts_;
    std::string name_;
    std::string section_;
    std::string parameter_;
    bool        optional_;
    bool        has_arg_;
    std::string arg_type_;
    std::string def_value_;
    std::string description_;

    KeyFormat (const char* shortopts, const char** longopts, const char* name, const char* section, const char* parameter, bool optional, bool has_arg, const char* arg_type, const char* def_value, const char* description);
};
typedef std::vector<KeyFormat> KeysFormat;

struct ArgFormat
{
    std::string name_;
    std::string type_;
    std::string description_;
    bool repeatable_;
    bool optional_;
    ArgFormat (const char* name, const char* type, const char* description, bool repeatable = false, bool optional = false)
        :
        name_ (name),
        type_ (type),
        description_ (description),
        repeatable_  (repeatable),
        optional_ (optional)
    {
    }
};
typedef std::vector<ArgFormat> ArgsFormat;

class CmdLine
{
    KeysFormat& keys_format_;
    ArgsFormat& args_format_;
    std::string error_report_;
    int first_optional_pos_;
    int last_optional_pos_;
    int repeatable_pos_;
    bool ok_;

    void validate_args_format (bool strict = true);
    void parse (int argc, char* argv[], bool strict = true);

public:
    typedef std::vector<std::string> svec;
    svec arguments_;
    typedef std::map<std::string, std::string> ssmap;
    ssmap keys_;

    CmdLine (KeysFormat& keys, ArgsFormat& args, int argc, char* argv [], bool strict = true);

    bool isOk ();
    void reportErrors (std::ostream& o);
    void printHelp (const char* progname, std::ostream& o, bool longform = false);
    void printArgs (std::ostream& o);

    bool hasKey (const char* name);
    const char* getValue (const char* name);
    int getFmtPos (int argno) const;
    int getArgno ();
    const char* getArg (int no);

    KeyFormat* keyFormat (const char* name);
    ArgFormat* argFormat (unsigned argno);
};


#endif
