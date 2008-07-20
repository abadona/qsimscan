
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

#ifndef __process_params_h__
#define __process_params_h__

#include "cmdline.h"
#include "parameters.h"
#include <ostream>

class InvalidParameter
{
public:
    InvalidParameter (const char* section, const char* name, const char* value, const char* explanation)
    :
    section_ (section),
    name_ (name),
    value_ (value),
    explanation_ (explanation)
    {
    }
    void report (std::ostream& o)
    {
        o << "Invalid parameter " << name_ << " in section " << section_ << " (value " << value_ << "). " << explanation_ << std::endl;
    }

    std::string section_;
    std::string name_;
    std::string value_;
    std::string explanation_;
};


class Process_params
{
    bool templates_prepared_;
    bool parameters_read_;
    bool verbose_;
    bool overwrite_;
    int  debug_level_;

protected:
    const char* header_;
    std::string procname_;
    std::string version_;
    bool help_mode_;

    CmdLine* cmdline_;
    // SCIDM::Session_ptr session_;
    std::string parfname_;
    std::string def_parfname_;

    KeysFormat keys_format_;
    ArgsFormat args_format_;

    static const char *parfname_default_;

    virtual bool interpreteCmdline ();

    // the following functions can be overloaded in derived classes:

    virtual bool createParamsObject ();   // parameters_ object sometimes has to be created by specialized constructor (like Persistent_params one)
    virtual void add_cmdline_srv ();
    virtual void add_cmdline_app ();
    virtual void add_cmdline_overwrite ();
    virtual void add_parameters_overwrite ();

public: // HACK for AppServer (uses interpreteParameters directly. TODO: change AppServer to use derivate?
    // the following overloaded function should BETTER ensure that Process_params's variant of that function is also called
    virtual bool prepareCmdlineFormat (); // fills in keys_format_ list.
    virtual bool prepareParameters ();    // fills in parameters_ list. Must call base class prepareParameters before returning
    virtual bool populateParameters ();   // populate parameters values from any non-cmdline sources (this includes reading cfgfile, database or whatever customization of base parameters). Be carful if not calling base version.
    virtual bool interpreteParameters (); // reads parameters object and assigns specific members of Process_params derivate; must call base version or nase options and arguments will not be in.
    virtual bool validate ();             // validates specific values; should better call base version at the end.
    virtual bool postProcessParams ();    // does any task needed after parameters are fully interpreted, like writing cfg file

public:

    Parameters* parameters_;

    Process_params (const char* header = NULL, const char* procname = NULL, const char* version = NULL);
    virtual ~Process_params ();

    bool help_mode () const;
    CmdLine* cmdline ();

    bool prepareTemplates ();
    bool parseCmdline (int argc, char* argv [], bool strict = true);
    bool process ();
    void writeParams ();
    const char* parfname () const {return parfname_.c_str ();}
    virtual const char* def_parfname () const {return def_parfname_.c_str ();}
    bool parameters_read () {return parameters_read_;}

    // get section
    bool verbose () const {return verbose_;}
    bool overwrite () const {return overwrite_;}
    int  debug () const {return debug_level_;}

    // set section
    void verbose (bool verbosity) {verbose_ = verbosity;}
    void overwrite (bool val) {overwrite_ = val;}
    void debug (int val) {debug_level_ = val;}

    // defaults section
    virtual const char* verbose_default () const;
    virtual const char* overwrite_default () const;
    virtual const char* debug_default () const;

    // help messages section
    virtual const char* verbose_help () const;
    virtual const char* overwrite_help () const;
    virtual const char* debug_help () const;
};


#endif // __process_params_h__