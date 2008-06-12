
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

#pragma warning (disable : 4786)

#ifndef __rerror_h__
#define __rerror_h__

#include <string>
#include <sstream>
#include "i64out.h"

class Rerror
{
protected:
    std::string msg_;

public:
    Rerror (const char* s = "") : msg_ ((const char*) s) {}
    operator const char* () const {return msg_.c_str ();}
friend class ErrorStream;
};

class Throw__
{
public:
    const char* fname_;
    int lno_;
    Rerror exception_;
    Throw__ (Rerror exception, const char* fname = NULL, int lno = 0)
    :
    exception_ (exception),
    fname_ (fname),
    lno_ (lno)
    {
    }
};

class ErrorStream : public std::ostringstream
{
public:
    void throw_ (Rerror& exception);
};

std::ostream& operator << (std::ostream& e, Throw__);
std::ostream& operator << (std::ostream&, const Rerror& rerror);


#ifndef __rerror_cpp__
extern ErrorStream ers;
extern Throw__ Throw_;
#endif

#define Throw Throw__(Rerror (),__FILE__, __LINE__)
#define ThrowEx(X) Throw__(X (),__FILE__, __LINE__)
#define ERRINFO " (module " << __FILE__ << ", line " << __LINE__ << ") "
#define ERR(x) ers << x << Throw
#define FATAL ERR
#define Error(C) ers << ThrowEx (C)
#if !defined (NOCONSTEST)
    #define CONSISTENCY_TEST(x) { if (!(x)) ers << "Consistency test failed for (" #x ")" << Throw; }
#else
    #define CONSISTENCY_TEST(x)
#endif


// standard error strings
#include "common_errors.h"


#endif // __rerror_h__
