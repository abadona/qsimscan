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

// #pragma warning (disable : 4786)

#ifndef __rerror_h__
#define __rerror_h__

#include <string>
#include <sstream>
#include <ctime>
#include "i64out.h"
#include "common_str.h"

class Rerror
{
protected:
    std::string msg_;

public:
    Rerror (const char* s = EMPTY_STR) : msg_ ((const char*) s) {}
    operator const char* () const {return msg_.c_str ();}
friend class ErrorStream;
};

template <typename ExceptionType = Rerror>
class Throw__
{
public:
    const char* fname_;
    int lno_;
    ExceptionType exception_;
    Throw__ (ExceptionType exception, const char* fname = NULL, int lno = 0)
    :
    fname_ (fname),
    lno_ (lno),
    exception_ (exception)
    {
    }
};

class ErrorStream : public std::ostringstream
{
public:
    template <typename ExceptionType>
    void throw_ (ExceptionType exception)
    {
        if (exception.msg_.length ())
            exception.msg_ += " : ";
        exception.msg_ += str ();
        str (EMPTY_STR);
        throw exception;
    }
};

template <typename ExceptionType>
inline ErrorStream& operator << (ErrorStream& e, Throw__<ExceptionType> t)
{
    time_t rawtime;
    struct tm* timeinfo;
    static const unsigned TMBUF_SZ = 24;
    char tmbuffer [TMBUF_SZ];

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (tmbuffer, TMBUF_SZ, "[%x %X %Z]", timeinfo);

    if (t.fname_) e << " (module " << t.fname_ << ", line " << t.lno_ << ")";
    e << tmbuffer;
    e.throw_ (t.exception_);
    return e;
}

template <typename OpType>
inline ErrorStream& operator << (ErrorStream& e, const OpType& t)
{
    ((std::ostringstream&) e) << t;
    return e;
}

// ostream manipulators support
inline ErrorStream& operator << (ErrorStream& e, std::ostream& (*op) (std::ostream&))
{
    ((std::ostringstream&) e) << op;
    return e;
}

std::ostream& operator << (std::ostream&, const Rerror& rerror);

#ifndef __rerror_cpp__
extern ErrorStream ers;
//extern Throw__ Throw_;
#endif

#define Throw Throw__<Rerror>(Rerror (),__FILE__, __LINE__)
#define ThrowEx(ErrorType) Throw__<ErrorType>(ErrorType (),__FILE__, __LINE__)
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
