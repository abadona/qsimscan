
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

#define __rerror_cpp__
#include "rerror.h"
#include <time.h>

ErrorStream ers;

void ErrorStream:: throw_ (Rerror& exception)
{
    if (exception.msg_.length ())
        exception.msg_ += " : ";
    exception.msg_ += str ();
    str ("");
    throw exception;
}

static const char TMBUF_SZ = 64;

std::ostream& operator << (std::ostream& e, Throw__ t)
{
    time_t rawtime;
    struct tm* timeinfo;
    char tmbuffer [TMBUF_SZ];

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (tmbuffer, TMBUF_SZ, "[%x %X %Z]", timeinfo);

    if (t.fname_) e << " (module " << t.fname_ << ", line " << t.lno_ << ")";
    e << tmbuffer;
    ((ErrorStream&) e).throw_ (t.exception_);
    return e;
}

std::ostream& operator << (std::ostream& o, const Rerror& rerror)
{
    o << (const char*) rerror;
    return o;
}
