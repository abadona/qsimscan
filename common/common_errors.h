
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

#ifndef __common_errors_h__
#define __common_errors_h__

#ifndef __common_errors_cpp__

extern const char* ERR_NoMemory;
extern const char* ERR_Internal;
extern const char* ERR_FileNotFound;
extern const char* ERR_OSError;
extern const char* ERR_OutOfBounds;


// synonyms
#define NOEMEM ERR_NoMemory
#define INTERNAL ERR_Internal

#define MAKE_ERROR_TYPE(C,N) class C : public Rerror\
{\
public:\
    C (const char* s = "")\
    : Rerror (N)\
    {\
        if (*s)\
        {\
            msg_ += ": ";\
            msg_ += s;\
        }\
    }\
};

MAKE_ERROR_TYPE (MemoryRerror, ERR_NoMemory);
MAKE_ERROR_TYPE (InternalRerror, ERR_Internal);
MAKE_ERROR_TYPE (FileNotFoundRerror, ERR_FileNotFound);
MAKE_ERROR_TYPE (OutOfBoundsRerror, ERR_OutOfBounds);

#endif

class OSRerror : public Rerror
{
    const char* get_err_str () const;
    const char* get_errno_str () const;
public:
    OSRerror (const char* s = "")
    : Rerror (ERR_OSError)
    {
        msg_ += ": errno ";
        msg_ += get_errno_str ();
        msg_ += ": ";
        msg_ += get_err_str ();
        if (*s)
        {
            msg_ += ": ";
            msg_ += s;
        }
    }
};


#endif // __common_errors_h__
