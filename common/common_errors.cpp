
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


#define __common_errors_cpp__

const char* ERR_NoMemory = "Unable to allocate memory";
const char* ERR_Internal = "Internal program error";
const char* ERR_FileNotFound = "File not found";
const char* ERR_OSError = "Operating system error";
const char* ERR_OutOfBounds = "Out of bounds error";

#include <cerrno>
#include <cstring>
#include <cstdio>
#include "rerror.h"

// Note that "common_errors.h" is nt included - rerror.h includes it at the end

const char* OSRerror::get_err_str () const
{
    return strerror (errno);
}

static char errno_str_buf [8];
const char* OSRerror::get_errno_str () const
{
    sprintf (errno_str_buf, "%d", errno);
    return errno_str_buf;
}

