//////////////////////////////////////////////////////////////////////////////
//// This software module is developed by SciDM (Scientific Data Management) in 1998-2015
//// 
//// This program is free software; you can redistribute, reuse,
//// or modify it with no restriction, under the terms of the MIT License.
//// 
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//// 
//// For any questions please contact Denis Kaznadzey at dkaznadzey@yahoo.com
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

