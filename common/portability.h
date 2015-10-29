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

#ifndef __portability_h__
#define __portability_h__

// long file support flavours for different platforms
#if defined (_MSC_VER)
    #include <io.h>
    #define sci_stat _stati64
    #define sci_stat_struc _stati64
    #define sci_open _open
    #define sci_read _read
    #define sci_write _write
    #define sci_close _close
    #define sci_lseek _lseeki64
    #define sci_tell _telli64
#elif defined (__CYGWIN__)
    #include <unistd.h>
    #define sci_stat stat
    #define sci_stat_struc stat
    #define sci_open open
    #define sci_read read
    #define sci_write write
    #define sci_close close
    #define sci_lseek lseek
    #define sci_tell tell
#elif defined (__APPLE__)
    #include <unistd.h>
    #define sci_stat stat
    #define sci_stat_struc stat
    #define sci_open open
    #define sci_read read
    #define sci_write write
    #define sci_close close
    #define sci_lseek lseek
    #define sci_tell tell
#else // plain unix :)
    #include <unistd.h>
    #define sci_stat stat64
    #define sci_stat_struc stat64
    #define sci_open open64
    #define sci_read read
    #define sci_write write
    #define sci_close close
    #define sci_lseek lseek64
    #define sci_tell tell64
#endif

#if defined (_MSC_VER)
    #define strncasecmp strnicmp
    #define strcasecmp stricmp
#endif

#ifndef O_BINARY
#define O_BINARY 0
#endif

#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif


#endif // __portability_h__
