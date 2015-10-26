
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

#ifndef __PLATFORM_H__
#define __PLATFORM_H__

#include "compile_time_macro.h"

#ifdef _MSC_VER
typedef __int64 longlong;
typedef unsigned __int64 ulonglong;
#else
typedef long long int longlong;
typedef long long unsigned int ulonglong;
#endif

#if defined (__x86__) || defined (__powerpc__)
#define SCIDM_LITTLE_ENDIAN
//#warning LITTLE ENDIAN
#elif defined (__mips__) || defined (__alpha__) || defined (__rx000__)
#define SCIDM_BIG_ENDIAN
#warning BIG ENDIAN
#else
#error Byte order unknown
#endif

// the following works for MSVC on 32-bit intel and gcc on 32 and 64 bit (LP64 speciication) intel.
// Could fail on other platforms / compilers
typedef unsigned char BYTE;
typedef char SBYTE;
typedef unsigned short WORD;
typedef short SWORD;
typedef unsigned int DWORD;
typedef int SDWORD;
typedef ulonglong QWORD;
typedef longlong SQWORD;

static void dummy ()
{
    ASSERT_EXACT_BITSIZE (BYTE, 8)
    ASSERT_EXACT_BITSIZE (SBYTE, 8)
    ASSERT_EXACT_BITSIZE (WORD, 16)
    ASSERT_EXACT_BITSIZE (SWORD, 16)
    ASSERT_EXACT_BITSIZE (DWORD, 32)
    ASSERT_EXACT_BITSIZE (SDWORD, 32)
    ASSERT_EXACT_BITSIZE (QWORD, 64)
    ASSERT_EXACT_BITSIZE (SQWORD, 64)
}

#ifdef SCIDM_LITTLE_ENDIAN

#define GET32_U(ptr) (*(const DWORD*)(ptr))
#define GET32_UR(ptr) ( (((DWORD) ((const BYTE*) ptr) [0]) << 24) | \
                        (((DWORD) ((const BYTE*) ptr) [1]) << 16) | \
                        (((DWORD) ((const BYTE*) ptr) [2]) << 8 ) | \
                        (((DWORD) ((const BYTE*) ptr) [3])      ) )

#define GET64_U(ptr) (*(const QWORD*)(ptr))
#define GET64_UR(ptr) ( (((QWORD) ((const BYTE*) ptr) [0]) << 56) | \
                        (((QWORD) ((const BYTE*) ptr) [1]) << 48) | \
                        (((QWORD) ((const BYTE*) ptr) [2]) << 40) | \
                        (((QWORD) ((const BYTE*) ptr) [3]) << 32) | \
                        (((QWORD) ((const BYTE*) ptr) [4]) << 24) | \
                        (((QWORD) ((const BYTE*) ptr) [5]) << 16) | \
                        (((QWORD) ((const BYTE*) ptr) [6]) << 8 ) | \
                        (((QWORD) ((const BYTE*) ptr) [7])      ) )

#else

#define GET32_U(ptr)  ( (((DWORD) ((const BYTE*) ptr) [0]) << 24) | \
                        (((DWORD) ((const BYTE*) ptr) [1]) << 16) | \
                        (((DWORD) ((const BYTE*) ptr) [2]) << 8 ) | \
                        (((DWORD) ((const BYTE*) ptr) [3])      ) )
#define GET32_UR(ptr) (*(const DWORD*)(ptr))


#define GET64_U(ptr)  ( (((QWORD) ((const BYTE*) ptr) [0]) << 56) | \
                        (((QWORD) ((const BYTE*) ptr) [1]) << 48) | \
                        (((QWORD) ((const BYTE*) ptr) [2]) << 40) | \
                        (((QWORD) ((const BYTE*) ptr) [3]) << 32) | \
                        (((QWORD) ((const BYTE*) ptr) [4]) << 24) | \
                        (((QWORD) ((const BYTE*) ptr) [5]) << 16) | \
                        (((QWORD) ((const BYTE*) ptr) [6]) << 8 ) | \
                        (((QWORD) ((const BYTE*) ptr) [7])      ) )
#define GET64_UR(ptr) (*(const QWORD*)(ptr))

#endif

bool cpu_simd ();

#if defined (_MSC_VER)
longlong atoll (const char* strval);
#endif

char* lltoa (longlong val, char* buf, int base);
char* ulltoa (ulonglong val, char* buf, int base);

#endif //__PLATFORM_H__

