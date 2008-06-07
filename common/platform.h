
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

#ifndef __PLATFORM_H__
#define __PLATFORM_H__

#ifdef _MSC_VER
typedef __int64 longlong;
typedef unsigned __int64 ulonglong;
#else
typedef long long int longlong;
typedef long long unsigned int ulonglong;
#endif

#if defined (__x86__) || defined (__powerpc__)
#define SCIDM_LITTLE_ENDIAN
#elif defined (__mips__) || defined (__alpha__) || defined (__rx000__)
#define SCIDM_BIG_ENDIAN
#else
#error Byte order unknown
#endif

// the following works for MSVC on 32-bit intel and gcc on 32 and 64 bit (LP64 speciication) intel.
// Could fail on other platforms / compilers
typedef unsigned char BYTE;
typedef char SBYTE;
typedef unsigned short WORD;
typedef short SWORD;
#ifndef _MSC_VER
typedef unsigned int DWORD;
#else
typedef unsigned long DWORD;
#endif
typedef int SDWORD;
typedef ulonglong QWORD;
typedef longlong SQWORD;

#ifdef SCIDM_LITTLE_ENDIAN
#define GET32_U(ptr) (*(int*)(ptr))
#define GET32_UR(ptr) ((ptr[0]<<24)|(ptr[1]<<16)|(ptr[2]<<8)|ptr[3])
#else
#define GET32_U(ptr) ((ptr[0]<<24)|(ptr[1]<<16)|(ptr[2]<<8)|ptr[3])
#define GET32_UR(ptr) (*(int*)(ptr))
#endif

bool cpu_simd ();

#if defined (_MSC_VER)
longlong atoll (const char* strval);
#endif

char* lltoa (longlong val, char* buf, int base);
char* ulltoa (ulonglong val, char* buf, int base);

#endif //__PLATFORM_H__

