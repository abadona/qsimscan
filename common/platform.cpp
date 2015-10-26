
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

#include <stdlib.h>
#include <algorithm>
#include "platform.h"

//SIMD support detection
//currently works on Intel Pentium+ and AMD CPUs only

#if defined (__x86__)

# if defined (_MSC_VER)

bool cpu_simd ()
{
    int simd_support = 0;
    __asm
    {
        mov     eax,            1           // Put a "1" in eax to tell CPUID to get the feature bits
        cpuid                               // Perform CPUID (puts processor feature info into EDX)
        test    edx,            02000000h   // Test bit 25, for Streaming SIMD Extensions existence.
        jz      NotFound
        mov     [simd_support], 1
        NotFound:
    }
    return simd_support != 0;
}

#elif defined (__GNUC__)

bool cpu_simd ()
{
    int ax, bx, cx, dx;
    int op = 1;

    asm ("cpuid"
        : "=a" (ax),
        "=b" (bx),
        "=c" (cx),
        "=d" (dx)
        : "a" (op)
        );
    return (dx & 0x02000000) != 0;
}

#else // unknown compiler, Intel plathform

bool cpu_simd ()
{
    return false; // unknown inline assembler
}

#endif

#else // Non-intel plathform

bool cpu_simd ()
{
    return false;
}

#endif

#if defined(_MSC_VER)
longlong atoll (const char* strval)
{
    return _atoi64 (strval);
}
char* lltoa (longlong val, char* buf, int base)
{
    return _i64toa (val, buf, base);
}
char* ulltoa (ulonglong val, char* buf, int base)
{
    return _ui64toa (val, buf, base);
}
#else
#include <algorithm>
char* lltoa (longlong val, char* buf, int base)
{
    // characters are 0-9A-Z
    int sign = 1;
    if (val < 0) val = -val, sign = -1;
    int pos = 0;
    do
    {
        char c = (char) (val % base);
        if (c <= 9) c += '0';
        else c += 'A';
        buf [pos ++] = c;
        val = val / base;
    }
    while (val > 0);
    if (sign == -1) buf [pos++] = '-';
    buf [pos] = 0;
    std::reverse (buf, buf + pos);
    return buf;
}

char* ulltoa (ulonglong val, char* buf, int base)
{
    // characters are 0-9A-Z
    int pos = 0;
    do
    {
        char c = (char) (val % base);
        if (c <= 9) c += '0';
        else c += 'A';
        buf [pos ++] = c;
        val = val / base;
    }
    while (val > 0);
    buf [pos] = 0;
    std::reverse (buf, buf + pos);
    return buf;
}

#endif

ulonglong atoull (const char* strval)
{
    ulonglong rv = 0;
    // skip leading whitespace
    while (*strval == ' ' || *strval == '\t') strval ++;
    // convert string to int
    while (*strval != '\0')
    {
        if (*strval >= '0' && *strval <= '9')
            rv = rv * 10 + *strval++ - '0';
        strval ++;
    }
    return rv;
}

