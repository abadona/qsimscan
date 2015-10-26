
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

#ifndef __bitops_h__
#define __bitops_h__

#include <ostream>
#include <climits>
#include "platform.h"
#include "tracer.h"
#include "common_str.h"

#define BITS_PER_BYTE CHAR_BIT
#define BITS_PER_BYTE_SHIFT 3

#define BITS_PER_QWORD (BITS_PER_BYTE * sizeof (QWORD))
#define BASES_PER_QWORD (BASES_PER_BYTE * sizeof (QWORD))


// helper class template for use in bitw function
template <typename T>
class BITW
{
public:
    enum { bitw = sizeof (T) * BITS_PER_BYTE };
};

// bit width of an variable of given type
template <typename T>
inline BYTE bitw ()
{
    return BITW<T>::bitw;
}


template <int SZ>
class EXBITW
{
public:
    enum { exbitw = EXBITW < (SZ+1)/2 >::exbitw + 1 };
};

template <>
class EXBITW < 1 >
{
public:
    enum { exbitw = 0 };
};

// exponential bit width - how many address bits are 'eaten up' by an instance of given type
template <typename T>
inline BYTE exbitw ()  // exponent bit width
{
    return EXBITW < sizeof (T) >::exbitw;
}


// writes the bit-wise representation of passed value to the given stream
// uses C-style byte layout: writes "most significant" component first.
// This is opposite to the memory layout on little-endian (intel-based) CPUs
// this template is layout-agnostic, it does not change order depending on 'endiness' of architecture
// To get byte order that matches little-endian memory layout, pass rev as 'false'
template <typename TYPE>
inline void print_bits (TYPE& value, std::ostream& ostr, bool rev = true)
{
    BYTE* base = (BYTE*) &value;
    unsigned begb = rev ?  sizeof (TYPE) : 1;
    unsigned endb = rev ?  0 : sizeof (TYPE) + 1;
    unsigned step = rev ? -1 : 1 ;
    for (unsigned bpos = begb; bpos != endb; bpos += step)
    {
        if (bpos != begb) ostr << SPACE_STR;
        print_bits (base [bpos-1], ostr, rev);
    }
}


// specialization of print_bits for 1-byte values
template <>
inline void print_bits <BYTE> (BYTE& val, std::ostream& o, bool rv)
{
    // assume all architectures are "bitwise little-endian" (Most Significant Bit is on the left)
    for (unsigned pos = BITS_PER_BYTE; pos != 0; --pos)
        o << ((val >> (pos-1)) & 1);
}

// writes the byte-wise (hexadecimal) representation of passed value to the given stream
// uses C-style byte layout: writes "most significant" component first.
// This is opposite to the memory layout on little-endian (intel-based) CPUs
// this template is layout-agnostic, it does not change order depending on 'endiness' of architecture
// To get byte order that matches little-endian memory layout, pass rv as 'false'
template <typename TYPE>
inline void print_bytes (TYPE& value, std::ostream& ostr, bool rev = true)
{
    BYTE* base = (BYTE*) &value;
    unsigned begb = rev ?  sizeof (TYPE) : 1;
    unsigned endb = rev ?  0 : sizeof (TYPE) + 1;
    unsigned step = rev ? -1 : 1 ;
    for (unsigned bpos = begb; bpos != endb; bpos += step)
    {
        if (bpos) ostr << SPACE_STR;
        ostr << std::hex << (unsigned) base [bpos-1];
    }
}


template <typename TYPE>
inline TYPE invert_bits (TYPE val)
{
    TYPE retval = 0;
    unsigned pos = (sizeof (TYPE) << 3);
    while (pos)
    {
        retval <<= 1;
        retval |= val & 1;
        val >>= 1;
        pos --;
    }
    return retval;
}

// helper class for printing out 32-bit masks
struct MASK
{
    public:
        DWORD mask_;
        BYTE l_;
        MASK (DWORD m, BYTE l)
        :
        mask_ (m),
        l_ (l)
        {
        }
        MASK ()
        :
        l_ (0)
        {
        }
};

inline std::ostream& operator << (std::ostream& ostr, const MASK& mask)
{
    DWORD mm = mask.mask_ & ((~ (DWORD) 0) >> (sizeof (DWORD) * 8 - mask.l_));
    print_bits (mm, ostr);
    return ostr;
}

inline Logger& operator << (Logger& logger, const MASK& mask)
{
    if (logger.enabled ())
        logger.o_ << mask;
    return logger;
}



#endif // __bitops_h__
