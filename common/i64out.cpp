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

#include "i64out.h"

#if defined (_MSC_VER)
std::ostream& operator << (std::ostream& o, longlong ii)
{
    char buf [64];
    lltoa (ii, buf, 10);
    o << buf;
    return o;
}
std::ostream& operator << (std::ostream& o, ulonglong ii)
{
    char buf [64];
    ulltoa (ii, buf, 10);
    o << buf;
    return o;
}
#endif


