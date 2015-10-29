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

#ifndef __SCIMINMAX_H__
#define __SCIMINMAX_H__

#ifndef _MSC_VER

	#include <algorithm>
	#define max_(a,b) std::max (a, b)
	#define min_(a,b) std::min (a, b)

#else

	#ifndef max_
		template <typename T>
		inline const T&
		max_ (const T& a, const T& b)
		{
			if (a < b)
				return b;
			else
				return a;
		}
	#endif

	#ifndef min_
		template <typename T>
		inline const T&
		min_ (const T& a, const T& b)
		{
			if (a < b)
				return a;
			else
				return b;
		}
	#endif

#endif

#if 0

    #ifndef max_
    #define max_(a,b) (((a)>(b))?(a):(b))
    #endif

    #ifndef min_
    #define min_(a,b) (((b)>(a))?(a):(b))
    #endif

#endif

#endif
