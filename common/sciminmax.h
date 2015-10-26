
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
