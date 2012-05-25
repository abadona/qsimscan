
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
// For any questions please contact SciDM team by email at team@scidm.org
//////////////////////////////////////////////////////////////////////////////

#ifndef __common_algo_h__
#define __common_algo_h__

#include <vector>

template < typename Element, typename OrderedPosItr >
void remove_positions (std::vector < Element > & target, OrderedPosItr beg, OrderedPosItr end)
{
    // remove the elements of a vector at the indices from beg to end
    // std::remove_if is not used since it operates on values, not indexes, and proedicate's operator () cannot be non-const (and thus cannot count positions) (because remove_if copies it)
    int orig_pos = 0, new_pos = 0, tot = target.size ();
    while (orig_pos != tot)
    {
        if ((beg != end) && (orig_pos == *beg))
            // if matches item from removal indexes, skip
            beg ++;
        else
        {
            // otherwise, copy (if positions differ)
            if (new_pos != orig_pos)
                target [new_pos] = target [orig_pos];
            new_pos ++;
        }
        orig_pos ++;
    }
    // shrink container if needed
    if (new_pos != orig_pos)
        target.erase (target.begin () + new_pos, target.end ());
}

template < typename Element >
void remove_mask (std::vector < Element > & target, BoolVect& mask)
{
    // remove the elements of a vector for whose position mask [position] is true
    // std::remove_if is not used since it operates on values, not indexes, and proedicate's operator () cannot be non-const (and thus cannot count positions) (because remove_if copies it)
    int orig_pos = 0, new_pos = 0, tot = target.size ();
    while (orig_pos != tot)
    {
        if (!mask [orig_pos])
        {
            if (new_pos != orig_pos)
                target [new_pos] = target [orig_pos];
            new_pos ++;
        }
        orig_pos ++;
    }
    // shrink container if needed
    if (new_pos != orig_pos)
        target.erase (target.begin () + new_pos, target.end ());
}

template < typename Integral >
void range (std::vector < Integral > & target, unsigned number, Integral beg = (Integral) 0,  Integral step = (Integral) 1)
{
    target.resize (number);
    Integral v = beg;
    for (unsigned p = 0; p < number; p ++)
    {
        target [p] = v;
        v += step;
    }
}

#endif
