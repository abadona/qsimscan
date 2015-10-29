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
