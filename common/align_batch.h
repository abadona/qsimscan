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

#ifndef __align_batch_h__
#define __align_batch_h__


#pragma pack(push, 1)
struct BATCH
{
    unsigned xpos;
    unsigned ypos;
    unsigned len;
    BATCH () : xpos ((unsigned) 0), ypos ((unsigned) 0), len ((unsigned) 0) {}
    BATCH (unsigned x, unsigned y, unsigned l) : xpos (x), ypos (y), len (l) {}
    bool intersects (const BATCH& other) const
    {
        // same diagonal and beg1 < end2 and beg2 < end1
        return ((int) ypos - (int) xpos == int (other.ypos) - int (other.xpos)) && (xpos < other.xpos + other.len) && (other.xpos < xpos + len);
    }
};
#pragma pack (pop)

#endif // __align_batch_h__
