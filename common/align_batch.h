
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
