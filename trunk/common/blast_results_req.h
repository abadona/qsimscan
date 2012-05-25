
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

#ifndef __blast_results_req_h__
#define __blast_results_req_h__

#include <platform.h>
#include "biosequence.h"

struct req
{
    longlong uid;
    char rev;
    req ()
    {
        uid = -1;
        rev = 0;
    }
    req (const NN_SEQ& seq)
    {
        uid = seq.uid;
        rev = seq.rev;
    }
    bool operator < (const req& other) const
    {
        if (uid < other.uid) return true;
        else if (uid > other.uid) return false;
        else if (rev < other.rev) return true;
        else return false;
    }
    bool operator == (const req& other) const
    {
        return ((uid == other.uid) && (rev == other.rev));
    }
    bool operator == (const NN_SEQ& seq) const
    {
        return ((uid == seq.uid) && (rev == seq.rev));
    }
    void operator = (NN_SEQ& seq)
    {
        uid = seq.uid;
        rev = seq.rev;
    }
};

#endif // __blast_results_req_h__
