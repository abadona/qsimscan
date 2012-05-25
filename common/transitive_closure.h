
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

#ifndef __transitive_closure_h__
#define __transitive_closure_h__

#include <vector>
#include <algorithm>
#include "common_typedefs.h"
#include "resource.h"

typedef std::vector <UIntVect> Clusters;

class Closure
{
    enum STAGE { INIT, FILL, COLLECT };
    MemWrapper <unsigned> closure_;
    unsigned sz_;
    unsigned clno_;
    unsigned singlet_no_;
    STAGE stage_;
public:
    Closure ();
    Closure (unsigned sz);
    void init (unsigned sz);
    void add (unsigned idx1, unsigned idx2);
    void finalize ();
    void fillClusters (Clusters& dest, bool include_singletons = false);
    Clusters* collect (bool include_singletons = false);
    unsigned clusters_number () const
    {
        return clno_;
    }
    unsigned size () const
    {
        return sz_;
    }
    unsigned singletons_number () const
    {
        return singlet_no_;
    }
    unsigned non_trivial_number () const
    {
        return sz_ - singlet_no_;
    }
    const unsigned* result ();
};

#endif
