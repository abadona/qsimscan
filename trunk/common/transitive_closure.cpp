
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
// For any questions please contact SciDM team by email at scidmteam@yahoo.com
//////////////////////////////////////////////////////////////////////////////

#pragma warning (disable: 4786)
#include "transitive_closure.h"
#include "rerror.h"
#include <string.h>

Closure::Closure (unsigned sz)
:
stage_ (INIT)
{
    init (sz);
}

Closure::Closure ()
:
stage_ (INIT)
{
    init (0);
}

void Closure::init (unsigned sz)
{
    sz_ = sz;
    clno_ = sz;
    singlet_no_ = sz;
    stage_ = FILL;
    if (sz_)
        closure_ = new unsigned [sz];
    else
        closure_ = NULL;
    for (unsigned i = 0; i < sz_; i ++)
        closure_[i] = i;
    stage_ = FILL;
}
void Closure::add (unsigned idx1, unsigned idx2)
{
    // check for call consistency
    if (stage_ != FILL)
        ERR("Internal - adding links to finalized closure");
    if (idx1 >= sz_ || idx2 >= sz_)
        ERR("Internal - out-of bounds clusure access");

    // walk the link chains starting from idx1 down to the end
    while (idx1 != closure_ [idx1])
        idx1 = closure_ [idx1];

    // walk the link chains starting from idx2 down to the end
    while (idx2 != closure_ [idx2])
        idx2 = closure_ [idx2];

    // make the link if walk lead to different cells:
    // cell at higher index points to cell at lower one; lower remains pointing to itself
    if (idx1 < idx2)
    {
        closure_ [idx2] = idx1;
        clno_ --; // decrement cluster number if we actually made a link (merged two clusters in one
    }
    else if (idx2 < idx1)
    {
        closure_ [idx1] = idx2;
        clno_ --; // decrement cluster number if we actually made a link (merged two clusters in one
    }
}
Clusters* Closure::collect (bool include_singletons)
{
    // makes the vector of vectors containing values clustered together
    // the caller is responsible for deallocating of the result
    if (stage_ != COLLECT)
        finalize (); // finalize if necessary
    // allocate the result (vector of vectors). Wrap it for automatic resource management
    ObjWrapper <Clusters> result = new Clusters (); // the proper space allocated in fillClusters
    // fill resut data
    fillClusters (*result, include_singletons);
    // release wrapper and return result
    return result.release ();
}

void Closure::finalize ()
{
    if (stage_ != FILL)
        ERR ("Internal (Closure): finalize attempted not at filling stage");
    // set the value of the cells belonging to each cluster
    // to the index of lowest cell of those cluster
    for (unsigned idx = 0; idx < sz_; idx ++)
        if (closure_ [idx] != closure_ [closure_ [idx]])
            closure_ [idx] = closure_ [closure_ [idx]];
    // now all cluster members set to same value (cluster_id == lowest element index)
    stage_ = COLLECT;
}

void Closure::fillClusters (Clusters& dest, bool include_singletons)
{
    // count clusters; reserve appropriate space in results sub-vectors; make direct-lookup index
    // make the index from cluster_id (== lowest element idx) into cluster number.
    // Uses direct addressing to avoid the use of stl map (slow sue to many allocations)

    if (stage_ != COLLECT)
        ERR ("Internal (Closure): attempting to collect clusters without finalization");

    // clear dest in case it has some data
    dest.clear ();

    // do not bother doing something if closure is empty
    if (!sz_)
        return;

    // first, count how many times each position appears as cluster id
    MemWrapper <unsigned> cluster_index (sz_);
    std::fill ((unsigned*) cluster_index, ((unsigned*) cluster_index) + sz_, 0);
    unsigned idx;
    for (idx = 0; idx < sz_; idx ++)
        cluster_index [closure_ [idx]] ++;

    // reserve space
    int result_size = clno_;
    // if singletons are not wanted, reduce result_space by number of singletons
    if (!include_singletons)
    {
        unsigned* ptr = cluster_index;
        unsigned* sent = cluster_index + sz_;
        while (ptr < sent)
            if (*(ptr++) == 1)
                result_size --;
    }
    // fill dest with empty vectors
    dest.resize (result_size);

    // now walk over index again, allocating dest sub-vectors
    // and replacing sizes in cluster_index by cluster numbers
    // This makes actual direct-address index,
    // so that cluster_number = cluster_index [idx_of_lowest_element]
    unsigned cluster_count = 0;
    unsigned min_capacity = include_singletons ? 0:1; // it is not portable to use include_singletons directly
    for (idx = 0; idx < sz_; idx ++)
    {
        if (cluster_index [idx] > min_capacity)
        {
            // reduce singlet count
            if (cluster_index [idx] > 1)
                singlet_no_ -= cluster_index [idx];
            // reserve space
            dest [cluster_count].reserve (cluster_index [idx]);
            // make index entry
            cluster_index [idx] = cluster_count ++;
        }
        else
        {
            // mark index entry with centinel value
            cluster_index [idx] = sz_;
        }
    }
    // now make final pass over clusters and fill dest's sub-vectors
    unsigned cluster_idx;
    for (idx = 0; idx < sz_; idx ++)
    {
        cluster_idx = cluster_index [closure_ [idx]];
        if (cluster_idx != sz_) // idx == sentinel for singletons in case they are unwanted
            dest [cluster_idx].push_back (idx);
    }
}

const unsigned* Closure::result ()
{
    if (stage_ != COLLECT)
        ERR ("Internal: closure result request before finalization");
    return (unsigned*) closure_;
}
