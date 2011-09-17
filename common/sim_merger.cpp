
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

#include "sim_merger.h"
#include <deque>
#include <sciminmax.h>
#include <transitive_closure.h>
#include <common_algo.h>
#include <assert.h>

bool SimMerger :: merge (ARVect& sims, const char* tseq)
{
    bool toR = false;
    // rep / threads : select best from overlap groups
    if (merge_rep_ || merge_thr_)
        if (merge_repeats (sims, tseq))
            toR = true;
    // dom : conenct the connectables in the best way (dynamic progr - defined)
    if (merge_dom_)
        if (merge_domains (sims, tseq))
            toR = true;
    return toR;
}

bool SimMerger :: merge_repeats (ARVect& sims, const char* tseq)
{
    // transitive closure for any overlapping for more then max_rep_orp_ / intersecting
    unsigned tot = sims.size ();
    Closure closure (tot);
    for (unsigned idx1 = 0; idx1 != tot; idx1 ++)
    {
        for (unsigned idx2 = 0; idx2 != tot; idx2 ++)
        {
            if (idx1 == idx2) continue;
            if ((merge_rep_ && sims [idx1].overlaps (sims [idx2], max_rep_orp_)) ||
                (merge_thr_ && sims [idx1].intersects (sims [idx2])))
                closure.add (idx1, idx2);
        }
    }
    closure.finalize ();
    Clusters result;
    closure.fillClusters (result, false); // do not include singletones
    // if no non-singletons, were done
    if (!result.size ())
        return false;
    // for each cluster, find best similarity; put others into toRemove list
    unsigned to_remove_sz = closure.non_trivial_number () - result.size ();
    MemWrapper <unsigned> to_remove (to_remove_sz);
    unsigned to_remove_idx = 0;
    for (Clusters::iterator citr = result.begin (); citr != result.end (); citr ++)
    {
        unsigned best_idx = find_best (*citr, sims);
        for (UIntVect::iterator mitr = (*citr).begin (); mitr != (*citr).end (); mitr ++)
        {
            if ((*mitr) != best_idx)
                to_remove [to_remove_idx ++] = *mitr;
        }
    }

    // sort toRemove list; remove them
    std::sort ((unsigned*) to_remove, ((unsigned*) to_remove) + to_remove_idx);
    remove_positions (sims, (unsigned*) to_remove, ((unsigned*) to_remove) + to_remove_idx);
    return true;
}

unsigned SimMerger :: find_best (const UIntVect& indices, const ARVect& sims) const
{
    unsigned best_idx = 0;
    int best_score = INT_MIN;
    for (UIntVect::const_iterator itr = indices.begin (); itr != indices.end (); itr ++)
        if (best_score < sims [*itr].score_)
        {
            best_score = sims [*itr].score_;
            best_idx = *itr;
        }
    return best_idx;
}


struct BestPrev
{
    int idx_;
    int score_;
    int gap_;
    BestPrev () : idx_ (-1), score_ (0), gap_ (0) {}
};

class CompareBestPrevIdxByScore
{
    BestPrev* arr_;
public:
    CompareBestPrevIdxByScore (BestPrev* arr) : arr_ (arr) {}
    bool operator () (int i1, int i2) const
    {
        int res = arr_ [i1].score_ - arr_ [i2].score_;
        if (res > 0) return true; // > to produce descending order
        if (res < 0) return false;
        return arr_ [i1].gap_ < arr_ [i2].gap_; // if scores are same, pick the one with smaller gap first
    }
};


bool SimMerger :: merge_domains (ARVect& sims, const char* tseq)
{
    return seg_aligner_.merge (sims, tseq);
}

