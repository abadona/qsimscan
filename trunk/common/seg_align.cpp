
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

#include "seg_align.h"
#include "weights.h"
#include <common_typedefs.h>
#include <common_algo.h>
#include <deque>
#include <assert.h>

#include <iostream>
#include <iomanip>

SegAlign::SegAlign (WMatrix& wm, float gip, float gep, float gcap, int max_ovl)
:
wm_ (wm),
gip_ (gip),
gep_ (gep),
gcap_ (gcap),
max_ovl_ (max_ovl)
{
}

class CompareSimIdxByPos
{
    const ARVect& sims_;
    int max_ovl_;
public:
    CompareSimIdxByPos (const ARVect& sims, int max_ovl) : sims_ (sims), max_ovl_ (max_ovl) {}
    bool operator () (int idx1, int idx2)
    {
        const BATCH& b1 = sims_ [idx1].batches_ [0];
        const BATCH& b2 = sims_ [idx2].batches_ [0];
        int vv = b1.xpos - b2.xpos;
        if (vv < 0) return true;
        if (vv > 0) return false;
        return b1.ypos < b2.ypos;
    }
};

static void make_pos_ordered_sim_index (ARVect& sims, int max_ovl, UIntVect& target)
{
    range (target, sims.size ());
    std::sort (target.begin (), target.end (), CompareSimIdxByPos (sims, max_ovl));
}

void SegAlign :: fill_trace (ARVect& sims, UIntVect& order, TraceVect& trace)
{
    // walk all similarities preceeding current by X axis
        // if continuation is possible
            // calculate and remember score; find best score / cutoff position for overlaps

    // For each position, right-to-left:
    for (UIntVect::iterator i1 = order.begin (); i1 != order.end (); i1 ++)
    {
        // select best of compatibles for each ending not later then current.beg + max_dom_ovl_
        for (UIntVect::iterator i2 = order.begin (); i2 != order.end (); i2 ++)
        {
            if (i1 == i2)
                continue;
            bool overrun;
            bool comp = sims [*i1].can_continue (sims [*i2], max_ovl_, overrun);
            if (overrun) break;
            if (!comp) continue;
            // compatible. calculate and remember score
            float cur_score = trace [*i2].score_ + sims [*i2].al_score_;
            float& prev_score = trace [*i1].score_;
            if (cur_score > prev_score)
            {
                prev_score = cur_score;
                trace [*i1].idx_ = *i2;
                trace [*i1].gap_ = trace [*i2].gap_ + sims [*i2].gap (sims [*i1]);
            }
            else if (cur_score == prev_score)
            {
                int cur_gap = trace [*i2].gap_ + sims [*i2].gap (sims [*i1]);;
                int& prev_gap = trace [*i1].gap_;
                if (cur_gap < prev_gap)
                {
                    // no need to assign score - it is same
                    trace [*i1].idx_ = *i2;
                    prev_gap = cur_gap;
                }
            }
        }
    }
}

unsigned SegAlign::find_continuations (TraceVect& trace, BoolVect& continuations)
{
    continuations.resize (trace.size ());
    std::fill (continuations.begin (), continuations.end (), false);
    for (int ii = 0; ii < trace.size (); ii ++)
        if (trace [ii].idx_ != -1)
            continuations [trace [ii].idx_] = true;

    return std::count (continuations.begin (), continuations.end (), false);
}

void SegAlign::make_score_order (TraceVect& trace, BoolVect& continuations, unsigned resno, UIntVect& target)
{
    unsigned p = 0;
    target.resize (resno);
    for (int ii = 0; ii < trace.size (); ii ++)
        if (!continuations [ii])
            target [p ++] = ii;
    std::sort (target.begin (), target.end (), CompareTraceIdxByScore (trace));
}

void SegAlign::backtrace_from (unsigned idx, TraceVect& trace, ARVect& sims, BoolVect& processed, BoolVect& to_remove)
{
    // backtrace:
    // while there is preceeding sim:
        // accumulate batches, removing / adjusting overlapping ones
        // remember the sim in to_remove list (except head, which does not have preceeding sim)
        // switch to preceeding sim
    // replace batches and adjust batch_no for head sim (last iterated)
    // adjust scores

    int iniidx = idx;

    // optimization for singlet - no need to merge
    // if (trace [idx].idx_ == -1 || processed [trace [idx].idx_])
    if (trace [iniidx].idx_ == -1 || processed [trace [iniidx].idx_])
        return;

    // merge all batches from merged sims;
    std::deque <BATCH> batches;
    int upd_sim_idx;
    int score = 0;
    float al_score = 0;
    double bitscore = 0;
    double evalue = 1e30;
    float chi2 = -1e30;
    int q_auto = 0;
    int t_auto = 0;
    while (idx != -1)
    {
        int next = trace [idx].idx_;
        // adjust tail of the added batch list so it does not overlap the head of allready added one.
        BATCH* nbp = sims [idx].batches_ + sims [idx].batch_no_ - 1;
        if (batches.size ())
        {
            BATCH& fb = batches.front ();
            while ((nbp->xpos >= fb.xpos) ||
                    (nbp->ypos >= fb.ypos))
                nbp --; // no check for getting below batch array start - relying upon prior compatibility check
            if (nbp < sims [idx].batches_)
                ERR (ERR_Internal);
            // assert (nbp >= sims [idx].batches_); // just in case prior compatibility is not working ...
            int xd = max_ (0, (int) (nbp->xpos + nbp->len) - (int) fb.xpos);
            int yd = max_ (0, (int) (nbp->ypos + nbp->len) - (int) fb.ypos);
            int adj = max_ (xd, yd);
            nbp->len -= adj;
        }
        // add adjusted batches to (front of) accumulator
        batches.insert (batches.begin (), (BATCH*) sims [idx].batches_, nbp + 1);
        al_score += sims [idx].score_; // - GAP cost ?
        score += sims [idx].score_;
        bitscore += sims [idx].bitscore_;
        evalue = min_ (evalue, sims [idx].evalue_);
        chi2 = max_ (chi2, sims [idx].chi2_);
        q_auto += sims [idx].q_auto_score_;
        t_auto += sims [idx].t_auto_score_;
        // save present sim into toRemove array (if it is not the last one)
        if (next != -1)
            to_remove [idx] = true;
        upd_sim_idx = idx;
        processed [idx] = true;
        idx = next;
    }
    // put new batches into sim being updated (the one at cur
    AlignResult& upd_sim = sims [upd_sim_idx];
    try {
        upd_sim.batches_ = new BATCH [batches.size ()];
    } catch (std::bad_alloc&) {upd_sim.batches_ = NULL; }
    if (!upd_sim.batches_) ERR(NOEMEM);
    
    upd_sim.batch_no_ = batches.size ();
    std::copy (batches.begin (), batches.end (), (BATCH*) upd_sim.batches_);
    // adjust the sw score (NOTE: gaps introduced at zero cost here! Also score is NOT adjusted for removed overlaps!)
    upd_sim.score_ = score;
    upd_sim.al_score_ = al_score;
    upd_sim.bitscore_ = bitscore;
    upd_sim.evalue_ = evalue;
    upd_sim.q_auto_score_ = q_auto;
    upd_sim.t_auto_score_ = t_auto;

}

// remember the sim in to_remove list (except head, which is being modified)
// dynamic programming
// find best set of connected paths over similarity segments

bool SegAlign::merge (ARVect& sims)
{
    if (sims.size () <= 1)
        return false; // nothing to do

    // Make list of sim indices ordered by xpos
    UIntVect order;

# if 0
    range (order, sims.size ());
    std::cerr << std::endl << "sims size : " << sims.size () << ", order size : " << order.size () << std::endl;
    for (int ii = 0; ii < sims.size (); ii ++)
    {
        std::cerr << std::setw (3) << ii << " : " << std::setw (3) << order [ii] << " : " << std::setw (3) << sims [ii].batch_no_ << " : ";
        for (int bi = 0; bi < sims [ii].batch_no_; bi ++)
        {
            BATCH& b = sims [ii].batches_ [bi];
            std::cerr << "(" << b.xpos << "," << b.ypos << "," << b.len << ")";
        }
        std::cerr << std::endl;
    }
    std::sort (order.begin (), order.end (), CompareSimIdxByPos (sims, max_ovl_));
#endif
    make_pos_ordered_sim_index (sims, max_ovl_, order);

    // Allocate array of prev_best control structures  (one per sim)
    TraceVect trace (sims.size ());

    // fill in the trace
    fill_trace (sims, order, trace);

    // mark indexes that have continuations; the ones that do not are merged groups tails
    BoolVect cont;
    unsigned resno = find_continuations (trace, cont);
    if (resno == sims.size ())
        return false; // nothing to do

    // make list of tail sim indices in the order descending scores
    UIntVect score_order;
    make_score_order (trace, cont, resno, score_order);

    // processed are seen sims; direct addressed to avoid log-time searches
    BoolVect processed (sims.size ()); std::fill (processed.begin (), processed.end (), false);
    // to_remove are not needed anymore; direct addressed to avoid log-time searches
    BoolVect to_remove (sims.size ()); std::fill (to_remove.begin (), to_remove.end (), false);

    // for each similarity without continuation (in score descending order):
    for (UIntVect::iterator titr = score_order.begin (); titr != score_order.end (); titr ++)
        backtrace_from (*titr, trace, sims, processed, to_remove);

    // remove sims at remove positions
    remove_mask (sims, to_remove);
    return true;

}
