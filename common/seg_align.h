
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

#ifndef __seg_align_h__
#define __seg_align_h__

#include <common_typedefs.h>
#include "align_result.h"
#include "weights.h"

// SegAlign merges local similarity segments between two sequences into best possible arrangement
// (with given gap cost cap and maximal segment overlap)
// It rearranges the passed in similarities by relocating and changing the batches and scores
// and removing the ones that were joined in.

class SegAlign
{
    WMatrix& wm_;
    float gip_;
    float gep_;
    float gcap_;
    int max_ovl_;

    struct BTrace
    {
        int idx_;
        float score_;
        int gap_;
        BTrace () : idx_ (-1), score_ (0), gap_ (0) {}
    };

    typedef std::vector <BTrace> BTraceVect;

    class CompareTraceIdxByScore
    {
        BTraceVect& trace_;
    public:
        CompareTraceIdxByScore (BTraceVect& trace) : trace_ (trace) {}
        bool operator () (int i1, int i2) const
        {
            float res = trace_ [i1].score_ - trace_ [i2].score_;
            if (res > 0) return true; // > to produce descending order
            if (res < 0) return false;
            return trace_ [i1].gap_ < trace_ [i2].gap_; // if scores are same, pick the one with smaller gap first
        }
    };
    class IsContinuation
    {
    public:
        bool operator () (const BTrace& tr) const
        {
            return (tr.idx_ != -1);
        }
    };
    IsContinuation isContinuation;

    void fill_trace (ARVect& sims, UIntVect& order, BTraceVect& trace);
    unsigned find_continuations (BTraceVect& trace, BoolVect& continuations);
    void make_score_order (BTraceVect& trace, BoolVect& continuations, unsigned resno, UIntVect& target);
    void backtrace_from (unsigned idx, const char* tseq, BTraceVect& trace, ARVect& sims, BoolVect& processed, BoolVect& to_remove);
public:
    SegAlign (WMatrix& wm, float gip, float gep, float gcap, int max_ovl);
    bool merge (ARVect& sims, const char* tseq);
};

#endif
