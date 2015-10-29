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
