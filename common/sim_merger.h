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

#ifndef __sim_merger_h__
#define __sim_merger_h__

#include <common_typedefs.h>
#include <sciminmax.h>
#include "sim_merger_base.h"
#include "seg_align.h"

class SimMerger : public SimMergerBase
{
    enum
    {
        MAX_REP_ORP_DEF_ = 40,
        MAX_DOM_OVL_DEF_ = 5,
    };

    bool merge_rep_;
    bool merge_dom_;
    bool merge_thr_;
    unsigned max_rep_orp_;

    SegAlign seg_aligner_;


    bool merge_repeats (ARVect& sims, const char* tseq);
    bool merge_domains (ARVect& sims, const char* tseq);

    int calc_gap (const AlignResult& sim1, const AlignResult& sim2) const;
    unsigned find_best (const UIntVect& indices, const ARVect& sims) const;

public:
    SimMerger (bool merge_thr,
               bool merge_rep,
               bool merge_dom,
               WMatrix& wm,
               float gip,
               float gep,
               float gcap,
               int max_dom_ovl = MAX_DOM_OVL_DEF_,
               int max_rep_orp = MAX_REP_ORP_DEF_)
    :
    merge_rep_ (merge_rep),
    merge_dom_ (merge_dom),
    merge_thr_ (merge_thr),
    max_rep_orp_ (max_rep_orp),
    seg_aligner_ (wm, gip, gep, gcap, max_dom_ovl)
    {
    }
    bool merge (ARVect& sims, const char* tseq);
};

#endif
