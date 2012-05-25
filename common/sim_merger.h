
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

    SegAlign seg_aligner_;

    bool merge_thr_;
    bool merge_rep_;
    bool merge_dom_;
    unsigned max_rep_orp_;

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
