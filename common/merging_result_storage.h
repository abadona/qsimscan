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

#ifndef __merging_result_storage_h__
#define __merging_result_storage_h__

#include <vector>
#include <map>
#include "align_result_storage.h"
#include "sim_merger_base.h"


class MergingResultStorage : public AlignResultStorage
{
    SimMergerBase& merger_;

    typedef std::vector < AlignResult > ARVect;
    typedef std::map < longlong, ARVect > QryResults;

    unsigned res_per_target_;
    QryResults accum_fwd_;
    QryResults accum_rev_;
    unsigned accum_no_;
    longlong cur_sid_;

    longlong tot_merged_;

    void clear_accum ();
public:
    MergingResultStorage (SimMergerBase& merger, unsigned capacity, unsigned res_per_tagert = 0);
    bool add_result (longlong qid, longlong sid, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, BATCH* batches, const char* binsubj = NULL, QWORD subjid = 0, DWORD subjlen = 0);
    void flush (const char* tseq);
    bool reset ();
    longlong tot_merged () const { return tot_merged_; }
};

#endif
