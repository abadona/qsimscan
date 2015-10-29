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

#ifndef __blast_results_batch_h__
#define __blast_results_batch_h__

#include "sequtil.h"
#include "biosequence.h"
#include "result_reciever_blast_batch.h"
#include "merging_result_storage.h"

class Search_helper_frag;

class BlastResultsBatch : public ResultReciever_blast_batch, public MergingResultStorage
{
    static SimMergerBase null_merger_;

    int total_found_;
    int passed_repeats_;

    int  rep_percent_;
    int  rep_len_;
    int* rep_buf_;

    bool triangle_only_;

public:
    BlastResultsBatch (int keep_per_query, int rep_len = 0, int rep_perc = 0, bool triangle_only = false, SimMergerBase& merger = null_merger_);
    virtual ~BlastResultsBatch () {}
    virtual bool match_found  (NN_SEQ& xseq, NN_SEQ& yseq, BATCH* batches, int batch_no, int matches);
    int totalFound () const {return total_found_;}
    int passedRepeats () const {return passed_repeats_;}
};


#endif // __blast_results_frag_batch_h__
