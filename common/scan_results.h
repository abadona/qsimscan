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

#ifndef __scan_results_h__
#define __scan_results_h__
#pragma warning (disable:4786)
#include <map>
#include <platform.h>
#include <pqueue.h>
#include "result_reciever.h"
#include "biosequence.h"
#include "align_result_storage.h"
#include "weights.h"

enum SEQ_SCAN_TYPE
{
    SST_NN,
    SST_NN_REV,
    SST_AA,
    SST_AA_REV,
    SST_AN,
    SST_NA,
    SST_NA_REV
};


struct BATCH;
class ALIGN;
class ScanResults : public ResultReciever, public AlignResultStorage
{
public:
    ScanResults (ALIGN* align, int keep_per_query, SEQ_SCAN_TYPE scantype, WEIGHTS<int, 24>* wm, double score_thresh = 0, int rep_len = 0, double rep_perc = 0);
    virtual ~ScanResults ();

    virtual void match_found (SEQ* query_seq,  SEQ* search_seq, int score);

    int totalFound    () {return _total_found;   }
    int passedRepeats () {return _passed_repeats;}

private:

    void process_match_nn (NN_SEQ* query_seq, NN_SEQ* search_seq, int score);
    void process_match_aa (AA_SEQ* query_seq, AA_SEQ* search_seq, int score);
    void process_match_an (AA_SEQ* query_seq, NA_SEQ* search_seq, int score);
    void process_match_na (NA_SEQ* query_seq, AA_SEQ* search_seq, int score);
    void process_match_na_rev (AA_SEQ* search_seq, NA_SEQ* query_seq, int score);

    SEQ_SCAN_TYPE _scan_type;

    // filters
    ALIGN* _aligner;
    #define MAX_BATCHES 200
    BATCH _batches [MAX_BATCHES];

    WEIGHTS<int, 24>* _wm;

    int  _score_thresh;
    int  _rep_percent;
    int  _rep_len;
    int* _rep_buf;

    longlong _cur_id;

    int _total_found;
    int _passed_repeats;
};



#endif
