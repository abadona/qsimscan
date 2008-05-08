#ifndef __blast_results_batch_h__
#define __blast_results_batch_h__

#include "sequtil.h"
#include "biosequence.h"
#include "result_reciever_blast_batch.h"
#include "align_result_storage.h"

#define MAX_BATCHES 400

class Search_helper_frag;

class BlastResultsBatch : public ResultReciever_blast_batch, public AlignResultStorage
{
    int total_found_;
    int passed_repeats_;

    int  rep_percent_;
    int  rep_len_;
    int* rep_buf_;

    bool triangle_only_;

public:
    BlastResultsBatch (int keep_per_query, int rep_len = 0, int rep_perc = 0, bool triangle_only = false);
    virtual ~BlastResultsBatch () {}
    virtual bool match_found  (NN_SEQ& xseq, NN_SEQ& yseq, BATCH* batches, int batch_no, int matches);
    int totalFound () const {return total_found_;}
    int passedRepeats () const {return passed_repeats_;}
};


#endif // __blast_results_frag_batch_h__
