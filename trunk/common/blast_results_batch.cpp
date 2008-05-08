#pragma warning (disable: 4786)
#include <platform.h>
#include "blast_results_batch.h"
#include "align.h"
#include "filters.h"


BlastResultsBatch :: BlastResultsBatch (int keep_per_query, int rep_len, int rep_perc, bool triangle_only)
:
AlignResultStorage (keep_per_query),
rep_percent_ (rep_perc),
rep_len_ (rep_len),
total_found_ (0),
passed_repeats_ (0),
triangle_only_ (triangle_only)
{
    if (rep_len_)
    {
        rep_buf_ = new int [rep_len_];
        if (!rep_buf_) FATAL(NOEMEM);
    }
}

bool BlastResultsBatch :: match_found  (NN_SEQ& xseq, NN_SEQ& yseq, BATCH* batches, int batch_no, int matches)
{
    // for precluster - type searches:
    // We use the fact that fragment ids are fed in order, as they were obtained from set, where they are kept in order.
    // While searching, the Ids will come in strictly increasing order, which gives the possibility to cut out
    // those similarities that belong to orphan (last) triangles in striped search.
    if (triangle_only_)
        if (xseq.uid >= yseq.uid)
            return false;

    total_found_ ++;

    // By stupid but established convention, X IS QUERY while processing the results
    // so we'll swap batch here
    int tot_blen = swap_batches (batches, batch_no);

    // recalculate the score
    float nscore = nu_score_b (yseq.seq, xseq.seq, batches, batch_no, &matches);
    float chi2score = nu_all_chi2_b (yseq.seq, xseq.seq, batches, batch_no);

    // filter results by repeats
    if (rep_len_)
    {
        int max_idx = nu_offs_score_b (yseq.seq, xseq.seq, batches, batch_no, rep_len_, rep_buf_);
        if ((rep_buf_ [max_idx] * 100) / matches > rep_percent_)
            return false;
    }
    passed_repeats_ ++;

    // save the result
    return add_result (yseq.uid, xseq.uid, yseq.rev?true:false, (int) nscore, matches, chi2score, 0, 0, tot_blen, tot_blen, batch_no, batches);
}
