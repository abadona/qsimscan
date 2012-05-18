
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

#include "p_kt_scan.h"
#include <rerror.h>
#include <sciminmax.h>
#include <transitive_closure.h>
#include <common_algo.h>

#include <assert.h>

//#define DEBUG_NEG_START

#ifdef DEBUG_NEG_START
#include <iostream>
#endif

// DEBUG puropses
// #define TESTING

#ifdef TESTING
#include <stdio.h>
#else
void skip__ (const char* fmt, ...)
{
}
#define printf skip__
#endif


#define check_consistency(x) if (!(x)) ers << "Internal:" << ERRINFO << Throw;

void BAND::merge (BAND* toMerge, DIAGONAL_ENTRY* diags, int band_idx, int toMerge_idx, int* hits)
{
    DIAGONAL_ENTRY *cd;
    int didx;

    leftmost_ = min_ (leftmost_, toMerge->leftmost_);
    rightmost_ = max_ (rightmost_, toMerge->rightmost_);
    min_off_ = min_ (min_off_, toMerge->min_off_);
    max_off_ = max_ (max_off_, toMerge->max_off_);

    for (didx = leftmost_; didx <= rightmost_; didx ++)
    {
        cd = diags + didx;
        if (cd->band_ == toMerge_idx)
            cd->band_ = band_idx;
    }
    if (toMerge->best_score_ > best_score_)
        best_score_ = toMerge->best_score_;
    if (toMerge->hit_idx_ != -1 && hit_idx_ == -1)
    {
        hit_idx_ = toMerge->hit_idx_;
        toMerge->hit_idx_ = -1;
        hits [hit_idx_] = band_idx;
    }
    else
        toMerge->skip_ = true;
}

void BAND::add (int diag_idx, int offset, int tuple_size, DIAGONAL_ENTRY* diags, int band_idx)
{
    leftmost_  = min_ (leftmost_, diag_idx);
    rightmost_ = max_( rightmost_, diag_idx);
    min_off_   = min_ (min_off_, offset);
    max_off_   = max_ (max_off_, offset + tuple_size);
    (diags + diag_idx)->band_ = band_idx;
}

bool BAND::overlaps (BAND* other, int widen, int extend)
{
    // non_ovl: e1 < b2 || e2 < b1
    if (rightmost_ + widen < other->leftmost_ || other->rightmost_ + widen < leftmost_)
        return false;
    if (max_off_ + extend < other->min_off_ || other->max_off_ + extend < min_off_)
        return false;
    return true;
}

int calc_tuple_index (char* tuple, int tuple_size, int alphabet_size = ALPHABET_SIZE)
{
    int result = 0;
    unsigned char cur;
    for (int idx = 0; idx < tuple_size; idx ++, tuple ++)
    {
        result *= alphabet_size;
        cur = *tuple;
        if (cur >= ALPHABET_SIZE) cur = ALPHABET_SIZE-2; // HACK: replace any unknown code with 'X'
        result += cur;
    }
    return result;
}

void decode_tuple (int tuple_idx, char* tuple, int tuple_size, int alphabet_size = ALPHABET_SIZE)
{
    for (int pos = tuple_size - 1; pos >= 0; pos--)
    {
        // good code optimizer should merge % and / into one assembly instruction. Check and replace by assembly code if neccessary
        tuple [pos] = tuple_idx % alphabet_size;
        tuple_idx /= alphabet_size;
    }
}

void calc_tuple_auto_scores (short* scores, int tuple_idx, int tuple_size, WMatrix& weight_matrix, int alphabet_size = ALPHABET_SIZE)
{
    short result = 0;
    int symbol;
    for (int pos = tuple_size - 1; pos >= 0; pos --)
    {
        // good code optimizer should merge % and / into one assembly instruction. Check and replace by assembly code if neccessary
        symbol = tuple_idx % alphabet_size;
        tuple_idx /= alphabet_size;
        result += weight_matrix.mx [symbol][symbol];
        scores [pos] = result;
    }
}


PKTSCAN::PKTSCAN (WMatrix* w, ResultReciever_pblast* res, int tuple_size, int max_queries, int max_tot_queries_len, int max_target_len, int max_band, int max_hit, int max_batch)
:
weight_matrix_ (w),
results_ (res),
tuple_size_ (tuple_size),
max_queries_ (max_queries),
max_tot_queries_len_ (max_tot_queries_len),
max_target_len_ (max_target_len),
max_band_ (max_band),
max_hit_ (max_hit),
max_batch_ (max_batch)
{
    // init allocatable ptrs
    diags_ = NULL;
    tuples_ = NULL;
    increments_ = NULL;
    queries_ = NULL;
    hits_ = NULL;
    bands_ = NULL;
    batches_ = NULL;
    entries_ = NULL;

    // check parameters for consistency
    if (tuple_size_ > MAX_TUPLE_SIZE) ers << "Tuple size is too long" << Throw;

    // initialize
    init_vars ();

    printf ("\nPKTSCAN created, tuple_size = %d, total %d tuples", tuple_size_, total_tuples_);
}

PKTSCAN::~PKTSCAN ()
{
    if (increments_) delete [] increments_;
    if (aligner_) delete aligner_;
    if (entries_) delete entries_;
    if (bands_) delete bands_;
    if (hits_) delete [] hits_;
    if (batches_) delete [] batches_;
    if (diags_) delete [] diags_;
    if (tuples_) delete [] tuples_;
    if (queries_) delete [] queries_;
}

void PKTSCAN::max_shift (int opt)
{
    //check_consistency (opt <= MAX_MAX_SHIFT && opt >= 0);
    max_shift_ = opt;
}
void PKTSCAN::init_vars ()
{
    queries_number_ = 0;
    accum_queries_len_ = 0;
    increments_ = NULL;
    entries_ = NULL;

    // allocate queries array
    queries_ = NULL;
    try {
        queries_ = new SEQUENCE_INFO [max_queries_];
    } catch (std::bad_alloc&) {}
    if (!queries_) ers << "unable to allocate queries array" << Throw;

    // count number of tuples (assume 24-letter alphabet)
    total_tuples_ = 1;
    for (int pos = 0; pos < tuple_size_; pos++) total_tuples_ *= ALPHABET_SIZE;

    // allocate tuples_
    tuples_ = NULL;
    try {
        tuples_ = new TUPLE_INFO [total_tuples_];
    } catch (std::bad_alloc&) {}
    if (!tuples_) ers << "unable to allocate tuples array" << Throw;
    memset (tuples_, 0, sizeof (TUPLE_INFO)*total_tuples_);

    // innit tuple auto-scores
    init_tuple_auto_scores ();

    // init diagonal info array
    diags_ = NULL;
    try {
        diags_ = new DIAGONAL_ENTRY [max_target_len_ + max_tot_queries_len_];
    } catch (std::bad_alloc&) {}
    if (!diags_) ers << "unable to allocate diagonal info array" << Throw;
    memset (diags_, 0xff, sizeof (DIAGONAL_ENTRY) * (max_target_len_ + max_tot_queries_len_));

    // init batches array
    batches_ = NULL;
    try {
        batches_ = new BATCH [max_batch_];
    } catch (std::bad_alloc&) {}
    if (!batches_) ers << "unable to allocate batches array" << Throw;

    // init bands array
    try {
        bands_ = new BAND [max_band_];
    } catch (std::bad_alloc&) {}
    if (!bands_) ers << "unable to allocate bands array" << Throw;

    // init hits array
    hits_ = NULL;
    try {
        hits_ = new int  [max_hit_];
    } catch (std::bad_alloc&) {}
    if (!hits_) ers <<  "unable to allocate hits array" << Throw;

    aligner_ = NULL;

    // assign default values for search parameters
    k_thresh_   = DEFAULT_K_THRESH;
    max_shift_  = DEFAULT_MAX_SHIFT;
    step_ = DEFAULT_STEP;
    widen_factor_ = DEFAULT_WIDEN_FACTOR;
    extend_factor_ = DEFAULT_EXTEND_FACTOR;
    distance_factor_ = DEFAULT_DIST_FACT;
}

void PKTSCAN::init_tuple_auto_scores ()
{
    TUPLE_INFO* sentinel = tuples_ + total_tuples_;
    // localize vars used in loop
    int tuple_size = tuple_size_;
    WMatrix& wm = *weight_matrix_;
    int tuple_idx;
    TUPLE_INFO* tuple;
    for (tuple = tuples_, tuple_idx = 0; tuple < sentinel; tuple ++, tuple_idx ++)
        calc_tuple_auto_scores (tuple->scores_, tuple_idx, tuple_size, wm);
}
/*
void PKTSCAN::set_tuple_weights (int* weights)
{
    TUPLE_INFO* tptr = tuples_;
    for (int tuple_idx = 0;
        tuple_idx < total_tuples_;
        tuple_idx ++, tptr ++, weights ++)
        tptr->weight_ = *weights;
}
*/

bool PKTSCAN::add_query (SEQ& qseq)
{
    if (queries_number_ > max_queries_) ers << "Too many queries, limit declared in parameters exceeded" << Throw;
    if (accum_queries_len_ + qseq.len > max_tot_queries_len_) ers << "Total queries length exceeded limit declared in parameters" << Throw;
    SEQUENCE_INFO& added_query = queries_ [queries_number_];
    added_query.ref_ = &qseq;
    added_query.query_id_ = qseq.uid;
    added_query.len_ = qseq.len;
    added_query.seq_ = qseq.seq;
    added_query.start_ = accum_queries_len_;
    accum_queries_len_ += qseq.len;
    queries_number_ ++;

    // fix possible invalid symbols in query
    char* last = added_query.seq_ + added_query.len_;
    for (char* symbol_ptr = added_query.seq_; symbol_ptr < last; symbol_ptr ++)
    {
        // coorrect invalid query symbols
        if (*symbol_ptr < 0 || *symbol_ptr >= ALPHABET_SIZE)
            *symbol_ptr = ALPHABET_SIZE-1;
    }
    return true;
}

void PKTSCAN::compute_lookup_space (double threshold)
{
    printf ("\ncomputing lookup space, thr = %f", threshold);

    printf ("\n%d queries, total %d aas", queries_number_, accum_queries_len_);

    int tot_hits = 0;
    total_diversified_tuples_ = 0;
    total_diversified_entries_ = 0;
    total_entries_ = 0;

    int skip = 1;
    for (int pos = tuple_size_ - 1; pos >= 0; pos --)
    {
        skips_ [pos] = skip;
        skip *= ALPHABET_SIZE;
    }

    printf ("\ncounting primary hits");
    diversity_threshold_ = threshold;
    // count primary (exact matching) tuples
    for (SEQUENCE_INFO* query = queries_; query < queries_ + queries_number_; query ++)
    {
        // for every position in query:
        char* last_tuple = query->seq_ + query->len_ - tuple_size_;
        for (char* tuple_ptr = query->seq_; tuple_ptr < last_tuple; tuple_ptr ++)
        {
            // increment tuples_[tuple_index].entries_no_
            tuples_ [calc_tuple_index (tuple_ptr, tuple_size_)].entries_no_ ++;
            tot_hits ++;
        }
    }
    printf ("\n%d primary hits recorded", tot_hits);

    // calculate the number of entries that has to be added to each tuple_info due to diversification
    increments_ = new int [total_tuples_];
    if (!increments_) ers << "unable to allocate tuple entries increment tables" << Throw;
    memset (increments_, 0, total_tuples_ * sizeof (int));

    printf ("\ncounting diversification increments");
    int self_score;
    int last_tuple_idx;
    TUPLE_INFO* cur_tuple;
    for (cur_tuple = tuples_, cur_tuple_idx_ = 0; cur_tuple_idx_ < total_tuples_; cur_tuple ++, cur_tuple_idx_++)
    {
        if (cur_tuple_idx_ % 100000 == 0) printf ("\r%d of %d: %d events", cur_tuple_idx_, total_tuples_, total_diversified_tuples_);

        if (cur_tuple->entries_no_ == 0)
            continue;
        partial_scores_ = cur_tuple->scores_;
        self_score = *partial_scores_;
        // do not diversify negative and zero self-scoring tuples
        if (self_score <= 0)
            continue;
        score_threshold_ = (int) (diversity_threshold_ * self_score);
        decode_tuple (cur_tuple_idx_, decoded_tuple_, tuple_size_);
        last_tuple_idx = count_diversity ();
        check_consistency (last_tuple_idx == total_tuples_);

    }
    printf ("\n%d diversifications, %d entries added\n", total_diversified_tuples_, total_diversified_entries_);

    total_entries_ = tot_hits + total_diversified_entries_;
}

void PKTSCAN::fill_lookup_table ()
{
    printf ("\nfilling lookup table");

    // allocate the arrays
    printf ("\nallocating space for entries");
    entries_ = new TUPLE_ENTRY [total_entries_];
    if (!entries_) ers <<  "unable to allocate the array for entries" << Throw;
    TUPLE_ENTRY* fill_level = entries_;
    TUPLE_ENTRY* max_fill = entries_ + total_entries_;

    TUPLE_INFO* cur_tuple;
    for (cur_tuple_idx_ = 0, cur_tuple = tuples_; cur_tuple_idx_ < total_tuples_; cur_tuple ++, cur_tuple_idx_ ++)
    {
        // allocate buffer for entries
        int entries_no = cur_tuple->entries_no_ + increments_ [cur_tuple_idx_];
        if (entries_no == 0)
            continue;

        cur_tuple->entries_ = fill_level;
        fill_level += entries_no;
        check_consistency (fill_level <= max_fill);
        // reset count to zero - it will be used for filling up
        // cur_tuple->pad_ = entries_no; // store alloc space in pad_ for debugging
        cur_tuple->entries_no_ = 0;
    }

    printf ("\nfilling primary (exact) hits");
    int tot_hits = 0;
    int tuple_size = tuple_size_;
    // fill in the primary (exact matching) tuples
    TUPLE_ENTRY* entry;
    for (SEQUENCE_INFO* query = queries_; query < queries_ + queries_number_; query ++)
    {
        // for every position in query:
        char* last_tuple_ptr = query->seq_ + query->len_ - tuple_size;
        for (char* tuple_ptr = query->seq_; tuple_ptr < last_tuple_ptr; tuple_ptr ++)
        {
            cur_tuple_idx_ = calc_tuple_index (tuple_ptr, tuple_size);
            cur_tuple = tuples_ + cur_tuple_idx_;
            // increment tuples_[tuple_index].entries_no_
            entry = cur_tuple->entries_  + cur_tuple->entries_no_;
            entry->query_idx_ = query - queries_;
            entry->offset_ = max_target_len_ + query->start_ + (tuple_ptr - query->seq_); // global!
            memcpy (entry->scores_, cur_tuple->scores_, sizeof (entry->scores_));
            // check_consistency (cur_tuple->entries_no_ < cur_tuple->pad_);
            cur_tuple->entries_no_ ++;
            tot_hits ++;
        }
    }
    printf ("\n%d primary hits recorded", tot_hits);

    // fill in the diversifications

    // reset increments - they will be used as current fill levels
    memset (increments_, 0, total_tuples_ * sizeof (int));

    printf ("\nfilling in diversifications");

    for (cur_tuple_idx_ = 0, cur_tuple = tuples_; cur_tuple_idx_ < total_tuples_; cur_tuple_idx_ ++, cur_tuple ++)
    {
        if (cur_tuple_idx_ % 100000 == 0) printf ("\r%d of %d", cur_tuple_idx_, total_tuples_);

        if (cur_tuple->entries_no_ == 0)
            continue;

        partial_scores_ = cur_tuple->scores_;
        int self_score = *partial_scores_;
        // do not diversify negative and zero self-scoring tuples
        if (self_score <= 0)
            continue;
        score_threshold_ = (int) (diversity_threshold_ * self_score);
        decode_tuple (cur_tuple_idx_, decoded_tuple_, tuple_size);
        int last_tuple_idx = fill_diversity ();
        check_consistency (last_tuple_idx == total_tuples_);
    }

    printf ("\nupdating entry counts");

    // update entry_no_ for every TUPLE_INFO to account for the increment
    int* incr_ptr = increments_;
    for (cur_tuple = tuples_; cur_tuple < tuples_ + total_tuples_; cur_tuple++, incr_ptr ++)
        cur_tuple->entries_no_ += *incr_ptr;

    // increments are not needed anymore - we may delete them (thow we may need to use them for debug :)
    delete [] increments_;
    increments_ = NULL;

    // it might make sense to sort the entry arrays here to improve locality (however time loss on sort may be ore then gain by locality) - need experimenting
}

int PKTSCAN::count_diversity (int position, int prefix_score, int matching_tuple_idx)
{
    TUPLE_INFO* source_tuple;
    // skip those letters at <position> that can not reach the score of (diversity_threshold_ * auto_score) - prefix_score
    int rest_score = prefix_score;
    if (position < tuple_size_ - 1)
        rest_score += partial_scores_ [position + 1];

    for (int symbol = 0; symbol < ALPHABET_SIZE; symbol ++)
    {
        int cur_symbol_score = weight_matrix_->mx [symbol][decoded_tuple_ [position]];
        int cur_best_score = cur_symbol_score + rest_score;
        if (cur_best_score < score_threshold_)
        {
            matching_tuple_idx += skips_ [position];
        }
        else if (position != tuple_size_ - 1)
        {
            matching_tuple_idx = count_diversity (position + 1, prefix_score + cur_symbol_score, matching_tuple_idx);
        }
        else
        {
            // avoid self-duplication
            if (matching_tuple_idx != cur_tuple_idx_)
            {
                // this tuple is similar enough to cur_tuple_idx_
                source_tuple = tuples_ + cur_tuple_idx_;
                increments_ [matching_tuple_idx] += source_tuple->entries_no_;
                total_diversified_entries_ += source_tuple->entries_no_;
                total_diversified_tuples_ ++;
            }
            matching_tuple_idx ++;
        }
    }
    return matching_tuple_idx;
}

int PKTSCAN::fill_diversity (int position, int prefix_score, int matching_tuple_idx)
{
    // skip those letters at <position> that can not reach the score of (diversity_threshold_ * auto_score) - prefix_score
    int pos;
    short *s;
    int accum;
    TUPLE_INFO* source_tuple;
    TUPLE_INFO* dest_tuple;
    TUPLE_ENTRY* source_entry;
    TUPLE_ENTRY* dest_entry;

    int rest_score = prefix_score;
    if (position < tuple_size_ - 1)
        rest_score += partial_scores_ [position + 1];
    for (int symbol = 0; symbol < ALPHABET_SIZE; symbol ++)
    {
        int cur_symbol_score = weight_matrix_->mx [symbol][decoded_tuple_ [position]];
        cur_match_scores_ [position] = cur_symbol_score;
        int cur_best_score = cur_symbol_score + rest_score;
        if (cur_best_score < score_threshold_)
        {
            matching_tuple_idx += skips_ [position];
        }
        else if (position != tuple_size_ - 1)
        {
            matching_tuple_idx = fill_diversity (position + 1, prefix_score + cur_symbol_score, matching_tuple_idx);
        }
        else
        {
            // avoid self-duplication
            if (matching_tuple_idx != cur_tuple_idx_)
            {
                // this tuple is similar enough to cur_tuple_idx_
                source_tuple = tuples_ + cur_tuple_idx_;
                dest_tuple = tuples_ + matching_tuple_idx;
                for (int entry_idx = 0; entry_idx < source_tuple->entries_no_; entry_idx ++)
                {
                    source_entry = source_tuple->entries_ + entry_idx;
                    dest_entry = dest_tuple->entries_ + dest_tuple->entries_no_ + increments_ [matching_tuple_idx] + entry_idx;
                    dest_entry->query_idx_ = source_entry->query_idx_;
                    dest_entry->offset_ = source_entry->offset_;
                    accum = 0;
                    s = dest_entry->scores_;
                    for (pos = tuple_size_ - 1; pos >= 0; pos --)
                    {
                        accum += cur_match_scores_ [pos];
                        s [pos] = accum;
                    }
                }
                increments_ [matching_tuple_idx] += source_tuple->entries_no_;
                // check_consistency (dest_tuple->entries_no_ + increments_ [matching_tuple_idx] <= dest_tuple->pad_)
            }
            matching_tuple_idx ++;
        }
    }
    return matching_tuple_idx;
}

void PKTSCAN::init_search ()
{
    printf ("\ninitializing search");
    int idx;

    // calc query sequence composition
    double composition [ALPHABET_SIZE];
    for (idx = 0; idx < ALPHABET_SIZE; idx ++) composition [idx] = 0.0;
    int max_query_len = 0;
    for (SEQUENCE_INFO* query_ptr = queries_; query_ptr < queries_ + queries_number_; query_ptr ++)
    {
        max_query_len = max_ (query_ptr->len_, max_query_len);
        char* last = query_ptr->seq_ + query_ptr->len_;
        for (char* symbol_ptr = query_ptr->seq_; symbol_ptr < last; symbol_ptr ++)
            composition [*symbol_ptr] += 1;
    }
    printf ("Total queries len = %d", accum_queries_len_);
    for (idx = 0; idx < ALPHABET_SIZE; idx ++) composition [idx] /= (double (accum_queries_len_) / ALPHABET_SIZE);

    // calc average self match and average mismatch scores
    ave_self_match_ = 0.0;
    min_self_match_ = 100000.0;
    for (idx = 0; idx < ALPHABET_SIZE; idx ++)
    {
        double self_match = composition [idx] * weight_matrix_->mx [idx][idx];
        ave_self_match_ += self_match;
        //printf ("self_match for %d = %f\n", idx, self_match);
        if (self_match > 0)
            min_self_match_ = min_ (min_self_match_, self_match);
    }
    ave_self_match_ /= ALPHABET_SIZE;
    check_consistency (min_self_match_ != 100000.0);

    //printf ("consistent\n");

    ave_mismatch_ = 0.0;
    for (int idx1 = 0; idx1 < ALPHABET_SIZE; idx1 ++)
        for (int idx2 = 0; idx2 < ALPHABET_SIZE; idx2 ++)
        {
            if (idx2 >= idx1) break;
            ave_mismatch_ += composition [idx1] * composition [idx2] *  weight_matrix_->mx [idx1][idx2];
        }
    // printf ("ave mism = %f\n", ave_mismatch_);
    ave_mismatch_ /= ((ALPHABET_SIZE * ALPHABET_SIZE) / 2 - ALPHABET_SIZE);
    check_consistency (ave_mismatch_ < 0);

    // printf ("\nAverage self-match score on query set is %f", ave_self_match_);
    // printf ("\nAverage mismatch score on query set is %f", ave_mismatch_);

    // calc diag_thresold_ from k_thresh_
    diag_threshold_ = ave_self_match_ * k_thresh_;

    // allocate the aligner knowing maximum query sequence size
    aligner_ = new ALIGN (weight_matrix_, max_target_len_, max_query_len*max_target_len_);

    target_ord_ = 0;

    extension_ = max_ (0, int ((diag_threshold_/ave_self_match_)*extend_factor_));

#ifdef DEBUG_NEG_START
    std::cerr << "init_search" << std::endl;
    std::cerr << "    tupe_size: " << tuple_size_ << std::endl;
    std::cerr << "    diversity_threshold: " << diversity_threshold_ << std::endl;
    std::cerr << "    k_thresh: " << k_thresh_ << std::endl;
    std::cerr << "    max_shift: " << max_shift_ << std::endl;
    std::cerr << "    ave_self_match: " << ave_self_match_ << std::endl;
    std::cerr << "    min_self_match: " << ave_mismatch_ << std::endl;
    std::cerr << "    ave_mismatch: " << ave_mismatch_ << std::endl;
    std::cerr << "    diag_threshold: " << diag_threshold_ << std::endl;
    std::cerr << "    distance_factor: " << distance_factor_ << std::endl;
    std::cerr << "    extend_factor: " << extend_factor_ << std::endl;
    std::cerr << "    widen_factor: " << widen_factor_ << std::endl;
    std::cerr << "    total_diversified_entries: " << total_diversified_entries_ << std::endl;
    std::cerr << "    total_diversified_tuples: " << total_diversified_tuples_ << std::endl;
    std::cerr << "    max_target_len: " << max_target_len_ << std::endl;
    std::cerr << "    max_query_len: " << max_query_len << std::endl;
    std::cerr << "    weight_matrix: gip " << weight_matrix_->gip << ", gep " << weight_matrix_->gep << std::endl;
    std::cerr << "                     minmx " << weight_matrix_->minmx << ", maxmx " << weight_matrix_->maxmx << std::endl;
    int x1, x2;
    for (x1 = 0; x1 < AANUM; x1 ++)
    {
        for (x2 = 0; x2 < AANUM; x2 ++)
            std::cerr << weight_matrix_->mx [x1][x2] << " ";
        std::cerr << std::endl;
    }
#endif
}


void PKTSCAN::search (SEQ& tseq)
{
    if (tseq.len > max_target_len_) ers << "target sequence length exceeds declared maximum length" << Throw;

    target_ = &tseq;

    // correect invalid target symbols
    char* last = tseq.seq + tseq.len;
    for (char* symbol_ptr = tseq.seq; symbol_ptr < last; symbol_ptr ++)
    {
        // coorrect invalid query symbols
        if (*symbol_ptr < 0 || *symbol_ptr >= ALPHABET_SIZE)
        {
#ifdef DEBUG_NEG_START
            std::cerr << "bad target symbol '" << *symbol_ptr << "' at " << symbol_ptr - tseq.seq << std::endl;
#endif
            *symbol_ptr = ALPHABET_SIZE-1;
        }
    }

    target_len_ = tseq.len;
    target_id_ = tseq.uid;
    target_seq_ = tseq.seq;

    diag_scanner ();

    target_ord_ ++;
}

void PKTSCAN::reset_y ()
{
    // reset diagonal info array
    // reset_diags (); - no need anymore - diag_scanner cleans up
    // unallocate tuple entries and reset tuple info
    for (TUPLE_INFO* tuple_ptr = tuples_; tuple_ptr < tuples_ + total_tuples_; tuple_ptr ++)
    {
        if (tuple_ptr->entries_)
        {
            tuple_ptr->entries_ = NULL;
            tuple_ptr->entries_no_ = 0;
        }
        check_consistency (tuple_ptr->entries_no_ == 0);
    }
    delete [] entries_;
    entries_ = NULL;
    total_entries_ = 0;

    accum_queries_len_ = 0;
    queries_number_ = 0;
}

void PKTSCAN::diag_scanner ()
{
    // alloc space on stack to avoid in-loop SP/BP manipulations
    int cur_tuple_idx;
    TUPLE_INFO* cur_tuple;
    TUPLE_ENTRY* cur_entry;
    TUPLE_ENTRY* sentinel_entry;
    int diag_idx;
    DIAGONAL_ENTRY* cur_diag;
    BAND* cur_band;
    BAND* adj_band;
    int test_band_idx;
    BAND* test_band;
    int query_idx;  // index of the current query
    int query_start;
    int offset;     // offset for the current entry
    int dist;       // distance from previous to current hit on the diagonal
    int cur_match_score;
    double score;
    int shift;
    double penalty;
    DIAGONAL_ENTRY* adj_diag;
    DIAGONAL_ENTRY* best_adj_diag;
    double adj_diag_score;
    double best_adj_diag_score;
    double new_score;
    //DIAGONAL_ENTRY* all_adj_diags [MAX_MAX_SHIFT*2];
    //int adj_diags_count;


    bool need_new_band;

    // localize some vars for faster (stack) access
    char* target_seq = target_seq_;
    int target_ord = target_ord_;
    int tuple_size = tuple_size_;
    TUPLE_INFO* tuples = tuples_;
    DIAGONAL_ENTRY* diags = diags_;
    BAND* bands = bands_;
    double ave_mismatch = ave_mismatch_;
    int max_shift = max_shift_;
    double gip = weight_matrix_->gip;
    double gep = weight_matrix_->gep;
    double diag_threshold = diag_threshold_;
    double distance_factor = distance_factor_;

    int max_target_len = max_target_len_;
    int max_diag_idx = max_target_len + accum_queries_len_;
    int max_band_idx = max_tot_queries_len_;
    int cur_band_idx;
    int bands_count = 0;

    //int didx;
    //DIAGONAL_ENTRY* cd;

    // reset hits
    hits_count_ = 0;

    // the flag indicating that pre-allocated resource reached the limit; causes skip of further processing
    bool over_limit = false;

    // for every position from 0 to (target_len_ - tuple_size_)
    int last_pos = target_len_ - tuple_size_;
    for (int target_pos = 0; (target_pos < last_pos) && !over_limit; target_pos += step_, target_seq += step_)
    {
        // make the tuple index
        cur_tuple_idx = calc_tuple_index (target_seq, tuple_size_);
        // find the tuple info
        cur_tuple = tuples + cur_tuple_idx;

        // for every entry of the tuple in the query:
        sentinel_entry = cur_tuple->entries_ + cur_tuple->entries_no_;
        for (cur_entry = cur_tuple-> entries_; cur_entry < sentinel_entry; cur_entry ++)
        {
            offset = cur_entry->offset_;
            query_idx = cur_entry->query_idx_;
            query_start = (queries_ + query_idx)->start_;

            // find the diagonal
            diag_idx = offset - target_pos;
            cur_diag = diags + diag_idx;
            cur_match_score = *cur_entry->scores_;
            score = cur_match_score;

            // find the band
            cur_band = NULL;

            // if diag belongs to band for current sequence pair
            cur_band_idx = cur_diag->band_;
            if (cur_band_idx != -1  && (cur_diag->target_idx_ == target_ord) && ((bands + cur_band_idx)->query_idx_ == query_idx))
            {
                dist = offset - cur_diag->offset_; // distance from prev hit on this diagonal

                if (dist < 0) // if overlaps with prev k-tuple match - add suffix score
                {
                    score = cur_diag->score_ + cur_entry->scores_ [-dist];
                    cur_band = bands + cur_band_idx;
                }
                else // if does not overlap with prev k-tuple match - penalize offset (by average mismatch score). Can not go lower then current match score.
                {
                    new_score = double (cur_match_score + cur_diag->score_) + ave_mismatch * dist * distance_factor;
                    if (new_score >= cur_match_score) // start new band
                        score = new_score, cur_band = bands + cur_band_idx;
                }
            }

            //look at previous matches on adjacent diagonals, pick the best score
            // for offs from 1 to max_shift
            best_adj_diag = NULL;
            best_adj_diag_score = 0;
            //adj_diags_count = 0;
            for (shift = 1, penalty = gip; shift <= max_shift; shift ++, penalty += gep)
            {
                // go down; transfer if more then current
                if (diag_idx + shift < max_diag_idx)
                {
                    adj_diag = cur_diag + shift;

                    test_band_idx = adj_diag->band_;
                    if (test_band_idx != -1)
                    {
                        test_band = bands_ + test_band_idx;

                        if ((test_band->query_idx_ == query_idx) && (adj_diag->target_idx_ == target_ord))
                        {
                            dist = offset - adj_diag->offset_;
                            if (dist <= -tuple_size)
                                adj_diag_score = adj_diag->score_ - penalty;
                            else if (dist < 0)
                                adj_diag_score = adj_diag->score_ + cur_entry->scores_ [-dist] - penalty;
                            else
                                adj_diag_score = adj_diag->score_ + cur_match_score + dist*ave_mismatch*distance_factor - penalty;

                            if (adj_diag_score > score && best_adj_diag_score < adj_diag_score)
                            {
                                best_adj_diag_score = adj_diag_score;
                                best_adj_diag = adj_diag;
                            }
                        }
                    }
                }

                // go left; transfer if more then current
                if (diag_idx - shift >= 0)
                {
                    adj_diag = cur_diag - shift;

                    test_band_idx = adj_diag->band_;
                    if (test_band_idx != -1)
                    {
                        test_band = bands_ + test_band_idx;

                        if ((test_band->query_idx_ == query_idx) && (adj_diag->target_idx_ == target_ord))
                        {
                            dist = offset - shift - adj_diag->offset_;
                            if (dist <= -tuple_size)
                                adj_diag_score = adj_diag->score_ - penalty;
                            else if (dist < 0)
                                adj_diag_score = adj_diag->score_ + cur_entry->scores_ [-dist] - penalty;
                            else
                                adj_diag_score = adj_diag->score_ + cur_match_score + dist*ave_mismatch*distance_factor - penalty;

                            if (adj_diag_score > score && best_adj_diag_score < adj_diag_score)
                            {
                                best_adj_diag_score = adj_diag_score;
                                best_adj_diag = adj_diag;
                            }
                        }
                    }
                }
            }

            if (best_adj_diag)
            {
                adj_band = bands_ + best_adj_diag->band_;
                if (cur_band && best_adj_diag->band_ != cur_band_idx)
                {
                    // merge bands if there is continuity on diag (see above)
                    // check_consistency (cur_band->best_score_ < diag_threshold || adj_band->best_score_ >= diag_threshold);
                    adj_band->merge (cur_band, diags, best_adj_diag->band_, cur_band_idx, hits_);
                }

                cur_band = adj_band;
                score = best_adj_diag_score;
            }

            cur_diag->target_idx_ = target_ord;
            cur_diag->score_ = score;
            cur_diag->offset_ = offset + tuple_size;

            if (!cur_band)
            {
                // allocate new band; init it
                if (bands_count == max_band_)
                {
                    std::cerr << "WARNING: Too many bands detected per target. Extra bands will be ignored. You may want to increase MAX_BAND." << ERRINFO << std::endl;
                    over_limit = true;
                    break;
                    // ers << "too many bands per target" << Throw;
                }
                else
                {
                    cur_diag->band_ = bands_count;
                    cur_band = bands + bands_count ++;

                    cur_band->query_idx_ = query_idx;
                    cur_band->rightmost_ = diag_idx;
                    cur_band->leftmost_ = diag_idx;
                    cur_band->min_off_ = offset;
                    cur_band->max_off_ = offset + tuple_size;
                    cur_band->best_score_ = 0;
                    cur_band->skip_ = false;
                    cur_band->hit_idx_ = -1;
                }
            }
            else
            {
                // update band with the new match
                check_consistency (!cur_band->skip_);
                cur_band->add (diag_idx, offset, tuple_size, diags, cur_band - bands);
            }

            if (cur_band->best_score_ < diag_threshold && score >= diag_threshold)
            {
                if (hits_count_ == max_hit_)
                {
                    // ers << "Too many hits per target" << Throw;
                    std::cerr << "WARNING: Too many hits detected per target. Extra hits will be ignored. You may want to increase MAX_HIT." << ERRINFO << std::endl;
                    over_limit = true;
                    break;
                }
                else
                {
                    cur_band->hit_idx_ = hits_count_;
                    hits_ [hits_count_ ++] = cur_band - bands_;
                }
            }
            cur_band->best_score_ = max_ (cur_band->best_score_, score);
        }
    }

    if (!hits_count_)
        return;

    int real_hits_count = 0;
    for (int hit_idx = 0; hit_idx < hits_count_; hit_idx ++)
    {
        cur_band = bands_ + hits_ [hit_idx];
        if (!cur_band->skip_)
            real_hits_count ++;
    }

    if (real_hits_count > 1)
    {
        // make transitive closure

        Closure closure (hits_count_);
        // add all pairs that could be merged
        for (unsigned i1 = 0; i1 < hits_count_; i1 ++)
        {
            BAND* band1 = bands + *(hits_ + i1);
            if (band1->skip_) continue;
            for (unsigned i2 = 0; i2 < i1; i2 ++)
            {
                BAND* band2 = bands + *(hits_ + i2);
                if (band2->skip_)
                    continue;
                if (compatible (band1, band2))
                    closure.add (i2, i1);
            }
        }
        closure.finalize ();
        Clusters clusters;
        closure.fillClusters (clusters, true);
        // now merge content of each cluster into first member, marking other members as 'skipped'
        for (Clusters::iterator clitr = clusters.begin (); clitr != clusters.end (); clitr ++)
        {
            if (bands_ [clitr->front ()].skip_)
                continue;
            merge_cluster (*clitr);
        }
    }

    // now walk the hits and flush good ones them;
    for (int hit_idx = 0; hit_idx < hits_count_; hit_idx ++)
    {
        cur_band = bands_ + hits_ [hit_idx];
        if (!cur_band->skip_)
            batch_assembler (cur_band);
    }
}


double p_score (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, WMatrix* wm, double* x_auto, double* y_auto)
{
    double pscore = 0;
    double xscore = 0;
    double yscore = 0;
    int prevxend = -1;
    int prevyend = -1;

    for (int idx = 0; idx < b_cnt; idx ++)
    {
        BATCH& batch = b_ptr [idx];

        // calc gaps
        if (prevxend != -1)
        {
            if (prevxend < batch.xpos)
                pscore -= (wm->gip + (batch.xpos - prevxend) * wm->gep);
            if (prevyend < batch.ypos)
                pscore -= (wm->gip + (batch.ypos - prevyend) * wm->gep);
        }
        // calc matching scores
        int xaa, yaa;
        for (int bpos = 0; bpos < batch.len; bpos ++)
        {
            xaa = xseq [bpos + batch.xpos];
            yaa = yseq [bpos + batch.ypos];

            xscore += wm->mx [xaa][xaa];
            yscore += wm->mx [yaa][yaa];
            pscore += wm->mx [yaa][xaa];
        }

        prevxend = batch.xpos + batch.len;
        prevyend = batch.ypos + batch.len;
    }
    if (x_auto) *x_auto = xscore;
    if (y_auto) *y_auto = yscore;
    return pscore;
}


void PKTSCAN::batch_assembler (BAND* band)
{
    // call aligner here - banded; call backtrace; wrap and report hit to results_
    int query_idx = band->query_idx_;
    SEQUENCE_INFO* query = queries_ + query_idx;
    SEQ* query_seq = query->ref_;

    int half_width = (band->rightmost_ - band->leftmost_ + 1) / 2;
    // int query_pos = band->min_off_ + (band->leftmost_ - (query->start_ + max_target_len_));
    int query_pos = band->min_off_ - (query->start_ + max_target_len_);
    // int target_pos = band->min_off_ - half_width;
    int target_pos = band->min_off_ - band->rightmost_ + half_width;
    int len  = band->max_off_ - band->min_off_;
    int width = max_ (1, half_width);
    width += widen_factor_;

    // extend start by average undetectible length
    // clip by seq start
    int down_ext = min_ (extension_, query_pos);
    down_ext = min_ (down_ext, target_pos);

    query_pos -= down_ext;
    target_pos -= down_ext;

    // extend to the end by tuple_size and clip len if needed
    len += 2*extension_ + tuple_size_;

    // clip by sequence boundary
    if (query_pos + len > query_seq->len) len = query_seq->len - query_pos;
    // if (target_pos + len > target_->len) len = target_->len - target_pos;

    // int score_a = aligner_->align (*query_seq, *target_);

    int score_a = aligner_->align_band (*query_seq, *target_, query_pos, target_pos, len, width);

    int batch_no = aligner_->backtrace (batches_, max_batch_, width);

    check_consistency (batch_no > 0);
    check_consistency (batches_->xpos >= 0);
    check_consistency (batches_->ypos >= 0);

    double query_auto;
    double target_auto;

    BATCH* lastb = batches_ + (batch_no-1);

    check_consistency (query->len_ >= lastb->xpos + lastb->len);
    //check_consistency (target_->len >= lastb->ypos + lastb->len);

    p_score (query_seq->seq, target_->seq, batches_, batch_no, weight_matrix_, &query_auto, &target_auto);

    results_->match_found (*query_seq, *target_, batches_, batch_no, score_a, query_auto, target_auto);

}

bool PKTSCAN::compatible (BAND* band1, BAND* band2)
{
    // checks if close enough:
    if (band1->query_idx_ != band2->query_idx_)
        return false;

    if (band1->overlaps (band2, widen_factor_, extension_))
        return true;
    // gap cost for distance between bands should be less then any band score
    // make sure ban2 is after band1 diagonal-wise
    if (band2->rightmost_ < band1->leftmost_)
        { BAND* t = band2; band2 = band1; band1 = t; }
    int hd = 0, vd = 0;
    if (band1->max_off_ < band2->min_off_)
    {
        hd = band2->min_off_ - band1->max_off_;
        vd = max_ (0, hd - (band1->rightmost_ - band2->leftmost_));
    }
    else if (band1->min_off_ > band2->max_off_)
    {
        hd = band1->min_off_ - band2->max_off_;
        vd = max_ (0, hd + band1->rightmost_ - band2->leftmost_);
    }
    else // overlap
    {

        if (band1->rightmost_ < band2->leftmost_)
            vd = band2->leftmost_ - band1->rightmost_; // actually 0-bonding is not needed if these are non overlapping bands
        else if (band2->rightmost_ < band1->leftmost_)
            vd = band1->leftmost_ - band2->rightmost_; // actually 0-bonding is not needed if these are non overlapping bands

        // vd = max_ (0, band1->rightmost_ - band2->leftmost_); // actually 0-bonding is not needed if these are non overlapping bands
    }
    double cost = 0;
    if (hd)
        cost += weight_matrix_->gip + weight_matrix_->gep * hd;
    if (vd)
        cost += weight_matrix_->gip + weight_matrix_->gep * vd;

    return cost <= min_ (band1->best_score_, band2->best_score_);
}

void PKTSCAN::merge_cluster (UIntVect& cluster)
{
    // merges content of each cluster into first member, marking other members as 'skipped'
    UIntVect::iterator itr = cluster.begin ();
    BAND* band = bands_ + *(hits_ + *itr);
    int& rightmost = band->rightmost_;
    int& leftmost = band->leftmost_;
    int& min_off = band->min_off_;
    int& max_off = band->max_off_;
    double& score = band->best_score_;
    itr ++;
    for (; itr != cluster.end (); itr ++)
    {
        band = bands_ + *(hits_ + *itr);;
        rightmost = max_ (band->rightmost_, rightmost);
        leftmost = min_ (band->leftmost_, leftmost);
        min_off = min_ (band->min_off_, min_off);
        max_off = max_ (band->max_off_, max_off);
        score += band->best_score_;
        band->skip_ = true;
    }
}

