
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

#ifndef P_KT_SCAN_H
#define P_KT_SCAN_H

#include "weights.h"
#include "biosequence.h"
#include "align.h"
#include "result_reciever_pblast.h"
#include "common_typedefs.h"

#define MAX_QUERIES_DEFAULT 5000
#define MAX_TOTAL_QUERIES_LEN_DEFAULT 1000000 // 200 aa/seq
#define MAX_TARGET_LEN_DEFAULT 10000  // are there longer proteins? (longer translations may appear if searching DNA)
#define TUPLE_SIZE_DEFAULT 5 // 7,962,624 tuples
#define ALPHABET_SIZE 24 // NCBI weight matrix alphabet
#define MAX_TUPLE_SIZE 5 // 7,962,624 tuples
#define DEFAULT_DIVERSITY 0.82
#define DEFAULT_TUPLE_WEIGHT 100
#define MAX_BAND_DEFAULT 2000000 // 2 mln bands == 64 Mb space.
#define MAX_BATCH_DEFAULT 10000
#define MAX_HIT_DEFAULT  400000 // maximum number of hits per target sequence
//#define MAX_MAX_SHIFT 100 // maximum diagonal deviation supported

#define DEFAULT_MAX_SHIFT 4     // gap size accounted for in intial scan
#define DEFAULT_K_THRESH 20.0   // in self-correlation averages on the query sequence
#define DEFAULT_STEP 1 // step over subject sequence (shift between adjuscent tuples for lookup)
#define DEFAULT_EXTEND_FACTOR 2.0 // number of tismes to extend the band by average undetectible gap size
#define DEFAULT_WIDEN_FACTOR 3  // number of diagonals to widen the band before calling banded alignment
#define DEFAULT_DIST_FACT 1.0 // multiplier for diagonal distance penalty

#pragma pack (push, 1)

struct TUPLE_ENTRY
{
    unsigned short query_idx_;      // query sequence index
    int offset_;                    // position of last hit in comparison matrix (along query axis) = global_query_pos+max_target_len_
    short scores_ [MAX_TUPLE_SIZE]; // full and partial weights (for suffixes). Assuming the weight does not go over 5461
                                    // for compactness, this could be array of char, but then matrices PAM300 and above, as well as cdi matrices will fail
    // char pad [32-MAX_TUPLE_SIZE*sizeof (short)-2*sizeof (int)];
    // no need to pad - fits into 16 bytes perfectly! (no more then 65535 queries at once)
};

struct TUPLE_INFO
{
    // int weight_;           // the statistical weght of this ktuple
    int entries_no_;       // size of entries_ array
    TUPLE_ENTRY* entries_; // array of TUPLE_ENTRIES

    short scores_ [MAX_TUPLE_SIZE];
    // int pad_ [32 - (sizeof (short) * MAX_TUPLE_SIZE + sizeof (int) * 2 + sizeof (TUPLE_ENTRY*))];

    // replace init by memset - use 0 weight as 100.
    // TUPLE_INFO (): weight_ (DEFAULT_TUPLE_WEIGHT), entries_no_ (0), entries_ (NULL) {}
    // do not delete - the entries memory is controlled externally
    // ~TUPLE_INFO () { if (entries_) delete [] entries_; }

};


struct DIAGONAL_ENTRY
{
    int band_;          // index of band in bands_ array
    int offset_;        // position of last hit in comparison matrix (along query axis) = global_query_pos+max_target_len_
    double score_;       // accumulated score on current band
    int target_idx_;    // target number this diagonal used with last time

    void init ()
    {
        band_ = -1;
        target_idx_ = -1;
        offset_ = 0;
        score_ = 0.0;
    }
};

struct BAND
{
    // int start_diag_; // the diagonal number where the band has started

    // flags
    bool skip_ : 1;     // this band has been joined to another; do not report it separate
    int  hit_idx_ : 31;  // if this band was recorded as a hit, the hit number - to keep track while merging

    int query_idx_;     // the query index for this band

    int rightmost_;     // the righttmost diagonal of a hit in this band
    int leftmost_;      // the leftmost diagonal of a hit in this band
    int min_off_;       // minimal offset (along query axis) of a hit in a band in similarity matrix
    int max_off_;       // maximal offset (along query axis) of a hit in a band in similarity matrix

    double best_score_;  // the maximal score achieved by this band

    void merge  (BAND* toMerge, DIAGONAL_ENTRY* diags, int band_idx, int toMerge_idx, int* hits);
    void add    (int diag_idx, int offset, int tuple_size, DIAGONAL_ENTRY* diags, int band_idx);
    bool overlaps (BAND* other, int widen, int extend);
};

struct SEQUENCE_INFO
{
    SEQ* ref_;
    longlong query_id_;    // sequence id
    int start_;           // start in the diag_info array (less max_x_len)
    int len_;             // sequence length
    char* seq_;           // actual sequence

    SEQUENCE_INFO (): query_id_ (-1L), start_ (0), len_ (0), seq_ (NULL), ref_ (NULL) {}
};

#pragma pack (pop)

// The prototype for the protein blast
class PKTSCAN
{
    WMatrix* weight_matrix_;
    int tuple_size_;

    int max_queries_;
    int max_tot_queries_len_;
    int queries_number_;
    int accum_queries_len_;
    SEQUENCE_INFO* queries_;

    int max_target_len_;
    DIAGONAL_ENTRY* diags_;

    int total_tuples_;
    TUPLE_INFO* tuples_;

    int total_entries_;
    TUPLE_ENTRY* entries_;

    double diversity_threshold_;
    int* increments_;
    char decoded_tuple_ [MAX_TUPLE_SIZE];
    short* partial_scores_;
    short cur_match_scores_ [MAX_TUPLE_SIZE]; // this is NOT a suffix scores array - just a matching scores array. To be turned into suffix scores at the moment of entry recording during diversification
    int cur_tuple_idx_;
    int score_threshold_;
    int skips_ [MAX_TUPLE_SIZE];
    int total_diversified_entries_;
    int total_diversified_tuples_;

    longlong target_id_;
    int target_len_;
    char* target_seq_;

    double k_thresh_;
    int max_shift_;
    unsigned step_;

    double ave_self_match_;
    double min_self_match_;
    double ave_mismatch_;
    double diag_threshold_;
    double distance_factor_;

    double extend_factor_;
    unsigned int  widen_factor_;
    int extension_;

    int target_ord_; // the ordinal number of currently processed target in this job

    SEQ* target_;

    BAND* bands_; // storage for bands
    int max_band_;
    int* hits_;   // indexes of bands_ with scores over the threshold
    int max_hit_;
    int hits_count_;

    // batch_assembler data
    BATCH* batches_;
    int max_batch_;
    ALIGN* aligner_;
    ResultReciever_pblast* results_;

    void init_vars ();
    void init_tuple_auto_scores ();
    // void reset_diags (); - no need anymore - diag_scanner cleans up

    int count_diversity (int position = 0, int prefix_score = 0, int matching_tuple_idx = 0);
    int fill_diversity  (int position = 0, int prefix_score = 0, int matching_tuple_idx = 0);

    void diag_scanner ();
    void batch_assembler (BAND* band);
    bool compatible (BAND* band1, BAND* band2);
    void merge_cluster (UIntVect& cluster);


public:
    PKTSCAN (WMatrix* w, ResultReciever_pblast* res,
            int tuple_size = TUPLE_SIZE_DEFAULT,
            int max_queries = MAX_QUERIES_DEFAULT,
            int max_total_queries_len = MAX_TOTAL_QUERIES_LEN_DEFAULT,
            int max_target_len = MAX_TARGET_LEN_DEFAULT,
            int max_bands = MAX_BAND_DEFAULT,
            int max_hits = MAX_HIT_DEFAULT,
            int max_batch = MAX_BATCH_DEFAULT);
    virtual ~PKTSCAN ();

    // void set_tuple_weights (int* tuple_weights = NULL);
    bool add_query     (SEQ& qseq);
    void compute_lookup_space (double threshold = DEFAULT_DIVERSITY);
    void fill_lookup_table ();
    void init_search ();
    void search    (SEQ& xseq);
    void reset_y ();           // prepares for new cycle of search

    void k_thresh (double opt) { k_thresh_ = opt; }
    void max_shift (int opt);
    void step (int opt) { step_ = (unsigned) opt; }
    void extend_factor (double opt) {extend_factor_ = opt;}
    void widen_factor (int opt) {widen_factor_ = (unsigned) opt;}
    void dist_fact (double opt) {distance_factor_ = opt;}

};

#endif
