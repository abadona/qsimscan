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

#ifndef __search_helper__
#define __search_helper__

#include "biosequence.h"
#include "align_result_storage.h"
#include "weights.h"
#include "fasta.h"

#include <limits>
#include <ostream>



class Parameters;

class Search_helper_files
{
protected:
    NN_SEQ n_xseq_;
    AA_SEQ p_xseq_;
    NA_SEQ na_xseq_;
    unsigned reserved_seqspace_;
    unsigned skipped_;
    unsigned target_no_;
    unsigned beg_;
    unsigned end_;

    bool iter_done_;
    bool iter_init_;

    bool rev_;
    FastaFile query_;
    FastaFile target_;
    std::string target_file_name_;

public:
    Search_helper_files ();

    // k-tuple distribution load from file
    int* load_ktuple_distribution (const char* objname, double sim_level, int& mismatch_score);
    // k-tuple distribution create
    int* create_ktuple_distribution (int ksize, double sim_level, int& mismatch_score);

    // load queries
    unsigned load_queries_nn (const char* fname, std::vector <NN_SEQ>& dest, unsigned min_len, unsigned max_len, unsigned* total = NULL, unsigned beg = 0, unsigned end = std::numeric_limits <unsigned>::max ());
    unsigned load_queries_aa (const char* fname, std::vector <AA_SEQ>& dest, unsigned min_len, unsigned max_len, unsigned* total = NULL, unsigned beg = 0, unsigned end = std::numeric_limits <unsigned>::max ());
    unsigned load_queries_na (const char* fname, std::vector <NA_SEQ>& dest, unsigned min_len, unsigned max_len, unsigned* total = NULL, unsigned beg = 0, unsigned end = std::numeric_limits <unsigned>::max ());

    // search; use SCIDM::INVALUD_ID for global search across all objects in db
    void init_nn_search (const char* fname, unsigned beg = 0, unsigned end = std::numeric_limits<unsigned>::max ());
    NN_SEQ* next_nn_seq (unsigned min_len, unsigned max_len);

    void init_aa_search (const char* fname, unsigned beg = 0, unsigned end = std::numeric_limits<unsigned>::max ());
    AA_SEQ* next_aa_seq (unsigned min_len, unsigned max_len);

    void init_na_search (const char* fname, unsigned beg = 0, unsigned end = std::numeric_limits<unsigned>::max ());
    NA_SEQ* next_na_seq (unsigned min_len, unsigned max_len);

    ulonglong current () const; // position, not record number
    ulonglong total () const; // file size, not number of records
    unsigned skipped () const;

    // results output
    bool output_results_nn (AlignResultStorage&, unsigned res_per_query, unsigned alignments_no, std::vector<NN_SEQ>& f_qry, std::vector<NN_SEQ>& r_qry, WMatrix* w, std::ostream& o, unsigned max_x_len);
    bool output_results_aa (AlignResultStorage&, unsigned res_per_query, unsigned alignments_no, std::vector<AA_SEQ>& qry, WMatrix* w, std::ostream& o, unsigned max_x_len);
    bool output_results_an (AlignResultStorage&, unsigned res_per_query, unsigned alignments_no, std::vector<AA_SEQ>& qry, WMatrix* w, std::ostream& o, unsigned max_x_len);
    bool output_results_na (AlignResultStorage&, unsigned res_per_query, unsigned alignments_no, std::vector<NA_SEQ>& f_qry, std::vector<NA_SEQ>& r_qry, WMatrix* w, std::ostream& o, unsigned max_x_len);

    bool output_results_m8_aa (AlignResultStorage& as, unsigned res_per_query, unsigned alignments_no, std::vector<AA_SEQ>& qry, WMatrix* w, std::ostream& o, const char* hdr = NULL);
    bool output_results_m8_nn (AlignResultStorage& as, unsigned res_per_query, unsigned alignments_no, std::vector<NN_SEQ>& f_qry, std::vector<NN_SEQ>& r_qry, WMatrix* w, std::ostream& o, const char* hdr = NULL);

    bool output_results_tab_aa (AlignResultStorage& as, unsigned res_per_query, unsigned alignments_no, std::vector<AA_SEQ>& qry, WMatrix* w, std::ostream& o, bool hdr = false);
    bool output_results_tab_nn (AlignResultStorage& as, unsigned res_per_query, unsigned alignments_no, std::vector<NN_SEQ>& f_qry, std::vector<NN_SEQ>& r_qry, WMatrix* w, std::ostream& o, bool hdr = false);
};




#endif	// __search_helper__
