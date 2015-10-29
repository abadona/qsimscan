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

#ifndef __nsimscan_h__
#define __nsimscan_h__

#include <fstream>
#include <platform.h>
#include <resource.h>
#include <process_thread.h>
#include <search_helper_files.h>
#include <weights.h>
#include "nsimscan_params.h"
#include <sim_merger_base.h>


#define MAX_SEARCH_NUM 64

class BlastResultsBatch;
class KT_SEARCH;


class Nsimscan : public Process
{
    enum SUBPHASE
    {
        NOTHING = -1,
        LOADING_QUERIES = 0,
        LOADING_TUPLE_WEIGHTS = 1,
        PREPARING_SEARCHER = 2,
        ADDING_QUERIES = 3,
        COUNTING_TUPLES = 4,
        FILLING_LOOKUP = 5,
        INIT_SEARCH = 6,
        SEARCHING = 7,
        WRITING_RESULTS = 8
    };

    Nsimscan_params* p_;

    Search_helper_files sh_;

    std::ofstream out_file_;

    unsigned query_set_length_;
    ulonglong search_set_length_;
    unsigned query_bps_;
    ulonglong searched_bps_;

    int mismatch_score_;
    MemWrapper <int> k_distr_;

    ObjWrapper <WMatrix> wm_;
    ObjWrapper <SimMergerBase> sim_merger_;
    ObjWrapper <BlastResultsBatch> results_;
    ObjWrapper <KT_SEARCH> searcher_;

    std::vector <NN_SEQ> f_yseqs_;
    std::vector <NN_SEQ> r_yseqs_;
    NN_SEQ* xseq_;

    unsigned queries_searched_;
    unsigned targets_searched_;
    unsigned skipped_;
    unsigned resno_;
    unsigned max_qry_len_;

    bool load_queries ();
    bool load_tuple_weights ();
    bool prepare_searcher ();
    bool add_queries ();
    bool count_tuples ();
    bool fill_lookup ();
    bool init_search ();
    bool search_next ();
    bool write_results ();

public:
    Nsimscan ();
    ~Nsimscan ();

    // handles
    virtual bool init_handler (Process_params* params);
    virtual bool next_handler ();
    virtual bool close_handler ();
    virtual bool report_results_handler ();

    // naming convention
    virtual const char* process_name ();
    virtual const char* description ();
    virtual const char* version ();
    virtual const char* processing_item_name ();
    virtual const char* result_item_name ();

    // process description
    virtual const char* phase_name (int phase_no);
    virtual int phases_total ();
    virtual const char* subphase_name (int phase_no, int subphase_no);
    virtual int subphases_total (int phase_no);

    virtual void report_state (std::ostream& o);
    virtual void report_self (std::ostream& o);
    virtual void report_final (std::ostream& o);
};


#endif // __nsimscan_h__
