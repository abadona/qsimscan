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

#ifdef _MSC_VER
#pragma warning (disable: 4786)
#endif
#include "psimscan.h"
#include <fileutils.h>
#include <print_batches.h>
#include <biosequence.h>
#include <resource.h>
#include <weights.h>
#include <pblast_results.h>
#include <p_kt_scan.h>
#include <i64out.h>
#include <iomanip>
#include <sim_merger.h>

const char* VERSION = "0.9.6 (November 2018)";

Process* process_factory ()
{
    return new Psimscan ();
}

const char* process_name ()
{
    return "PSIMSCAN";
}

Psimscan::Psimscan ()
:
wm_ (NULL),
results_ (NULL),
searcher_ (NULL),
query_set_length_ (0),
search_set_length_ (0),
query_aas_ (0),
searched_aas_ (0),
queries_searched_ (0),
targets_searched_ (0),
skipped_ (0),
resno_ (0),
max_qry_len_ (0)
{
}

Psimscan::~Psimscan ()
{
    close_handler ();
}

bool Psimscan::init_handler (Process_params* params)
{
    p_ = (Psimscan_params*) params;
    // make sure the output is ready to receive data
    checkCreateOutputFile (p_->output_name (), p_->overwrite (), p_->append ());

    // check for existence of query set
    if (!file_exists (p_->query_name ()))
        ers << "ERROR: File " << p_->query_name () << " not found" << Throw;

    // check for existence of search set (if specified)
    if (!file_exists (p_->search_name ()))
        ers << "ERROR: File " << p_->search_name () << " not found" << Throw;

    total (0);
    current (0);

    query_set_length_ = 0;
    search_set_length_ = 0;
    query_aas_ = 0;
    searched_aas_ = 0;
    queries_searched_ = 0;
    targets_searched_ = 0;
    skipped_ = 0;
    resno_ = 0;
    max_qry_len_ = 0;

    subphase (LOADING_QUERIES);
    return true;
}

bool Psimscan::next_handler ()
{
    switch (subphase ())
    {
    case LOADING_QUERIES:        return load_queries ();
    case LOADING_TUPLE_WEIGHTS:  return load_tuple_weights ();
    case PREPARING_SEARCHER:     return prepare_searcher ();
    case ADDING_QUERIES:         return add_queries ();
    case COUNTING_TUPLES:        return count_tuples ();
    case FILLING_LOOKUP:         return fill_lookup ();
    case INIT_SEARCH:            return init_search ();
    case SEARCHING:              return search_next ();
    case WRITING_RESULTS:        return write_results ();
    default: ERR("Phasing error");
    }
    return false;
}

bool Psimscan::load_queries ()
{
    // load queries
    queries_searched_ = sh_.load_queries_aa (p_->query_name (), yseqs_, p_->min_seq_len (), p_->max_seq_len (), &query_set_length_, (unsigned) p_->q_beg (), (unsigned) p_->q_end ());

    // count number of bases in queries
    query_aas_ = 0;
    max_qry_len_ = 0;
    for (int iii = 0; iii < yseqs_.size (); iii ++)
    {
        int clen = yseqs_ [iii].len;
        if (max_qry_len_ < clen)
            max_qry_len_ = clen;
        query_aas_ += clen;

    }

    // report skipped queries
    if (p_->verbose () && queries_searched_ != query_set_length_)
        std::clog << query_set_length_ - queries_searched_ << " sequences were skipped due to length constraints" << std::flush;

    total (1);
    current (0);
    subphase (LOADING_TUPLE_WEIGHTS);
    return true;
}
bool Psimscan::load_tuple_weights ()
{
    subphase (PREPARING_SEARCHER);
    return true;
}

bool Psimscan::prepare_searcher ()
{
    // create weight matrix and result storage

    wm_ = new WEIGHTS<int, 24> (p_->matrix_name (), p_->matrix_name (), p_->gip (), p_->gep ());
    sim_merger_ = new SimMerger (true, p_->merge_repeats (), p_->merge_domains (), *wm_, p_->gip (), p_->gep (), p_->gap_cap ()*p_->gip (), p_->max_dom_ovl (), p_->max_rep_orp ());
    results_ = new PblastResults (p_->min_score(), p_->min_len (), wm_, p_->eval_eval(), p_->res_per_query (), p_->res_per_target (), *sim_merger_);

    // create the searcher
    searcher_ = new PKTSCAN (wm_, results_, p_->k_size (), p_->max_queries_number (), query_aas_, p_->max_seq_len ());

    // set search parameters
    searcher_->k_thresh(p_->k_thresh ());
    searcher_->max_shift (p_->max_shift ());
    searcher_->step (p_->step ());
    searcher_->extend_factor (p_->extend_band ());
    searcher_->widen_factor (p_->widen_band ());
    searcher_->dist_fact (p_->dist_fact ());

    subphase (ADDING_QUERIES);
    return true;
}


bool Psimscan::add_queries ()
{
    // fill y memory
    for (int yno = 0; yno < yseqs_.size (); yno ++)
    {
        // load forward
        if (!searcher_->add_query (yseqs_ [yno]))
            ERR("Error loading Y sequence");
    }
    subphase (COUNTING_TUPLES);
    return true;
}

bool Psimscan::count_tuples ()
{
    searcher_->compute_lookup_space (p_->approx ());
    subphase (FILLING_LOOKUP);
    return true;
}


bool Psimscan::fill_lookup ()
{
    searcher_->fill_lookup_table ();
    subphase (INIT_SEARCH);
    return true;
}

bool Psimscan::init_search ()
{
    // init search
    sh_.init_aa_search (p_->search_name (), (unsigned) p_->t_beg (), (unsigned) p_->t_end ());
    searcher_->init_search ();

    if (p_->t_end () == -1)
        total (sh_.total ());
    else
        total (p_->t_end () - p_->t_beg ());
    current (0);
    skipped_ = 0;
    resno_ = 0;
    searched_aas_ = 0;

    subphase (SEARCHING);
    return true;
}

bool Psimscan::search_next ()
{
    xseq_ = sh_.next_aa_seq (p_->min_seq_len (), p_->max_seq_len ());
    if (!xseq_)
    {
        subphase (WRITING_RESULTS);
        return true;
    }

    searched_aas_ += xseq_->len;
    searcher_->search (*xseq_);
    results_->flush (xseq_->seq);

    if (skipped_ != sh_.skipped ())
        skipped (skipped_ = sh_.skipped ());
    if (resno_ != results_->totalStored ())
        res_no (resno_ = results_->totalStored ());

    targets_searched_ ++;
    if (p_->t_end () == -1)
        current (sh_.current ());
    else
        current (targets_searched_ + skipped_);

    return true;
}

bool Psimscan::write_results ()
{
    // results_->flush ();
    if (resno_ != results_->totalStored ())
        res_no (resno_ = results_->totalStored ());

    out_file_.open (p_->output_name ());
    if (!out_file_.is_open ())
        ers << "Unable to open output file " << p_->output_name () << Throw;

    total (1);
    current (resno_);

    switch (p_->out_mode ())
    {
        case Search_params::TAB_OUT:
            sh_.output_results_tab_aa (*results_, p_->res_per_query (), p_->print_aligns (), yseqs_, wm_, out_file_, false);
            break;
        case Search_params::TABX_OUT:
            sh_.output_results_tab_aa (*results_, p_->res_per_query (), p_->print_aligns (), yseqs_, wm_, out_file_, true);
            break;
        case Search_params::NCBI_M8_OUT:
            sh_.output_results_m8_aa (*results_, p_->res_per_query (), p_->print_aligns (), yseqs_, wm_, out_file_);
            break;
        case Search_params::NCBI_M9_OUT:
            {
            std::string header = process_name ();
            header += " ";
            header += VERSION;
            sh_.output_results_m8_aa (*results_, p_->res_per_query (), p_->print_aligns (), yseqs_, wm_, out_file_, header.c_str ());
            }
            break;
        default:
            report_header (out_file_);
            report_self (out_file_);
            sh_.output_results_aa (*results_, p_->res_per_query (), p_->print_aligns (), yseqs_, wm_, out_file_, p_->max_seq_len ());
            out_file_.close ();
    }
    return false;
}

bool Psimscan::close_handler ()
{
    return true;
}

bool Psimscan::report_results_handler ()
{
    return true;
}

// naming convention
const char* Psimscan::process_name () {return ::process_name ();}
const char* Psimscan::description () {return "Fast protein similarity search, based on PSimScan algorithm";}
const char* Psimscan::version () {return VERSION;}
const char* Psimscan::processing_item_name ()
{
    if (p_->t_end () == -1)
        return "bytes";
    else
        return "sequences";
}
const char* Psimscan::result_item_name () {return "similarities";}

// process description
struct Phase
{
    const char* name_;
    const char** subphases_;
};

const char* subphases [] = {"Loading queries", "Loading tuple distribution", "Preparing search engine", "Adding queries", "Counting tuples", "Filling lookup table", "Initializing search", "Searching similarities", "Writing results"};

Phase phases [] =
{
    {
        "Processing",
        subphases
    }
};

const char* Psimscan::phase_name (int phase_no)
{
    return phases [phase_no].name_;
}

int Psimscan::phases_total ()
{
    return sizeof (phases) / sizeof (Phase);
}

const char* Psimscan::subphase_name (int phase_no, int subphase_no)
{
    return phases [phase_no].subphases_ [subphase_no];
}

int Psimscan::subphases_total (int phase_no)
{
    return sizeof (subphases) / sizeof (const char*);
}

void Psimscan::report_state (std::ostream& o)
{
    // Process::report_state (o);
    o << "\r" << subphase_name (phase (), subphase ());
    if (subphase () == SEARCHING)
    {
        time_t cur_t = cur_time_sec ();
        double speed = double (targets_searched_) / (cur_t - subphase_start_time_sec () + 1);
        double kaaps = double (searched_aas_) / (double (cur_t - subphase_start_time_sec () + 1)*1000);
        o << ", " << targets_searched_ << " sequences searched ("
                << std::setprecision (2) << std::fixed << speed << "/sec), "
                << res_no () << " " << result_item_name () << " found."
                << " (" << std::setprecision (2) << std::fixed << percent_done () << "% of target size; "
                << kaaps << " Kaa/sec, " << results_->totalFound () << " hits)" << std::flush;
    }
}

void Psimscan::report_self (std::ostream& o)
{
    o << "Performing similarity search:" << std::endl;

    o << "Query file: " << p_->query_name ();
    if (query_set_length_) o << ", " << query_set_length_  << " sequences selected for search";
    o << std::endl;

    o << "Target file: " << p_->search_name ();
    if (search_set_length_) o << " is " << search_set_length_ / 1024 << " Kbytes long";
    if (p_->t_end () != -1) o << ", " << p_->t_end () - p_->t_beg () << " sequences selected for search";
    o << std::endl;

    if      (p_->parameters_read ())   o << "Parameters read from configuration file " << p_->in_parfname () << " :" << std::endl;
    else  o << "Default parameters used :" << std::endl;
    p_->parameters_->write (o);
    o << "Output stored to file " << p_->output_name () << std::endl << std::flush;
}

void Psimscan::report_final (std::ostream& o)
{
    time_t totaltime = time (NULL) - start_time_sec ();
    o << std::endl << "Processing took " << totaltime/3600 << " hours, " << totaltime%3600/60 << " minutes, " << totaltime%60 << " seconds." << std::endl;
    o << query_set_length_ << " sequences ("  << query_aas_ << " aa) from query file " << p_->query_name () << " were compared to" << std::endl;
    o << targets_searched_ + skipped_ << " sequences (" << searched_aas_ << " aa) from target file " << p_->search_name () << std::endl;
    if (queries_searched_ != query_set_length_)
        o << "   " << query_set_length_ - queries_searched_ << " queries were skipped due to length constraints" << std::endl;
    if (skipped_)
        o << "   " << skipped_ << " targets were skipped due to length constraints" << std::endl;
    o << resno_ << " found similarites saved to file " << p_->output_name ();
    Process::report_final (o);
}

