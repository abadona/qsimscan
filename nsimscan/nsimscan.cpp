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

#pragma warning (disable: 4786)
#include "nsimscan.h"
#include <fileutils.h>
#include <print_batches.h>
#include <biosequence.h>
#include <resource.h>
#include <weights.h>
#include <blast_results_batch.h>
#include <kt_srch_gx.h>
#include <i64out.h>
#include <iomanip>
#include <sim_merger.h>

const char* VERSION = "1.1.76 (March 2016)";

Process* process_factory ()
{
    return new Nsimscan ();
}

const char* process_name ()
{
    return "NSIMSCAN";
}

Nsimscan::Nsimscan ()
:
wm_ (NULL),
results_ (NULL),
searcher_ (NULL),
query_set_length_ (0),
search_set_length_ (0),
query_bps_ (0),
searched_bps_ (0),
queries_searched_ (0),
targets_searched_ (0),
skipped_ (0),
resno_ (0),
max_qry_len_ (0)
{
}

Nsimscan::~Nsimscan ()
{
    close_handler ();
}

bool Nsimscan::init_handler (Process_params* params)
{
    p_ = (Nsimscan_params*) params;
    // make sure the output is ready to receive data
    checkCreateOutputFile (p_->output_name (), p_->overwrite (), p_->append ());

    // check for existence of query set
    if (!file_exists (p_->query_name ()))
        ers << "ERROR: File " << p_->query_name () << " not found" << Throw;

    // check for existence of search set (if specified)
    if (!file_exists (p_->search_name ()))
        ers << "ERROR: File " << p_->search_name () << " not found" << Throw;

    if (!p_->search_forward () && !p_->search_reverse ())
        ers << "ERROR: both directions are skipped, Nothing to do" << Throw;

    total (0);
    current (0);

    query_set_length_ = 0;
    search_set_length_ = 0;
    query_bps_ = 0;
    searched_bps_ = 0;
    queries_searched_ = 0;
    targets_searched_ = 0;
    skipped_ = 0;
    resno_ = 0;
    max_qry_len_ = 0;

    subphase (LOADING_QUERIES);
    return true;
}

bool Nsimscan::next_handler ()
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

bool Nsimscan::load_queries ()
{
    // load queries
    queries_searched_ = sh_.load_queries_nn (p_->query_name (), f_yseqs_, p_->min_seq_len (), p_->max_seq_len (), &query_set_length_, (unsigned) p_->q_beg (), (unsigned) p_->q_end ());

    // load reverse if needed; count number of bases in queries
    query_bps_ = 0;
    max_qry_len_ = 0;
    for (int iii = 0; iii < f_yseqs_.size (); iii ++)
    {
        int clen = f_yseqs_ [iii].len;
        if (max_qry_len_ < clen)
            max_qry_len_ = clen;
        query_bps_ += clen;
    }
    if (p_->search_reverse ())
    {
        reverse_seqs (f_yseqs_, r_yseqs_);
    }
    query_bps_ *= (p_->search_reverse () + p_->search_forward ());

    // report skipped queries
    if (p_->verbose () && queries_searched_ != query_set_length_)
        std::clog << query_set_length_ - queries_searched_ << " sequences were skipped due to length constraints" << std::flush;

    total (1);
    current (0);
    subphase (LOADING_TUPLE_WEIGHTS);
    return true;
}


bool Nsimscan::load_tuple_weights ()
{
    // load the k-tuple distribution
    if (!*p_->kdistr ())
        k_distr_ = sh_.create_ktuple_distribution (p_->k_size (), p_->sim_level () * 100, mismatch_score_);
    else
        k_distr_ = sh_.load_ktuple_distribution (p_->kdistr (), p_->sim_level () * 100, mismatch_score_);
    subphase (PREPARING_SEARCHER);
    return true;
}


bool Nsimscan::prepare_searcher ()
{
    // create the weight matrix
    wm_ = new WMatrix (int (p_->sim_level () * 100), p_->gip (), p_->gep ());

    // do we need score threshold here?
    sim_merger_ = new SimMerger (!p_->nomerge_threads (), p_->merge_repeats (), p_->merge_domains (), *wm_, p_->gip (), p_->gep (), p_->gap_cap ()*p_->gip (), p_->max_dom_ovl (), p_->max_rep_orp ());
    results_ = new BlastResultsBatch (p_->res_per_query (), p_->res_per_target (), p_->rep_len (), int (p_->rep_percent ()), false, *sim_merger_);

    // create the searcher
    int total_queries_count = query_set_length_ * (p_->search_reverse () ? 2 : 1);
    searcher_ = new KT_SEARCH (p_->k_size (), k_distr_, total_queries_count, query_bps_, p_->max_seq_len (), results_);

    // set search parameters
    int gecost = (int) p_->gep (); if (gecost < 1) gecost = 1; gecost *= mismatch_score_;
    int gicost = (int) p_->gip (); if (gicost < 1) gicost = 1; gicost *= mismatch_score_;

    searcher_->rep_del       = p_->rep_lookup ();
    searcher_->apprx_match   = p_->approx ();
    searcher_->k_thresh      = p_->k_thresh ();
    searcher_->b_thresh_s    = p_->min_thresh ();
    searcher_->b_thresh_l    = p_->max_thresh ();
    searcher_->l_thresh      = p_->min_len ();
    searcher_->max_offs      = p_->max_shift ();
    searcher_->g_period      = p_->gap_period ();
    searcher_->xstep         = p_->step ();
    searcher_->k_gep         = gecost;
    searcher_->k_gip         = gicost;

    subphase (ADDING_QUERIES);
    return true;
}

bool Nsimscan::add_queries ()
{
    // fill y memory
    for (int yno = 0; yno < f_yseqs_.size (); yno ++)
    {
        // load forward
        if (p_->search_forward ())
            if (!searcher_->add_yseq (f_yseqs_ [yno]))
                ERR("Error loading forward Y sequence");
        // load revers if needed
        if (p_->search_reverse ())
            if (!searcher_->add_yseq (r_yseqs_ [yno]))
                ERR("Error loading reverse Y sequence");
    }
    subphase (COUNTING_TUPLES);
    return true;
}

bool Nsimscan::count_tuples ()
{
    subphase (FILLING_LOOKUP);
    return true;
}

bool Nsimscan::fill_lookup ()
{
    subphase (INIT_SEARCH);
    return true;
}

bool Nsimscan::init_search ()
{
    // init search
    if (!searcher_->init ()) ERR ("Unable to initialize KT scanner");
    sh_.init_nn_search (p_->search_name (), (unsigned) p_->t_beg (), (unsigned) p_->t_end ());

    if (p_->t_end () == -1)
        total (sh_.total ());
    else
        total (p_->t_end () - p_->t_beg ());

    current (0);
    skipped_ = 0;
    resno_ = 0;
    searched_bps_ = 0;

    subphase (SEARCHING);
    return true;
}

bool Nsimscan::search_next ()
{
    xseq_ = sh_.next_nn_seq (p_->min_seq_len (), p_->max_seq_len ());
    if (!xseq_)
    {
        subphase (WRITING_RESULTS);
        return true;
    }

    searched_bps_ += xseq_->len;
    if (!searcher_->search (*xseq_))
        ERR ("Search engine error");
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

bool Nsimscan::write_results ()
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
            sh_.output_results_tab_nn (*results_, p_->res_per_query (), p_->print_aligns (), f_yseqs_, r_yseqs_, wm_, out_file_, false);
            break;
        case Search_params::TABX_OUT:
            sh_.output_results_tab_nn (*results_, p_->res_per_query (), p_->print_aligns (), f_yseqs_, r_yseqs_, wm_, out_file_, true);
            break;
        case Search_params::NCBI_M8_OUT:
            sh_.output_results_m8_nn (*results_, p_->res_per_query (), p_->print_aligns (), f_yseqs_, r_yseqs_, wm_, out_file_);
            break;
        case Search_params::NCBI_M9_OUT:
            {
            std::string header = process_name ();
            header += " ";
            header += VERSION;
            sh_.output_results_m8_nn (*results_, p_->res_per_query (), p_->print_aligns (), f_yseqs_, r_yseqs_, wm_, out_file_, header.c_str ());
            }
            break;
        default:
            report_header (out_file_);
            report_self (out_file_);
            sh_.output_results_nn (*results_, p_->res_per_query (), p_->print_aligns (), f_yseqs_, r_yseqs_, wm_, out_file_, p_->max_seq_len ());
            out_file_.close ();
    }
    return false;
}

bool Nsimscan::close_handler ()
{
    return true;
}

bool Nsimscan::report_results_handler ()
{
    return true;
}

// naming convention
const char* Nsimscan::process_name () {return ::process_name ();}
const char* Nsimscan::description () {return "Fast nucleotide sequence similarity search, based on YABLAST algorithm";}
const char* Nsimscan::version () {return VERSION;}
const char* Nsimscan::processing_item_name ()
{
    if (p_->t_end () == -1)
        return "bytes";
    else
        return "sequences";
}
const char* Nsimscan::result_item_name () {return "similarities";}

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

const char* Nsimscan::phase_name (int phase_no)
{
    return phases [phase_no].name_;
}

int Nsimscan::phases_total ()
{
    return sizeof (phases) / sizeof (Phase);
}

const char* Nsimscan::subphase_name (int phase_no, int subphase_no)
{
    return phases [phase_no].subphases_ [subphase_no];
}

int Nsimscan::subphases_total (int phase_no)
{
    return sizeof (subphases) / sizeof (const char*);
}

void Nsimscan::report_state (std::ostream& o)
{
    // Process::report_state (o);
    o << "\r" << subphase_name (phase (), subphase ());
    if (subphase () == SEARCHING)
    {
        time_t cur_t = cur_time_sec ();
        double speed = double (targets_searched_) / (cur_t - subphase_start_time_sec () + 1);
        double kbpps = double (searched_bps_) / (double (cur_t - subphase_start_time_sec () + 1)*1000);
        o << ", " << targets_searched_ << " sequences searched (" << std::setprecision (2) << std::fixed << percent_done () << "%, "
                << std::setprecision (2) << std::fixed << speed << " /s, " << kbpps << " Kb/s) "
                << res_no () << " " << result_item_name () << " found "
                << "[" << searcher_->k_matches << " diag, " << results_->totalFound () << " hits";
        if (p_->rep_len ())
            o << ", " << results_->passedRepeats () << " nr";
        o << ", "<< results_->tot_merged () << " merged]" << std::flush;
    }
}

void Nsimscan::report_self (std::ostream& o)
{
    o << "Performing similarity search:" << std::endl;

    o << "Query file: " << p_->query_name ();
    if (query_set_length_) o << ", " << query_set_length_  << " sequences selected for search";
    o << std::endl;

    o << "Target file: " << p_->search_name ();
    if (search_set_length_) o << " is " << search_set_length_ / 1024 << " Kbytes long";
    if (p_->t_end () != -1) o << ", " << p_->t_end () - p_->t_beg () << " sequences selected for search";
    o << std::endl;

    if (p_->parameters_read ())   o << "Parameters base values are loaded from the configuration file " << p_->in_parfname () <<  std::endl;
    o << "Parameters used :" << std::endl;
    p_->parameters_->write (o);
    o << "Output stored to file " << p_->output_name () << std::endl << std::flush;
}

void Nsimscan::report_final (std::ostream& o)
{
    time_t totaltime = time (NULL) - start_time_sec ();
    o << std::endl << "Processing took " << totaltime/3600 << " hours, " << totaltime%3600/60 << " minutes, " << totaltime%60 << " seconds." << std::endl;
    o << query_set_length_ << " sequences ";
    if (!p_->search_reverse ())
        o << ", forward chains ";
    else if (!p_->search_forward ())
        o << ", reverse chains ";
    else
        o << ", forward and reverse chains ";
    o << "("  << query_bps_ << " bp) from query file " << p_->query_name () << " were compared to" << std::endl;
    o << targets_searched_ + skipped_ << " sequences (" << searched_bps_ << " bp) from target file " << p_->search_name () << std::endl;
    if (queries_searched_ != query_set_length_)
        o << "   " << query_set_length_ - queries_searched_ << " queries were skipped due to length constraints" << std::endl;
    if (skipped_)
        o << "   " << skipped_ << " targets were skipped due to length constraints" << std::endl;
    o << resno_ << " found similarites saved to file " << p_->output_name () << std::endl;
    o << "    " << searcher_->k_matches << " candidate locations with score over " << p_->k_thresh () << " seen" << std::endl;
    o << "    " << searcher_->b_matches << " high-scoring alignment segments found" << std::endl;
    if (p_->rep_len ())
        o << "    " << results_->passedRepeats () << " passed repeats / redundancy filter" << std::endl;
    o << "    " << results_->tot_merged () << " similarites after segment / thread merging and repeat selection" << std::endl
      << "    " << resno_ << " best similarites selected" << std::endl;
    Process::report_final (o);
}


