
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

#include "psimscan_params.h"
#include <acc_str.h>

extern const char* VERSION;
static const char* HEADER = "Tool for searching similarities in protein databases.\nBased on QSIMSCAN algorithm by SciDM";

// defaults

static const char* MIN_SEQ_LEN_DEFAULT = "10";    // 10 aa.
static const char* MAX_SEQ_LEN_DEFAULT = "50000"; // 50 Kaa.
static const char* MAX_QRY_LEN_DEFAULT = "50000"; // 50 Kaa.

static const char MAX_QUERIES_NUMBER_DEFAULT [] = "6000";
static const char MAX_TOTAL_QUERIES_LEN_DEFAULT [] = "6000000";

static const char K_SIZE_DEFAULT [] = "5";        // tuple size
static const char APPROX_DEFAULT [] = "0.84";     // maximal diversification score - fraction of self-score
static const char K_THRESH_DEFAULT [] = "18.0";   // initial hit threshold (in self-correlation averages)
static const char MAX_SHIFT_DEFAULT [] = "3";     // maximum detectible gap size
static const char KDISTR_NAME_DEFAULT [] = "";    // name for the object containing ktuple distribution. Unused for now
static const char EXTEND_BAND_DEFAULT [] = "1.5"; // the multiplier for max. undetectible gap size for calculating band extension
static const char WIDEN_BAND_DEFAULT [] = "3";    // the number of diagonals to the sides to leftmost and rightmost hits to check in band alignment
static const char DIST_FACT_DEFAULT [] = "0.05";  // multiplier for inter-match distance cost

static const char MATRIX_NAME_DEFAULT [] = "BLOSUM62";
static const char GEP_DEFAULT [] = "0.3";         // in matrix averages
static const char GIP_DEFAULT [] = "2.5";         // in matrix averages

static const char RES_PER_QUERY_DEFAULT [] = "5000";
static const char PRINT_ALIGNS_DEFAULT [] = "50";
static const char DB_OUT_DEFAULT [] = "FALSE";
static const char APPEND_DEFAULT [] = "False";

static const char EVAL_EVAL_DEFAULT [] = "True";
static const char MIN_SCORE_DEFAULT [] = "0.0";   // the absolute score
static const char MIN_LEN_DEFAULT [] = "0";      // aminoacids


// Help

static const char MAX_QUERIES_NUMBER_HELP [] = "Maximum expected number of queries in single search";
static const char MAX_TOTAL_QUERIES_LEN_HELP [] = "Maximum expected sum of query sequence lengths";

static const char K_SIZE_HELP [] = "tuple size";
static const char APPROX_HELP [] = "Maximum diversification (fraction of self-score)";
static const char K_THRESH_HELP [] = "initial hit threshold (in self-correlation averages)";
static const char MAX_SHIFT_HELP [] = "maximum detectible gap size";
static const char KDISTR_NAME_HELP [] = "name for the object containing ktuple distribution";
static const char EXTEND_BAND_HELP [] = "band extension (multiplier for maximal undetectible gap size)";
static const char WIDEN_BAND_HELP [] = "number of diagonals to the sides to leftmost and rightmost hits to check in band alignment";
static const char DIST_FACT_HELP [] = "multiplier for inter-match distance cost";

static const char MATRIX_NAME_HELP [] = "BLOSUM62";
static const char GEP_HELP [] = "Gap extension penalty (multiplier for matrix average match)";
static const char GIP_HELP [] = "Gap initiation penalty (multiplier for matrix average match)";         // in matrix averages

static const char EVAL_EVAL_HELP [] = "Evalue evaluation using Karlin-Altschul statistics";
static const char MIN_SCORE_HELP [] = "Minimal similarity score";
static const char MIN_LEN_HELP [] = "Minimal similarity length";


static const char QUERY_NAME_HELP [] = "Query : set of protein sequences (workspace name for NSeq{1} object)";
static const char SEARCH_NAME_HELP [] = "Target : set of protein sequences (workspace name for NSeq{1} object)";
static const char OUTPUT_NAME_HELP [] = "Output file name or output name in database (workspace name for Sim{1} object)";


// section names

const char* KTSEARCH_SECTNAME = "KTSEARCH";
const char* POST_FILTERS_SECTNAME = "POST_FILTERS";
const char* SWSEARCH_SECTNAME = "SWSEARCH";

Process_params* process_params_factory ()
{
    return new Psimscan_params ();
}

Psimscan_params::Psimscan_params ()
:
Search_params (HEADER, NULL, VERSION)
{
}

static const char* lomaxq []     = {"maxq", NULL};
static const char* lototqlen []  = {"totqlen", NULL};
static const char* loapprox []   = {"approx", NULL};
static const char* loksize []    = {"ksize", NULL};
static const char* lokthresh []  = {"kthresh", NULL};
static const char* lomxshift []  = {"mxshift", NULL};
static const char* loext []      = {"ext", NULL};
static const char* lowiden []    = {"widen", NULL};
static const char* lodistfact [] = {"distfact", NULL};
static const char* lonoeval []   = {"noeval", NULL};
static const char* lominlen []   = {"minlen", NULL};
static const char* lominsc []    = {"minsc", NULL};
static const char* logep []      = {"gep", NULL};
static const char* logip []      = {"gip", NULL};
static const char* lomatrix []   = {"matrix", NULL};

bool Psimscan_params::prepareCmdlineFormat ()
{
    bool toR = Search_params::prepareCmdlineFormat ();
    if (toR)
    {
        keys_format_.push_back (KeyFormat ("", lomaxq,    "maxq",    PRE_FILTERS_SECTNAME,  "MAX_QUERIES_NUMBER",  true, true, INTEGER_STR, max_queries_number_default (), max_queries_number_help ()));
        keys_format_.push_back (KeyFormat ("", lototqlen, "totqlen", PRE_FILTERS_SECTNAME,  "MAX_TOTAL_QUERIES_LEN",true,true, INTEGER_STR, max_total_queries_len_default (),  max_total_queries_len_help ()));

        keys_format_.push_back (KeyFormat ("", loapprox,  "approx",  KTSEARCH_SECTNAME,     "APPROX",              true, true, FLOAT_STR,   approx_default (), approx_help ()));
        keys_format_.push_back (KeyFormat ("", loksize,   "ksize",   KTSEARCH_SECTNAME,     "K_SIZE",              true, true, INTEGER_STR, k_size_default (), k_size_help ()));
        keys_format_.push_back (KeyFormat ("", lokthresh, "kthresh", KTSEARCH_SECTNAME,     "K_THRESH",            true, true, FLOAT_STR,   k_thresh_default (), k_thresh_help ()));
        keys_format_.push_back (KeyFormat ("", lomxshift, "mxshift", KTSEARCH_SECTNAME,     "MAX_SHIFT",           true, true, INTEGER_STR, max_shift_default (), max_shift_help ()));
        keys_format_.push_back (KeyFormat ("", loext,     "ext",     KTSEARCH_SECTNAME,     "EXTEND_BAND",         true, true, FLOAT_STR,   extend_band_default (), extend_band_help ()));
        keys_format_.push_back (KeyFormat ("", lowiden,   "widen",   KTSEARCH_SECTNAME,     "WIDEN_BAND",          true, true, INTEGER_STR, widen_band_default (), widen_band_help ()));
        keys_format_.push_back (KeyFormat ("", lodistfact,"distfact",KTSEARCH_SECTNAME,     "DIST_FACT",           true, true, FLOAT_STR,   dist_fact_default (), dist_fact_help ()));
        keys_format_.push_back (KeyFormat ("", lonoeval,  "noeval",  POST_FILTERS_SECTNAME, "EVAL_EVAL",           true, false,BOOLEAN_STR, eval_eval_default (), eval_eval_help ()));
        keys_format_.push_back (KeyFormat ("", lominlen,  "minlen",  POST_FILTERS_SECTNAME, "MIN_LEN",             true, true, INTEGER_STR, min_len_default (), min_len_help ()));
        keys_format_.push_back (KeyFormat ("", lominsc,   "minsc",   POST_FILTERS_SECTNAME, "MIN_SCORE",           true, true, FLOAT_STR,   min_score_default (), min_score_help ()));
        keys_format_.push_back (KeyFormat ("", logep,     "gep",     SWSEARCH_SECTNAME,     "GEP",                 true, true, FLOAT_STR,   gep_default (), gep_help ()));
        keys_format_.push_back (KeyFormat ("", logip,     "gip",     SWSEARCH_SECTNAME,     "GIP",                 true, true, FLOAT_STR,   gip_default (), gip_help ()));
        keys_format_.push_back (KeyFormat ("", lomatrix,  "matrix",  SWSEARCH_SECTNAME,     "MATRIX_NAME",         true, true, STRING_STR,  matrix_name_default (), matrix_name_help ()));
    }
    return toR;
}

bool Psimscan_params::prepareParameters ()
{
    static Parameter_descr PRE_SCAN_FILTERS_SECTION [] =
    {
        {"MAX_QUERIES_NUMBER", INTEGER_STR, max_queries_number_default (), max_queries_number_help ()},
        {"MAX_TOTAL_QUERIES_LEN", INTEGER_STR, max_total_queries_len_default (), max_total_queries_len_help ()}
    };
    static Parameter_descr KTSEARCH_SECTION [] =
    {
        {"APPROX", FLOAT_STR, approx_default (), approx_help ()},
        {"K_SIZE", INTEGER_STR, k_size_default (), k_size_help ()},
        {"K_THRESH", FLOAT_STR, k_thresh_default (), k_thresh_help ()},
        {"MAX_SHIFT", INTEGER_STR, max_shift_default (), max_shift_help ()},
        {"KDISTR_NAME", STRING_STR, kdistr_name_default (), kdistr_name_help ()},
        {"EXTEND_BAND", FLOAT_STR, extend_band_default (), extend_band_help ()},
        {"WIDEN_BAND", INTEGER_STR, widen_band_default (), widen_band_help ()},
        {"DIST_FACT", FLOAT_STR, dist_fact_default (), dist_fact_help ()}
    };
    static Parameter_descr POST_SCAN_FILTERS_SECTION [] =
    {
        {"EVAL_EVAL", BOOLEAN_STR, eval_eval_default (), eval_eval_help ()},
        {"MIN_LEN", INTEGER_STR, min_len_default (), min_len_help ()},
        {"MIN_SCORE", FLOAT_STR, min_score_default (), min_score_help ()},
    };
    static Parameter_descr SWSEARCH_SECTION [] =
    {
        {"GIP",         FLOAT_STR,    gip_default (),            gip_help ()},
        {"GEP",         FLOAT_STR,    gep_default (),            gep_help ()},
        {"MATRIX_NAME", STRING_STR,   matrix_name_default (),    matrix_name_help ()}
    };

    bool toR = Search_params::prepareParameters ();
    if (toR)
    {
        parameters_->addSection (PRE_FILTERS_SECTNAME,     "Pre-scan sequence filtering",      PRE_SCAN_FILTERS_SECTION,  sizeof (PRE_SCAN_FILTERS_SECTION) / sizeof (Parameter_descr));
        parameters_->addSection (KTSEARCH_SECTNAME,        "Tuple lookup search parameters",   KTSEARCH_SECTION,          sizeof (KTSEARCH_SECTION) / sizeof (Parameter_descr));
        parameters_->addSection (POST_FILTERS_SECTNAME,    "Similarity filters",               POST_SCAN_FILTERS_SECTION, sizeof (POST_SCAN_FILTERS_SECTION) / sizeof (Parameter_descr));
        parameters_->addSection (SWSEARCH_SECTNAME,        "Edit distance parameters",         SWSEARCH_SECTION,          sizeof (SWSEARCH_SECTION) / sizeof (Parameter_descr));
    }
    return toR;
}

bool Psimscan_params::interpreteParameters ()
{
    bool toR = Search_params::interpreteParameters ();
    if (toR)
    {
        max_queries_number (parameters_->getInteger (PRE_FILTERS_SECTNAME, "MAX_QUERIES_NUMBER"));
        max_total_queries_len (parameters_->getInteger (PRE_FILTERS_SECTNAME,  "MAX_TOTAL_QUERIES_LEN"));

        approx      (parameters_->getFloat     (KTSEARCH_SECTNAME, "APPROX"));
        k_size      (parameters_->getInteger   (KTSEARCH_SECTNAME, "K_SIZE"));
        k_thresh    (parameters_->getFloat     (KTSEARCH_SECTNAME, "K_THRESH"));
        max_shift   (parameters_->getInteger   (KTSEARCH_SECTNAME, "MAX_SHIFT"));
        kdistr_name (parameters_->getParameter (KTSEARCH_SECTNAME, "KDISTR_NAME"));
        extend_band (parameters_->getFloat     (KTSEARCH_SECTNAME, "EXTEND_BAND"));
        widen_band  (parameters_->getInteger   (KTSEARCH_SECTNAME, "WIDEN_BAND"));
        dist_fact   (parameters_->getFloat     (KTSEARCH_SECTNAME, "DIST_FACT"));

        min_len     (parameters_->getInteger   (POST_FILTERS_SECTNAME, "MIN_LEN"));
        min_score   (parameters_->getFloat     (POST_FILTERS_SECTNAME, "MIN_SCORE"));
        eval_eval   (parameters_->getBoolean   (POST_FILTERS_SECTNAME, "EVAL_EVAL"));

        matrix_name  (parameters_->getParameter (SWSEARCH_SECTNAME, "MATRIX_NAME"));
        gip          (parameters_->getFloat     (SWSEARCH_SECTNAME, "GIP"));
        gep          (parameters_->getFloat     (SWSEARCH_SECTNAME, "GEP"));
    }
    return toR;
}

const char* Psimscan_params::max_queries_number_default () const
{
    return MAX_QUERIES_NUMBER_DEFAULT;
}

const char* Psimscan_params::max_total_queries_len_default () const
{
    return MAX_TOTAL_QUERIES_LEN_DEFAULT;
}

const char* Psimscan_params::approx_default () const
{
    return APPROX_DEFAULT;
}

const char* Psimscan_params::k_size_default () const
{
    return K_SIZE_DEFAULT;
}

const char* Psimscan_params::k_thresh_default () const
{
    return K_THRESH_DEFAULT;
}

const char* Psimscan_params::max_shift_default () const
{
    return MAX_SHIFT_DEFAULT;
}

const char* Psimscan_params::kdistr_name_default () const
{
    return KDISTR_NAME_DEFAULT;
}

const char* Psimscan_params::extend_band_default () const
{
    return EXTEND_BAND_DEFAULT;
}

const char* Psimscan_params::widen_band_default () const
{
    return WIDEN_BAND_DEFAULT;
}

const char* Psimscan_params::dist_fact_default () const
{
    return DIST_FACT_DEFAULT;
}


const char* Psimscan_params::matrix_name_default () const
{
    return MATRIX_NAME_DEFAULT;
}

const char* Psimscan_params::gip_default () const
{
    return GIP_DEFAULT;
}

const char* Psimscan_params::gep_default () const
{
    return GEP_DEFAULT;
}

const char* Psimscan_params::eval_eval_default () const
{
    return EVAL_EVAL_DEFAULT;
}

const char* Psimscan_params::min_score_default () const
{
    return MIN_SCORE_DEFAULT;
}

const char* Psimscan_params::min_len_default () const
{
    return MIN_LEN_DEFAULT;
}

// HELP

const char* Psimscan_params::max_queries_number_help () const
{
    return MAX_QUERIES_NUMBER_HELP;
}

const char* Psimscan_params::max_total_queries_len_help () const
{
    return MAX_TOTAL_QUERIES_LEN_HELP;
}

const char* Psimscan_params::approx_help () const
{
    return APPROX_HELP;
}

const char* Psimscan_params::k_size_help () const
{
    return K_SIZE_HELP;
}

const char* Psimscan_params::k_thresh_help () const
{
    return K_THRESH_HELP;
}

const char* Psimscan_params::max_shift_help () const
{
    return MAX_SHIFT_HELP;
}

const char* Psimscan_params::kdistr_name_help () const
{
    return KDISTR_NAME_HELP;
}

const char* Psimscan_params::extend_band_help () const
{
    return EXTEND_BAND_HELP;
}

const char* Psimscan_params::widen_band_help () const
{
    return WIDEN_BAND_HELP;
}

const char* Psimscan_params::dist_fact_help () const
{
    return DIST_FACT_HELP;
}

const char* Psimscan_params::matrix_name_help () const
{
    return MATRIX_NAME_HELP;
}

const char* Psimscan_params::gip_help () const
{
    return GIP_HELP;
}

const char* Psimscan_params::gep_help () const
{
    return GEP_HELP;
}

const char* Psimscan_params::eval_eval_help () const
{
    return EVAL_EVAL_HELP;
}

const char* Psimscan_params::min_score_help () const
{
    return MIN_SCORE_HELP;
}

const char* Psimscan_params::min_len_help () const
{
    return MIN_LEN_HELP;
}

const char* Psimscan_params::query_name_help () const
{
    return QUERY_NAME_HELP;
}

const char* Psimscan_params::search_name_help () const
{
    return SEARCH_NAME_HELP;
}

const char* Psimscan_params::output_name_help () const
{
    return OUTPUT_NAME_HELP;
}

const char* Psimscan_params::min_seq_len_default () const
{
    return MIN_SEQ_LEN_DEFAULT;
}

const char* Psimscan_params::max_seq_len_default () const
{
    return MAX_SEQ_LEN_DEFAULT;
}

const char* Psimscan_params::max_qry_len_default () const
{
    return MAX_QRY_LEN_DEFAULT;
}

