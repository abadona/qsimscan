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

#include "ktsearch_params.h"
#include <acc_str.h>
static const char* SEARCH_FWD_DEFAULT = TRUE_STR;
static const char* SEARCH_REV_DEFAULT = TRUE_STR;
static const char* REP_LOOKUP_DEFAULT = "2";
static const char* K_SIZE_DEFAULT = "11";
static const char* APPROX_DEFAULT = FALSE_STR;
static const char* K_THRESH_DEFAULT = "310";
static const char* MAX_SHIFT_DEFAULT = "3";
static const char* GAP_PERIOD_DEFAULT = "10";
static const char* MIN_THRESH_DEFAULT = "70";
static const char* MAX_THRESH_DEFAULT = "55";
static const char* MIN_LEN_DEFAULT = "45";
static const char* STEP_DEFAULT = "1";
static const char* KDISTR_DEFAULT = EMPTY_STR;

static const char* GIP_DEFAULT = "2.0";
static const char* GEP_DEFAULT = "0.3";
static const char* SIM_LEVEL_DEFAULT = "0.60";

static const char* MIN_SCORE_DEFAULT = "0.0";
static const char* REP_LEN_DEFAULT = "4";
static const char* REP_PERCENT_DEFAULT = "50.0";

static const char* SEARCH_FWD_HELP = "Search forward chains";
static const char* SEARCH_REV_HELP = "Search reverse chains";
static const char* SEARCH_FWD_HELP_CMDL = "Do not search forward chain";
static const char* SEARCH_REV_HELP_CMDL = "Do not search reverse chain";

static const char* REP_LOOKUP_HELP = "K-tuple redundancy lookup offset";
static const char* K_SIZE_HELP = "K-tuple size";
static const char* APPROX_HELP = "Use approximate tuple match algorithm ";
static const char* K_THRESH_HELP = "Diagonal score threshold for primary scan";
static const char* MAX_SHIFT_HELP = "Maximum diagonal shift for primary scan";
static const char* GAP_PERIOD_HELP = "Minimal period between gaps";
static const char* MIN_THRESH_HELP = "Detection threshold at MIN_LEN";
static const char* MAX_THRESH_HELP = "Detection threshold at infine length";
static const char* MIN_LEN_HELP = "Minimal length of detected similarity";
static const char* STEP_HELP = "Number of positions in subject sequence to step over at a time";
static const char* KDISTR_HELP = "Tuple distribution file (flat by default)";

static const char* GIP_HELP = "Gap initiation penalty";
static const char* GEP_HELP = "Gap extention penalty";
static const char* SIM_LEVEL_HELP = "Mismatch cost";

static const char* MIN_SCORE_HELP = "Minimal detected score";
static const char* REP_LEN_HELP = "Maximum repeat lengths for repeat filter";
static const char* REP_PERCENT_HELP = "Score threshold for repeat filter";

const char* KTSEARCH_SECTNAME = "KTSEARCH";
const char* ALIGN_SECTNAME = "ALIGN";
const char* FILTERS_SECTNAME = "FILTERS";



void KTSearch_params::add_parameters_seq_dir ()
{
    Parameter_descr PRE_SCAN_FILTERS_SECTION [] =
    {
        {"FORWARD",     BOOLEAN_STR, search_fwd_default (), search_fwd_help ()},
        {"REVERSE",     BOOLEAN_STR, search_rev_default (), search_rev_help ()}
    };
    parameters_->addSection (PRE_FILTERS_SECTNAME,     NULL,                             PRE_SCAN_FILTERS_SECTION,sizeof (PRE_SCAN_FILTERS_SECTION) / sizeof (Parameter_descr));
}
void KTSearch_params::add_parameters_ktsearch ()
{
    Parameter_descr KTSEARCH_SECTION [] =
    {
        {"REP_LOOKUP",  INTEGER_STR, rep_lookup_default (), rep_lookup_help ()},
        {"APPROX",      BOOLEAN_STR, approx_default (),     approx_help ()},
        {"K_SIZE",      INTEGER_STR, k_size_default (),     k_size_help ()},
        {"K_THRESH",    INTEGER_STR, k_thresh_default (),   k_thresh_help ()},
        {"MAX_SHIFT",   INTEGER_STR, max_shift_default (),  max_shift_help ()},
        {"GAP_PERIOD",  INTEGER_STR, gap_period_default (), gap_period_help ()},
        {"MIN_THRESH",  INTEGER_STR, min_thresh_default (), min_thresh_help ()},
        {"MAX_THRESH",  INTEGER_STR, max_thresh_default (), max_thresh_help ()},
        {"MIN_LEN",     INTEGER_STR, min_len_default (),    min_len_help ()},
        {"STEP",        INTEGER_STR, step_default (),       step_help ()},
        {"KDISTR",      STRING_STR,  kdistr_default (),     kdistr_help ()}
    };
    parameters_->addSection (KTSEARCH_SECTNAME,        "Tuple lookup search parameters", KTSEARCH_SECTION,       sizeof (KTSEARCH_SECTION) / sizeof (Parameter_descr));
}

void KTSearch_params::add_parameters_align ()
{
    Parameter_descr ALIGN_SECTION [] =
    {
        {"GIP",         FLOAT_STR,   gip_default (),        gip_help ()},
        {"GEP",         FLOAT_STR,   gep_default (),        gep_help ()},
        {"SIM_LEVEL",   FLOAT_STR,   sim_level_default (),  sim_level_help ()}
    };
    parameters_->addSection (ALIGN_SECTNAME,           "Alignment parameters",           ALIGN_SECTION,          sizeof (ALIGN_SECTION) / sizeof (Parameter_descr));
}
void KTSearch_params::add_parameters_filters ()
{
    Parameter_descr FILTERS_SECTION [] =
    {
        {"MIN_SCORE",   FLOAT_STR,   min_score_default (),  min_score_help ()},
        {"REP_LEN",     INTEGER_STR, rep_len_default (),    rep_len_help ()},
        {"REP_PERCENT", FLOAT_STR,   rep_percent_default (),rep_percent_help ()}
    };
    parameters_->addSection (FILTERS_SECTNAME,         "Filters",                        FILTERS_SECTION,        sizeof (FILTERS_SECTION) / sizeof (Parameter_descr));
}

bool KTSearch_params::prepareParameters ()
{
    bool toRet = Search_params::prepareParameters ();
    if (toRet)
    {
        add_parameters_seq_dir ();
        add_parameters_ktsearch ();
        add_parameters_align ();
        add_parameters_filters ();
    }
    return toRet;
}

bool KTSearch_params::interpreteParameters ()
{
    bool toR = Search_params::interpreteParameters ();
    if (toR)
    {
        search_forward (parameters_->getBoolean (PRE_FILTERS_SECTNAME, "FORWARD"));
        search_reverse (parameters_->getBoolean (PRE_FILTERS_SECTNAME, "REVERSE"));

        rep_lookup  (parameters_->getInteger  (KTSEARCH_SECTNAME, "REP_LOOKUP"));
        approx      (parameters_->getBoolean  (KTSEARCH_SECTNAME, "APPROX"));
        k_size      (parameters_->getInteger  (KTSEARCH_SECTNAME, "K_SIZE"));
        k_thresh    (parameters_->getInteger  (KTSEARCH_SECTNAME, "K_THRESH"));
        max_shift   (parameters_->getInteger  (KTSEARCH_SECTNAME, "MAX_SHIFT"));
        gap_period  (parameters_->getInteger  (KTSEARCH_SECTNAME, "GAP_PERIOD"));
        min_thresh  (parameters_->getInteger  (KTSEARCH_SECTNAME, "MIN_THRESH"));
        max_thresh  (parameters_->getInteger  (KTSEARCH_SECTNAME, "MAX_THRESH"));
        min_len     (parameters_->getInteger  (KTSEARCH_SECTNAME, "MIN_LEN"));
        step        ((unsigned) parameters_->getInteger  (KTSEARCH_SECTNAME, "STEP"));
        kdistr      (parameters_->getParameter(KTSEARCH_SECTNAME, "KDISTR"));

        gip         (parameters_->getFloat    (ALIGN_SECTNAME, "GIP"));
        gep         (parameters_->getFloat    (ALIGN_SECTNAME, "GEP"));
        sim_level   (parameters_->getFloat    (ALIGN_SECTNAME, "SIM_LEVEL"));

        min_score   (parameters_->getFloat    (FILTERS_SECTNAME, "MIN_SCORE"));
        rep_len     (parameters_->getInteger  (FILTERS_SECTNAME, "REP_LEN"));
        rep_percent (parameters_->getFloat    (FILTERS_SECTNAME, "REP_PERCENT"));
    }
    return toR;
}
static const char* lonofwd  []   = {"nofwd", NULL};
static const char* lonorev  []   = {"norev", NULL};
static const char* lorepk   []   = {"kred", NULL};
static const char* loapprox []   = {"approx", NULL};
static const char* loksize  []   = {"ksize", NULL};
static const char* lokthresh[]   = {"kthresh", NULL};
static const char* lomxshift[]   = {"mxshift", NULL};
static const char* logapper []   = {"gap_period", NULL};
static const char* lominthr []   = {"it", "minthr", NULL};
static const char* lomaxthr []   = {"xt", "maxthr", NULL};
static const char* lominlen []   = {"il", "minlen", NULL};
static const char* lostep   []   = {"step", NULL};
static const char* lokdistr []   = {"kdistr", NULL};
static const char* logep    []   = {"gep", NULL};
static const char* logip    []   = {"gip", NULL};
static const char* losimlev []   = {"simlev", NULL};
static const char* lominscore [] = {"minscore", NULL};
static const char* loreplen []   = {"replen", NULL};
static const char* loreplev []   = {"replev", NULL};

void KTSearch_params::add_cmdline_seq_dir ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR,  lonofwd,    "nofwd",    PRE_FILTERS_SECTNAME,  "FORWARD",  true, false, BOOLEAN_STR, inverse_bs (search_fwd_default ()), SEARCH_FWD_HELP_CMDL));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  lonorev,    "norev",    PRE_FILTERS_SECTNAME,  "REVERSE",  true, false, BOOLEAN_STR, inverse_bs (search_rev_default ()), SEARCH_REV_HELP_CMDL));
}
void KTSearch_params::add_cmdline_ktsearch ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR,  lorepk,     "replook",  KTSEARCH_SECTNAME,  "REP_LOOKUP",  true, true, INTEGER_STR, rep_lookup_default (), rep_lookup_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  loapprox,   "approx",   KTSEARCH_SECTNAME,  "APPROX",      true, false,BOOLEAN_STR, inverse_bs (approx_default ()), approx_help ()));
    keys_format_.push_back (KeyFormat ("k",        loksize,    "ksize",    KTSEARCH_SECTNAME,  "K_SIZE",      true, true, INTEGER_STR, k_size_default (), k_size_help ()));
    keys_format_.push_back (KeyFormat ("t",        lokthresh,  "kthr",     KTSEARCH_SECTNAME,  "K_THRESH",    true, true, INTEGER_STR, k_thresh_default (), k_thresh_help ()));

    keys_format_.push_back (KeyFormat (EMPTY_STR,  lomxshift,  "mxshift",  KTSEARCH_SECTNAME,  "MAX_SHIFT",   true, true, INTEGER_STR, max_shift_default (), max_shift_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  logapper,   "gapper",   KTSEARCH_SECTNAME,  "GAP_PERIOD",  true, true, INTEGER_STR, gap_period_default (), gap_period_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  lominthr,   "min_thr",  KTSEARCH_SECTNAME,  "MIN_THRESH",  true, true, INTEGER_STR, min_thresh_default (), min_thresh_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  lomaxthr,   "max_thr",  KTSEARCH_SECTNAME,  "MAX_THRESH",  true, true, INTEGER_STR, max_thresh_default (), max_thresh_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  lominlen,   "min_len",  KTSEARCH_SECTNAME,  "MIN_LEN",     true, true, INTEGER_STR, min_len_default (), min_len_help ()));

    keys_format_.push_back (KeyFormat ("q",        lostep,     "step",     KTSEARCH_SECTNAME,  "STEP",        true, true, INTEGER_STR, step_default (), step_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  lokdistr,   "kdistr",   KTSEARCH_SECTNAME,  "KDISTR",      true, true, STRING_STR,  kdistr_default (), kdistr_help ()));
}
void KTSearch_params::add_cmdline_align ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR,  logip,      "gip",      ALIGN_SECTNAME,  "GIP",            true, true, FLOAT_STR, gip_default (), gip_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  logep,      "gep",      ALIGN_SECTNAME,  "GEP",            true, true, FLOAT_STR, gep_default (), gep_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  losimlev,   "simlev",   ALIGN_SECTNAME,  "SIM_LEVEL",      true, true, FLOAT_STR, sim_level_default (), sim_level_help ()));
}
void KTSearch_params::add_cmdline_filters ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR,  lominscore, "min_score",FILTERS_SECTNAME,  "MIN_SCORE",    true, true, FLOAT_STR, min_score_default (), min_score_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  loreplen,   "rep_len",  FILTERS_SECTNAME,  "REP_LEN",      true, true, INTEGER_STR, rep_len_default (), rep_len_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR,  loreplev,   "rep_lev",  FILTERS_SECTNAME,  "REP_PERCENT",  true, true, FLOAT_STR, rep_percent_default (), rep_percent_help ()));
}

bool  KTSearch_params::prepareCmdlineFormat ()
{
    bool res = Search_params::prepareCmdlineFormat ();
    if (res)
    {
        add_cmdline_seq_dir ();
        add_cmdline_ktsearch ();
        add_cmdline_align ();
        add_cmdline_filters ();
	}
    return res;
}

const char* KTSearch_params::search_fwd_default () const
{
    return SEARCH_FWD_DEFAULT;
}

const char* KTSearch_params::search_rev_default () const
{
    return SEARCH_REV_DEFAULT;
}

const char* KTSearch_params::rep_lookup_default () const
{
    return REP_LOOKUP_DEFAULT;
}

const char* KTSearch_params::approx_default () const
{
    return APPROX_DEFAULT;
}

const char* KTSearch_params::k_size_default () const
{
    return K_SIZE_DEFAULT;
}

const char* KTSearch_params::k_thresh_default () const
{
    return K_THRESH_DEFAULT;
}

const char* KTSearch_params::max_shift_default () const
{
    return MAX_SHIFT_DEFAULT;
}

const char* KTSearch_params::gap_period_default () const
{
    return GAP_PERIOD_DEFAULT;
}

const char* KTSearch_params::min_thresh_default () const
{
    return MIN_THRESH_DEFAULT;
}

const char* KTSearch_params::max_thresh_default () const
{
    return MAX_THRESH_DEFAULT;
}

const char* KTSearch_params::min_len_default () const
{
    return MIN_LEN_DEFAULT;
}

const char* KTSearch_params::step_default () const
{
    return STEP_DEFAULT;
}

const char* KTSearch_params::kdistr_default () const
{
    return KDISTR_DEFAULT;
}

const char* KTSearch_params::gip_default () const
{
    return GIP_DEFAULT;
}

const char* KTSearch_params::gep_default () const
{
    return GEP_DEFAULT;
}

const char* KTSearch_params::sim_level_default () const
{
    return SIM_LEVEL_DEFAULT;
}

const char* KTSearch_params::min_score_default () const
{
    return MIN_SCORE_DEFAULT;
}

const char* KTSearch_params::rep_len_default () const
{
    return REP_LEN_DEFAULT;
}

const char* KTSearch_params::rep_percent_default () const
{
    return REP_PERCENT_DEFAULT;
}

//
const char* KTSearch_params::search_fwd_help () const
{
    return SEARCH_FWD_HELP;
}

const char* KTSearch_params::search_rev_help () const
{
    return SEARCH_REV_HELP;
}

const char* KTSearch_params::rep_lookup_help () const
{
    return REP_LOOKUP_HELP;
}

const char* KTSearch_params::approx_help () const
{
    return APPROX_HELP;
}

const char* KTSearch_params::k_size_help () const
{
    return K_SIZE_HELP;
}

const char* KTSearch_params::k_thresh_help () const
{
    return K_THRESH_HELP;
}

const char* KTSearch_params::max_shift_help () const
{
    return MAX_SHIFT_HELP;
}

const char* KTSearch_params::gap_period_help () const
{
    return GAP_PERIOD_HELP;
}

const char* KTSearch_params::min_thresh_help () const
{
    return MIN_THRESH_HELP;
}

const char* KTSearch_params::max_thresh_help () const
{
    return MAX_THRESH_HELP;
}

const char* KTSearch_params::min_len_help () const
{
    return MIN_LEN_HELP;
}

const char* KTSearch_params::step_help () const
{
    return STEP_HELP;
}

const char* KTSearch_params::kdistr_help () const
{
    return KDISTR_HELP;
}

const char* KTSearch_params::gip_help () const
{
    return GIP_HELP;
}

const char* KTSearch_params::gep_help () const
{
    return GEP_HELP;
}

const char* KTSearch_params::sim_level_help () const
{
    return SIM_LEVEL_HELP;
}

const char* KTSearch_params::min_score_help () const
{
    return MIN_SCORE_HELP;
}

const char* KTSearch_params::rep_len_help () const
{
    return REP_LEN_HELP;
}

const char* KTSearch_params::rep_percent_help () const
{
    return REP_PERCENT_HELP;
}
