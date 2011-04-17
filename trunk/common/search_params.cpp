
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

#include "search_params.h"
#include <acc_str.h>
#include <string.h>

#define __SEARCH_PARAMS_CPP__

extern const char* volatile_section_name;
static const char* MIN_SEQ_LEN_DEFAULT = "40"; // 40 bases. Could require overwriting for protein and SW searches
static const char* MAX_SEQ_LEN_DEFAULT = "10000000"; // 10 Mbases. Could require overwriting for protein and SW searches
static const char* MAX_QRY_LEN_DEFAULT = "10000000"; // 10 Mbases. Could require overwriting for protein and SW searches
static const char* NOMERGE_THREADS_DEFAULT = TRUE_STR;
static const char* MERGE_DOMAINS_DEFAULT = FALSE_STR;
static const char* MERGE_REPEATS_DEFAULT = FALSE_STR;
static const char* MAX_DOM_OVL_DEFAULT = "20";
static const char* MAX_REP_ORP_DEFAULT = "30";
static const char* GAP_CAP_DEFAULT = "4.0";
static const char* RES_PER_QUERY_DEFAULT = "500";
static const char* PRINT_ALIGNS_DEFAULT = "50";
static const char* OUT_MODE_DEFAULT = "TEXT";
static const char* APPEND_DEFAULT = FALSE_STR;
static const char* Q_BEG_DEFAULT = ZERO_STR;
static const char* Q_END_DEFAULT = MINUS_ONE_STR;
static const char* T_BEG_DEFAULT = ZERO_STR;
static const char* T_END_DEFAULT = MINUS_ONE_STR;
static const char* QUERY_NAME_DEFAULT = EMPTY_STR;
static const char* SEARCH_NAME_DEFAULT = EMPTY_STR;
static const char* OUTPUT_NAME_DEFAULT = EMPTY_STR;

static const char* MIN_SEQ_LEN_HELP = "Minimal length of input sequence (skip shorter ones)";
static const char* MAX_SEQ_LEN_HELP = "Maximal length of target sequence (skip longer ones)";
static const char* MAX_QRY_LEN_HELP = "Maximal length of query sequence (skip longer ones)";
static const char* NOMERGE_THREADS_HELP = "Do not merge close alternative alignments";
static const char* MERGE_DOMAINS_HELP = "Merge distant non-overlapping similarities for a sequence pair";
static const char* MERGE_REPEATS_HELP = "Select only one best representative per group of repeatitive similarities";
static const char* MAX_DOM_OVL_HELP = "Maximum overlap of merged domains";
static const char* MAX_REP_ORP_HELP = "Maximum orphan for reduced repeats";
static const char* GAP_CAP_HELP = "Gap cost limiting factor for long-range similarity merge: multiplier for gap initiaition cost";
static const char* RES_PER_QUERY_HELP = "Number of best matches to keep per query";
static const char* PRINT_ALIGNS_HELP = "Number of alignments to print in the report";
static const char* OUT_MODE_HELP = "Output mode. Valid modes are TEXT, TAB, TABX, M8, M9";
static const char* APPEND_HELP = "Append results to existing output object (or file)";
static const char* Q_BEG_HELP = "First query sequence number in input set to search";
static const char* Q_END_HELP = "Query sequence number where search stops";
static const char* T_BEG_HELP = "First target sequence number in input set to search";
static const char* T_END_HELP = "Target sequence number where search stops";
static const char* QUERY_NAME_HELP = "Query fasta file name";
static const char* SEARCH_NAME_HELP = "Target fasta file name";
static const char* OUTPUT_NAME_HELP = "Results file name";

const char* PRE_FILTERS_SECTNAME = "PRE_FILTERS";
const char* SEARCH_OUTPUT_SECTNAME = "SEARCH_OUTPUT";
const char* SIM_MERGE_SECTNAME = "SIM_MERGE";


Search_params::Search_params (const char* header, const char* procname, const char* version)
:
Process_params (header, procname, version)
{
}

static const char* loqbeg   []     = {"qbeg", NULL};
static const char* loqend   []     = {"qend", NULL};
static const char* lotbeg   []     = {"tbeg", NULL};
static const char* lotend   []     = {"tend", NULL};
static const char* lomxlen  []     = {"maxslen", NULL};
static const char* lomxqlen []     = {"maxqlen", NULL};
static const char* lomnlen  []     = {"minslen", NULL};
static const char* lorespq  []     = {"rpq", "res_per_qry", NULL};
static const char* lopral   []     = {"apq", "print_align", NULL};
static const char* looutmode[]     = {"om", "omode", "outmode", NULL};
static const char* loappend []     = {"ap", "append", NULL};
static const char* lonmthr   []     = {"nt", "nothr", NULL};
static const char* lomdom   []     = {"md", "mdom", NULL};
static const char* lomrep   []     = {"mr", "mrep", NULL};
static const char* lodovl   []     = {"do", "dovl", NULL};
static const char* lororp   []     = {"ro", "rorp", NULL};
static const char* logcap   []     = {"gc", "gcap", NULL};

void Search_params::add_cmdline_trg_bounds ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR, lotbeg,    "tbeg",    PRE_FILTERS_SECTNAME,  "T_BEG",               true, true, INTEGER_STR, t_beg_default (),      t_beg_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, lotend,    "tend",    PRE_FILTERS_SECTNAME,  "T_END",               true, true, INTEGER_STR, t_end_default (),      t_end_help ()));
}

void Search_params::add_cmdline_qry_bounds ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR, loqbeg,    "qbeg",    PRE_FILTERS_SECTNAME,  "Q_BEG",               true, true, INTEGER_STR, q_beg_default (),      q_beg_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, loqend,    "qend",    PRE_FILTERS_SECTNAME,  "Q_END",               true, true, INTEGER_STR, q_end_default (),      q_end_help ()));
}
void Search_params::add_cmdline_set_bounds ()
{
	add_cmdline_qry_bounds ();
	add_cmdline_trg_bounds ();
}
void Search_params::add_cmdline_seqlen_cutoffs ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR, lomxlen,   "maxslen",  PRE_FILTERS_SECTNAME,  "MAX_SEQ_LEN",            true, true, INTEGER_STR, max_seq_len_default (),max_seq_len_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, lomnlen,   "minslen",  PRE_FILTERS_SECTNAME,  "MIN_SEQ_LEN",            true, true, INTEGER_STR, min_seq_len_default (),min_seq_len_help ()));
    add_cmdline_qrylen_cutoffs ();
}
void Search_params::add_cmdline_qrylen_cutoffs ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR, lomxqlen,  "maxqlen",  PRE_FILTERS_SECTNAME,  "MAX_QRY_LEN",            true, true, INTEGER_STR, max_qry_len_default (),max_qry_len_help ()));
}
void Search_params::add_cmdline_args ()
{
    args_format_.push_back (ArgFormat ("QUERY_SET", "object_name or ID", query_name_help (), false));
    args_format_.push_back (ArgFormat ("SEARCH_SET", "object_name or ID", search_name_help (), false));
    args_format_.push_back (ArgFormat ("OUTPUT", "object_name", output_name_help (), false));
}
void Search_params::add_cmdline_append_control ()
{
    keys_format_.push_back (KeyFormat ("a", loappend,  "append",       SEARCH_OUTPUT_SECTNAME,  "APPEND",           true, false, BOOLEAN_STR, inverse_bs (append_default ()), append_help ()));
}
void Search_params::add_cmdline_dbout_control ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR, looutmode,   "outmode",       SEARCH_OUTPUT_SECTNAME,  "OUT_MODE",           true, true, STRING_STR, out_mode_default (), out_mode_help ()));
}
void Search_params::add_cmdline_res_filters ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR, lorespq,   "res_per_qry",  SEARCH_OUTPUT_SECTNAME,  "RES_PER_QUERY",    true, true, INTEGER_STR, res_per_query_default (), res_per_query_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, lopral,    "print_aligns", SEARCH_OUTPUT_SECTNAME,  "PRINT_ALIGNS",     true, true, INTEGER_STR, print_aligns_default (), print_aligns_help ()));

}
void Search_params::add_cmdline_merge ()
{
    keys_format_.push_back (KeyFormat (EMPTY_STR, lonmthr, "nomergethr",     SIM_MERGE_SECTNAME,  "NOMERGE_THR",      true, false, BOOLEAN_STR, inverse_bs (nomerge_threads_default ()), nomerge_threads_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, lomdom,  "mergedom",       SIM_MERGE_SECTNAME,  "MERGE_DOM",      true, false, BOOLEAN_STR, inverse_bs (merge_domains_default ()), merge_domains_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, lomrep,  "mergerep",       SIM_MERGE_SECTNAME,  "MERGE_REP",      true, false, BOOLEAN_STR, inverse_bs (merge_repeats_default ()), merge_repeats_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, lodovl,  "domovl",         SIM_MERGE_SECTNAME,  "MAX_DOM_OVL",    true, true,  INTEGER_STR, max_dom_ovl_default (), max_dom_ovl_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, lororp,  "reporp",         SIM_MERGE_SECTNAME,  "MAX_REP_ORP",    true, true,  INTEGER_STR, max_rep_orp_default (), max_rep_orp_help ()));
    keys_format_.push_back (KeyFormat (EMPTY_STR, logcap,  "gapcap",         SIM_MERGE_SECTNAME,  "GAP_CAP",        true, false, FLOAT_STR,   gap_cap_default (), gap_cap_help ()));
}

bool Search_params::prepareCmdlineFormat ()
{
    bool rv = Process_params::prepareCmdlineFormat ();
    if (rv)
    {
        add_cmdline_set_bounds ();
        add_cmdline_seqlen_cutoffs ();
        add_cmdline_merge ();
        add_cmdline_append_control ();
        add_cmdline_dbout_control ();
        add_cmdline_res_filters ();
        add_cmdline_args ();
    }
    return rv;
}

void Search_params::add_parameters_trg_bounds ()
{
    Parameter_descr PRE_SCAN_FILTERS_SECTION [] =
    {
        {"T_BEG",   INTEGER_STR, t_beg_default (), t_beg_help ()},
        {"T_END",   INTEGER_STR, t_end_default (), t_end_help ()},
    };
    parameters_->addSection (PRE_FILTERS_SECTNAME,         "Pre-scan filters",             PRE_SCAN_FILTERS_SECTION,   sizeof (PRE_SCAN_FILTERS_SECTION) / sizeof (Parameter_descr));
}
void Search_params::add_parameters_qry_bounds ()
{
    Parameter_descr PRE_SCAN_FILTERS_SECTION [] =
    {
        {"Q_BEG",   INTEGER_STR, q_beg_default (), q_beg_help ()},
        {"Q_END",   INTEGER_STR, q_end_default (), q_end_help ()},
    };
    parameters_->addSection (PRE_FILTERS_SECTNAME,         "Pre-scan filters",             PRE_SCAN_FILTERS_SECTION,   sizeof (PRE_SCAN_FILTERS_SECTION) / sizeof (Parameter_descr));
}
void Search_params::add_parameters_set_bounds ()
{
	add_parameters_qry_bounds ();
	add_parameters_trg_bounds ();
}
void Search_params::add_parameters_seqlen_cutoffs ()
{
    add_parameters_qrylen_cutoffs ();
    Parameter_descr PRE_SCAN_FILTERS_SECTION [] =
    {
        {"MAX_SEQ_LEN", INTEGER_STR, max_seq_len_default (), max_seq_len_help ()},
        {"MIN_SEQ_LEN", INTEGER_STR, min_seq_len_default (), min_seq_len_help ()},
    };
    parameters_->addSection (PRE_FILTERS_SECTNAME,         "Pre-scan filters",             PRE_SCAN_FILTERS_SECTION,   sizeof (PRE_SCAN_FILTERS_SECTION) / sizeof (Parameter_descr));
}
void Search_params::add_parameters_qrylen_cutoffs ()
{
    Parameter_descr PRE_SCAN_FILTERS_SECTION [] =
    {
        {"MAX_QRY_LEN", INTEGER_STR, max_qry_len_default (), max_qry_len_help ()},
    };
    parameters_->addSection (PRE_FILTERS_SECTNAME,         "Pre-scan filters",             PRE_SCAN_FILTERS_SECTION,   sizeof (PRE_SCAN_FILTERS_SECTION) / sizeof (Parameter_descr));
}
void Search_params::add_parameters_append_control ()
{
    Parameter_descr SEARCH_OUTPUT_SECTION [] =
    {
        {"APPEND",          BOOLEAN_STR,  append_default (), append_help ()}
    };
    parameters_->addSection (SEARCH_OUTPUT_SECTNAME,       "Search output parameters",     SEARCH_OUTPUT_SECTION,      sizeof (SEARCH_OUTPUT_SECTION) / sizeof (Parameter_descr));
}
void Search_params::add_parameters_dbout_control ()
{
    Parameter_descr SEARCH_OUTPUT_SECTION [] =
    {
        {"OUT_MODE",          STRING_STR,  out_mode_default (), out_mode_help ()},
    };
    parameters_->addSection (SEARCH_OUTPUT_SECTNAME,       "Search output parameters",     SEARCH_OUTPUT_SECTION,      sizeof (SEARCH_OUTPUT_SECTION) / sizeof (Parameter_descr));
}
void Search_params::add_parameters_res_filters ()
{
    Parameter_descr SEARCH_OUTPUT_SECTION [] =
    {
        {"RES_PER_QUERY",   INTEGER_STR,  res_per_query_default (),  res_per_query_help ()},
        {"PRINT_ALIGNS",    INTEGER_STR,  print_aligns_default (), print_aligns_help ()},
    };
    parameters_->addSection (SEARCH_OUTPUT_SECTNAME, "Search output parameters", SEARCH_OUTPUT_SECTION, sizeof (SEARCH_OUTPUT_SECTION) / sizeof (Parameter_descr));
}
void Search_params::add_parameters_merge ()
{
    Parameter_descr SIM_MERGE_SECTION [] =
    {
        {"NOMERGE_THR", BOOLEAN_STR,  nomerge_threads_default (), nomerge_threads_help ()},
        {"MERGE_DOM",   BOOLEAN_STR,  merge_domains_default (), merge_domains_help ()},
        {"MERGE_REP",   BOOLEAN_STR,  merge_repeats_default (), merge_repeats_help ()},
        {"MAX_DOM_OVL", INTEGER_STR,  max_dom_ovl_default (), max_dom_ovl_help ()},
        {"MAX_REP_ORP", INTEGER_STR,  max_rep_orp_default (), max_rep_orp_help ()},
        {"GAP_CAP",     FLOAT_STR,    gap_cap_default (), gap_cap_help ()},
    };
    parameters_->addSection (SIM_MERGE_SECTNAME, "Search output parameters", SIM_MERGE_SECTION, sizeof (SIM_MERGE_SECTION) / sizeof (Parameter_descr));
}
void Search_params::add_parameters_args ()
{
    Parameter_descr ARGUMENTS_SECTION [] =
    {
        {"QUERY_SET",   STRING_STR, query_name_default (),  query_name_help ()},
        {"SEARCH_SET",  STRING_STR, search_name_default (), search_name_help ()},
        {"OUTPUT",      STRING_STR, output_name_default (), output_name_help ()}
    };
    parameters_->addSection (volatile_section_name,  "Program arguments", ARGUMENTS_SECTION, sizeof (ARGUMENTS_SECTION) / sizeof (Parameter_descr));
}

bool Search_params::prepareParameters ()
{
    bool rv = Process_params::prepareParameters ();
    if (rv)
    {
        add_parameters_set_bounds ();
        add_parameters_seqlen_cutoffs ();
        add_parameters_append_control ();
        add_parameters_dbout_control ();
        add_parameters_merge ();
        add_parameters_res_filters ();
        add_parameters_args ();
    }
    return rv;
}

bool Search_params::interpreteParameters ()
{
    if (!Process_params::interpreteParameters ()) return false;
    max_seq_len (parameters_->getInteger (PRE_FILTERS_SECTNAME, "MAX_SEQ_LEN"));
    max_qry_len (parameters_->getInteger (PRE_FILTERS_SECTNAME, "MAX_QRY_LEN"));
    min_seq_len (parameters_->getInteger (PRE_FILTERS_SECTNAME, "MIN_SEQ_LEN"));
    q_beg (parameters_->getInteger (PRE_FILTERS_SECTNAME, "Q_BEG"));
    q_end (parameters_->getInteger (PRE_FILTERS_SECTNAME, "Q_END"));
    t_beg (parameters_->getInteger (PRE_FILTERS_SECTNAME, "T_BEG"));
    t_end (parameters_->getInteger (PRE_FILTERS_SECTNAME, "T_END"));
    nomerge_threads  (parameters_->getBoolean (SIM_MERGE_SECTNAME, "NOMERGE_THR"));
    merge_domains  (parameters_->getBoolean (SIM_MERGE_SECTNAME, "MERGE_DOM"));
    merge_repeats  (parameters_->getBoolean (SIM_MERGE_SECTNAME, "MERGE_REP"));
    max_dom_ovl    (parameters_->getInteger (SIM_MERGE_SECTNAME, "MAX_DOM_OVL"));
    max_rep_orp    (parameters_->getInteger (SIM_MERGE_SECTNAME, "MAX_REP_ORP"));
    gap_cap        (parameters_->getFloat (SIM_MERGE_SECTNAME, "GAP_CAP"));
    res_per_query  (parameters_->getInteger (SEARCH_OUTPUT_SECTNAME, "RES_PER_QUERY"));
    print_aligns   (parameters_->getInteger (SEARCH_OUTPUT_SECTNAME, "PRINT_ALIGNS"));

    const char* omode = parameters_->getParameter (SEARCH_OUTPUT_SECTNAME, "OUT_MODE");
    if (strcmp (omode, "TAB") == 0)
        out_mode (TAB_OUT);
    else if (strcmp (omode, "TABX") == 0)
        out_mode (TABX_OUT);
    else if (strcmp (omode, "M8") == 0)
        out_mode (NCBI_M8_OUT);
    else if (strcmp (omode, "M9") == 0)
        out_mode (NCBI_M9_OUT);
    else
        out_mode (TEXT_OUT);

    append         (parameters_->getBoolean (SEARCH_OUTPUT_SECTNAME, "APPEND"));
    query_name     (parameters_->getParameter (volatile_section_name, "QUERY_SET"));
    search_name    (parameters_->getParameter (volatile_section_name, "SEARCH_SET"));
    output_name    (parameters_->getParameter (volatile_section_name, "OUTPUT"));

    return true;
}

const char* Search_params::max_seq_len_default () const
{
    return MAX_SEQ_LEN_DEFAULT;
}

const char* Search_params::max_qry_len_default () const
{
    return MAX_QRY_LEN_DEFAULT;
}

const char* Search_params::min_seq_len_default () const
{
    return MIN_SEQ_LEN_DEFAULT;
}

const char* Search_params::res_per_query_default () const
{
    return RES_PER_QUERY_DEFAULT;
}

const char* Search_params::print_aligns_default () const
{
    return PRINT_ALIGNS_DEFAULT;
}

const char* Search_params::out_mode_default () const
{
    return OUT_MODE_DEFAULT;
}

const char* Search_params::append_default () const
{
    return APPEND_DEFAULT;
}

const char* Search_params::nomerge_threads_default () const
{
    return NOMERGE_THREADS_DEFAULT;
}

const char* Search_params::merge_domains_default () const
{
    return MERGE_DOMAINS_DEFAULT;
}

const char* Search_params::merge_repeats_default () const
{
    return MERGE_REPEATS_DEFAULT;
}

const char* Search_params::max_dom_ovl_default () const
{
    return MAX_DOM_OVL_DEFAULT;
}

const char* Search_params::max_rep_orp_default () const
{
    return MAX_REP_ORP_DEFAULT;
}

const char* Search_params::gap_cap_default () const
{
    return GAP_CAP_DEFAULT;
}

const char* Search_params::query_name_default () const
{
    return QUERY_NAME_DEFAULT;
}

const char* Search_params::search_name_default () const
{
    return SEARCH_NAME_DEFAULT;
}

const char* Search_params::output_name_default () const
{
    return OUTPUT_NAME_DEFAULT;
}

const char* Search_params::q_beg_default () const
{
    return Q_BEG_DEFAULT;
}

const char* Search_params::q_end_default () const
{
    return Q_END_DEFAULT;
}

const char* Search_params::t_beg_default () const
{
    return T_BEG_DEFAULT;
}

const char* Search_params::t_end_default () const
{
    return T_END_DEFAULT;
}


//
const char* Search_params::max_seq_len_help () const
{
    return MAX_SEQ_LEN_HELP;
}

const char* Search_params::max_qry_len_help () const
{
    return MAX_QRY_LEN_HELP;
}

const char* Search_params::min_seq_len_help () const
{
    return MIN_SEQ_LEN_HELP;
}

const char* Search_params::res_per_query_help () const
{
    return RES_PER_QUERY_HELP;
}

const char* Search_params::print_aligns_help () const
{
    return PRINT_ALIGNS_HELP;
}

const char* Search_params::out_mode_help () const
{
    return OUT_MODE_HELP;
}

const char* Search_params::append_help () const
{
    return APPEND_HELP;
}

const char* Search_params::nomerge_threads_help () const
{
    return NOMERGE_THREADS_HELP;
}

const char* Search_params::merge_domains_help () const
{
    return MERGE_DOMAINS_HELP;
}

const char* Search_params::merge_repeats_help () const
{
    return MERGE_REPEATS_HELP;
}

const char* Search_params::gap_cap_help () const
{
    return GAP_CAP_HELP;
}

const char* Search_params::max_dom_ovl_help () const
{
    return MAX_DOM_OVL_HELP;
}

const char* Search_params::max_rep_orp_help () const
{
    return MAX_REP_ORP_HELP;
}

const char* Search_params::query_name_help () const
{
    return QUERY_NAME_HELP;
}

const char* Search_params::search_name_help () const
{
    return SEARCH_NAME_HELP;
}

const char* Search_params::output_name_help () const
{
    return OUTPUT_NAME_HELP;
}

const char* Search_params::q_beg_help () const
{
    return Q_BEG_HELP;
}

const char* Search_params::q_end_help () const
{
    return Q_END_HELP;
}

const char* Search_params::t_beg_help () const
{
    return T_BEG_HELP;
}

const char* Search_params::t_end_help () const
{
    return T_END_HELP;
}


