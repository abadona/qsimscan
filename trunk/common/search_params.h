
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

#ifndef __search_params_h__
#define __search_params_h__

#include <process_params.h>

class Search_params : public Process_params
{
public:
    enum OUT_MODE
    {
        TEXT_OUT,
        TAB_OUT,
        TABX_OUT,
        NCBI_M8_OUT,
        NCBI_M9_OUT
    };

    int  max_seq_len_;
    int  min_seq_len_;
    int  max_qry_len_;
    int  q_beg_;
    int  q_end_;
    int  t_beg_;
    int  t_end_;

    bool nomerge_threads_;
    bool merge_domains_;
    bool merge_repeats_;
    int  max_dom_ovl_;
    int  max_rep_orp_;
    float gap_cap_;

    OUT_MODE out_mode_;
    bool append_;
    int  res_per_query_;
    int  print_aligns_;

    std::string query_name_;
    std::string search_name_;
    std::string output_name_;

protected:

    virtual void add_cmdline_trg_bounds ();
    virtual void add_cmdline_qry_bounds ();
    virtual void add_cmdline_set_bounds ();
    virtual void add_cmdline_seqlen_cutoffs ();
    virtual void add_cmdline_qrylen_cutoffs ();
    virtual void add_cmdline_append_control ();
    virtual void add_cmdline_dbout_control ();
    virtual void add_cmdline_res_filters ();
    virtual void add_cmdline_merge ();
    virtual void add_cmdline_args ();

    virtual void add_parameters_trg_bounds ();
    virtual void add_parameters_qry_bounds ();
    virtual void add_parameters_set_bounds ();
    virtual void add_parameters_seqlen_cutoffs ();
    virtual void add_parameters_qrylen_cutoffs ();
    virtual void add_parameters_append_control ();
    virtual void add_parameters_dbout_control ();
    virtual void add_parameters_res_filters ();
    virtual void add_parameters_merge ();
    virtual void add_parameters_args ();


    virtual bool prepareCmdlineFormat ();
    virtual bool prepareParameters ();
    virtual bool interpreteParameters ();

public:
    Search_params (const char* header = NULL, const char* procname = NULL, const char* version = NULL);

    int  min_seq_len () const {return min_seq_len_;}
    int  max_seq_len () const {return max_seq_len_;}
    int  max_qry_len () const {return max_qry_len_;}
    bool nomerge_threads () const {return nomerge_threads_;}
    bool merge_domains () const {return merge_domains_;}
    bool merge_repeats () const {return merge_repeats_;}
    int  max_dom_ovl () const {return max_dom_ovl_;}
    int  max_rep_orp () const {return max_rep_orp_;}
    float gap_cap () const {return gap_cap_;}
    int  res_per_query () const {return res_per_query_;}
    int  print_aligns () const {return print_aligns_;}
    OUT_MODE out_mode () const {return out_mode_;}
    bool append () const {return append_;}
    int q_beg () const {return q_beg_;}
    int q_end () const {return q_end_;}
    int t_beg () const {return t_beg_;}
    int t_end () const {return t_end_;}
    const char* query_name () const {return query_name_.c_str ();}
    const char* search_name () const {return search_name_.c_str ();}
    const char* output_name () const {return output_name_.c_str ();}

    void min_seq_len (int opt) {min_seq_len_ = opt;}
    void max_seq_len (int opt) {max_seq_len_ = opt;}
    void max_qry_len (int opt) {max_qry_len_ = opt;}
    void nomerge_threads (bool opt) {nomerge_threads_ = opt;}
    void merge_domains (bool opt) {merge_domains_ = opt;}
    void merge_repeats (bool opt) {merge_repeats_ = opt;}
    void max_dom_ovl (int opt) {max_dom_ovl_ = opt;}
    void max_rep_orp (int opt) {max_rep_orp_ = opt;}
    void gap_cap (float opt) {gap_cap_ = opt;}
    void res_per_query (int opt) {res_per_query_ = opt;}
    void print_aligns (int opt) {print_aligns_ = opt;}
    void out_mode (OUT_MODE opt) {out_mode_ = opt;}
    void append (bool opt) {append_ = opt;}
    void q_beg (int opt) {q_beg_ = opt;}
    void q_end (int opt) {q_end_ = opt;}
    void t_beg (int opt) {t_beg_ = opt;}
    void t_end (int opt) {t_end_ = opt;}
    void query_name (const char* opt) {query_name_ = opt;}
    void search_name (const char* opt) {search_name_ = opt;}
    void output_name (const char* opt) {output_name_ = opt;}

    virtual const char* min_seq_len_default () const;
    virtual const char* max_seq_len_default () const;
    virtual const char* max_qry_len_default () const;
    virtual const char* nomerge_threads_default () const;
    virtual const char* merge_domains_default () const;
    virtual const char* merge_repeats_default () const;
    virtual const char* max_dom_ovl_default () const;
    virtual const char* max_rep_orp_default () const;
    virtual const char* gap_cap_default () const;
    virtual const char* res_per_query_default () const;
    virtual const char* print_aligns_default () const;
    virtual const char* out_mode_default () const;
    virtual const char* append_default () const;
    virtual const char* q_beg_default () const;
    virtual const char* q_end_default () const;
    virtual const char* t_beg_default () const;
    virtual const char* t_end_default () const;
    virtual const char* query_name_default () const;
    virtual const char* search_name_default () const;
    virtual const char* output_name_default () const;

    virtual const char* min_seq_len_help () const;
    virtual const char* max_seq_len_help () const;
    virtual const char* max_qry_len_help () const;
    virtual const char* nomerge_threads_help () const;
    virtual const char* merge_domains_help () const;
    virtual const char* merge_repeats_help () const;
    virtual const char* max_dom_ovl_help () const;
    virtual const char* max_rep_orp_help () const;
    virtual const char* gap_cap_help () const;
    virtual const char* res_per_query_help () const;
    virtual const char* print_aligns_help () const;
    virtual const char* out_mode_help () const;
    virtual const char* append_help () const;
    virtual const char* q_beg_help () const;
    virtual const char* q_end_help () const;
    virtual const char* t_beg_help () const;
    virtual const char* t_end_help () const;
    virtual const char* query_name_help () const;
    virtual const char* search_name_help () const;
    virtual const char* output_name_help () const;

};

#ifndef __SEARCH_PARAMS_CPP__
extern const char* PRE_FILTERS_SECTNAME;
extern const char* SEARCH_OUTPUT_SECTNAME;
#endif


#endif // __search_params_h__
