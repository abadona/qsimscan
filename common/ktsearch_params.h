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

#ifndef __ktsearch_params_h__
#define __ktsearch_params_h__

#include "search_params.h"

class KTSearch_params : public Search_params
{
    bool search_forward_;
    bool search_reverse_;

    int rep_lookup_;
    bool approx_;
    int k_size_;
    int k_thresh_;
    int max_shift_;
    int gap_period_;
    int min_thresh_;
    int max_thresh_;
    int min_len_;
    unsigned step_;
    std::string kdistr_;

    double gip_;
    double gep_;
    double sim_level_;

    double min_score_;
    int rep_len_;
    double rep_percent_;

protected:

    virtual void add_parameters_seq_dir ();
    virtual void add_parameters_ktsearch ();
    virtual void add_parameters_align ();
    virtual void add_parameters_filters ();

    virtual void add_cmdline_seq_dir ();
    virtual void add_cmdline_ktsearch ();
    virtual void add_cmdline_align ();
    virtual void add_cmdline_filters ();

    virtual bool prepareCmdlineFormat ();
    virtual bool prepareParameters ();
    virtual bool interpreteParameters ();

public:
    KTSearch_params (const char* header = NULL, const char* procname = NULL, const char* version = NULL)
    : Search_params (header, procname, version) {}

    bool search_forward () const {return search_forward_;}
    bool search_reverse () const {return search_reverse_;}

    int  rep_lookup () const {return rep_lookup_;}
    bool approx () const {return approx_;}
    int  k_size () const {return k_size_;}
    int  k_thresh () const {return k_thresh_;}
    int  max_shift () const {return max_shift_;}
    int  gap_period () const {return gap_period_;}
    int  min_thresh () const {return min_thresh_;}
    int  max_thresh () const {return max_thresh_;}
    int  min_len () const {return min_len_;}
    unsigned step () const {return step_;}
    const char*  kdistr () const {return kdistr_.c_str ();}
    double gip () const {return gip_;}
    double gep () const {return gep_;}
    double sim_level () const {return sim_level_;}
    double min_score () const {return min_score_;}
    int rep_len () const {return rep_len_;}
    double rep_percent () const {return rep_percent_;}

    void search_forward (bool opt) {search_forward_ = opt;}
    void search_reverse (bool opt) {search_reverse_ = opt;}
    void rep_lookup (int opt) { rep_lookup_ = opt;}
    void approx (bool opt) { approx_ = opt;}
    void k_size (int opt) { k_size_ = opt;}
    void k_thresh (int opt) { k_thresh_ = opt;}
    void max_shift (int opt) { max_shift_ = opt;}
    void gap_period (int opt) { gap_period_ = opt;}
    void min_thresh (int opt) { min_thresh_ = opt;}
    void max_thresh (int opt) { max_thresh_ = opt;}
    void min_len (int opt) { min_len_ = opt;}
    void step (unsigned opt) { step_ = opt;}
    void kdistr (const char* opt) { kdistr_ = opt;}
    void gip (double opt) { gip_ = opt;}
    void gep (double opt) { gep_ = opt;}
    void sim_level (double opt) { sim_level_ = opt;}
    void min_score (double opt) { min_score_ = opt;}
    void rep_len (int opt) { rep_len_ = opt;}
    void rep_percent (double opt) { rep_percent_ = opt;}

    virtual const char* search_fwd_default () const;
    virtual const char* search_rev_default () const;
    virtual const char* rep_lookup_default () const;
    virtual const char* approx_default () const;
    virtual const char* k_size_default () const;
    virtual const char* k_thresh_default () const;
    virtual const char* max_shift_default () const;
    virtual const char* gap_period_default () const;
    virtual const char* min_thresh_default () const;
    virtual const char* max_thresh_default () const;
    virtual const char* min_len_default () const;
    virtual const char* step_default () const;
    virtual const char* kdistr_default () const;
    virtual const char* gip_default () const;
    virtual const char* gep_default () const;
    virtual const char* sim_level_default () const;
    virtual const char* min_score_default () const;
    virtual const char* rep_len_default () const;
    virtual const char* rep_percent_default () const;

    virtual const char* search_fwd_help () const;
    virtual const char* search_rev_help () const;
    virtual const char* rep_lookup_help () const;
    virtual const char* approx_help () const;
    virtual const char* k_size_help () const;
    virtual const char* k_thresh_help () const;
    virtual const char* max_shift_help () const;
    virtual const char* gap_period_help () const;
    virtual const char* min_thresh_help () const;
    virtual const char* max_thresh_help () const;
    virtual const char* min_len_help () const;
    virtual const char* step_help () const;
    virtual const char* kdistr_help () const;
    virtual const char* gip_help () const;
    virtual const char* gep_help () const;
    virtual const char* sim_level_help () const;
    virtual const char* min_score_help () const;
    virtual const char* rep_len_help () const;
    virtual const char* rep_percent_help () const;

};

#ifndef __KTSEARCH_PARAMS_CPP__
extern const char* KTSEARCH_SECTNAME;
extern const char* ALIGN_SECTNAME;
extern const char* FILTERS_SECTNAME;
#endif

#endif
