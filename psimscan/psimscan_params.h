
//////////////////////////////////////////////////////////////////////////////
// This software module is developed by SCIDM team in 1999-2012.
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
// For any questions please contact SciDM team by email at team@scidm.org
//////////////////////////////////////////////////////////////////////////////

#ifndef __psimscan_params_h__
#define __psimscan_params_h__

#include <search_params.h>

class Psimscan_params : public Search_params
{
    int max_queries_number_;
    int max_total_queries_len_;

    double approx_;
    int k_size_;
    double k_thresh_;
    int step_;
    int max_shift_;
    std::string kdistr_name_;
    double extend_band_;
    int widen_band_;
    double dist_fact_;

    // swsearch
    double gip_;
    double gep_;
    std::string matrix_name_;

    // sim filters
    double min_score_;
    int min_len_;
    bool eval_eval_;


protected:
    bool prepareCmdlineFormat ();
    bool prepareParameters ();
    bool interpreteParameters ();

public:
    Psimscan_params ();

    // getters
    int max_queries_number () const {return max_queries_number_;}
    int max_total_queries_len () const {return max_total_queries_len_;}

    double approx () const {return approx_;}
    int k_size () const {return k_size_;}
    double k_thresh () const {return k_thresh_;}
    int max_shift () const {return max_shift_;}
    int step () const {return step_;}
    const char* kdistr_name () const {return kdistr_name_.c_str ();}
    double extend_band () const { return extend_band_;}
    int widen_band () const { return widen_band_;}
    double dist_fact () const { return dist_fact_;}

    const char* matrix_name () const {return matrix_name_.c_str ();}
    double gip () const {return gip_;}
    double gep () const {return gep_;}

    bool eval_eval () const {return eval_eval_;}
    double min_score () const {return min_score_;}
    int min_len () const {return min_len_;}

    // setters
    void max_queries_number (const int opt) {max_queries_number_ = opt;}
    void max_total_queries_len (const int opt) { max_total_queries_len_ = opt;}

    void approx (double opt) { approx_ = opt;}
    void k_size (int opt) { k_size_ = opt;}
    void k_thresh (double opt) { k_thresh_ = opt;}
    void max_shift (int opt) { max_shift_ = opt;}
    void step (int opt) { step_ = opt;}
    void kdistr_name (const char* opt) { kdistr_name_ = opt;}
    void extend_band (double opt) { extend_band_ = opt;}
    void widen_band (int opt) { widen_band_ = opt;}
    void dist_fact (double opt) { dist_fact_ = opt;}

    void matrix_name (const char* opt) {matrix_name_ = opt;}
    void gip (double opt) {gip_ = opt;}
    void gep (double opt) {gep_ = opt;}

    void eval_eval (bool opt) { eval_eval_ = opt;}
    void min_score (double opt) { min_score_ = opt;}
    void min_len (int opt) { min_len_ = opt;}

    // default providers
    virtual const char* max_queries_number_default () const;
    virtual const char* max_total_queries_len_default () const;

    virtual const char* approx_default () const;
    virtual const char* k_size_default () const;
    virtual const char* k_thresh_default () const;
    virtual const char* max_shift_default () const;
    virtual const char* step_default () const;
    virtual const char* kdistr_name_default () const;
    virtual const char* extend_band_default () const;
    virtual const char* widen_band_default () const;
    virtual const char* dist_fact_default () const;

    virtual const char* matrix_name_default () const;
    virtual const char* gip_default () const;
    virtual const char* gep_default () const;

    virtual const char* eval_eval_default () const;
    virtual const char* min_score_default () const;
    virtual const char* min_len_default () const;

    virtual const char* max_seq_len_default () const; // overwrites Search_params one
    virtual const char* min_seq_len_default () const; // overwrites Search_params one
    virtual const char* max_qry_len_default () const; // overwrites Search_params one

    // helpers

    virtual const char* max_queries_number_help () const;
    virtual const char* max_total_queries_len_help () const;

    virtual const char* approx_help () const;
    virtual const char* k_size_help () const;
    virtual const char* k_thresh_help () const;
    virtual const char* max_shift_help () const;
    virtual const char* step_help () const;
    virtual const char* kdistr_name_help () const;
    virtual const char* extend_band_help () const;
    virtual const char* widen_band_help () const;
    virtual const char* dist_fact_help () const;

    virtual const char* matrix_name_help () const;
    virtual const char* gip_help () const;
    virtual const char* gep_help () const;

    virtual const char* eval_eval_help () const;
    virtual const char* min_score_help () const;
    virtual const char* min_len_help () const;

    virtual const char* query_name_help () const;
    virtual const char* search_name_help () const;
    virtual const char* output_name_help () const;
};

#endif // __psimscan_params_h__
