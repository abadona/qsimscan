
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

#include "search_helper_files.h"
#include "print_batches.h"
#include "parameters.h"
#include "fasta.h"
#include "resource.h"
#include <unistd.h>
#include <iomanip>

Search_helper_files::Search_helper_files ()
:
iter_init_ (false),
iter_done_ (false),
skipped_ (0)
{
}

int* Search_helper_files::load_ktuple_distribution (const char* fname, double sim_level, int& mismatch_score)
{
    FILE* q = fopen (fname, "r");
    if (!q)
        ers << "Tuple distribution File " << fname << "does not exist." << Throw;

    longlong distr_size = fseek (q, 0, SEEK_END) / 4;
    fseek (q, 0, SEEK_SET);

    // calculate distr size; check consistency
    int ksize = 0;
    int exp = 1;
    while (exp < distr_size)
    {
        ksize ++;
        exp *= 4;
    }
    if (exp != distr_size)
        ers << "Number of entries in tuple distribution file (" << distr_size <<") is not a power of 4" << Throw;

    MemWrapper <int> result  (distr_size);
    if (!result) ERR (NOEMEM);

    fread (result, 4, distr_size, q);
    fclose (q);

    // calculate the average ktuple score and mismatch score for gip/gep evaluation
    longlong avg = 0;
    for (int ckt = 0; ckt < distr_size; ckt ++)
        avg += result [ckt];
    avg /= distr_size;
    mismatch_score = (int) ((avg / ksize) * sim_level) / 100;

    return result.release ();
}

int* Search_helper_files::create_ktuple_distribution (int ksize, double sim_level, int& mismatch_score)
{
    int  distr_size = 1; for (int pwr = 0; pwr < ksize; pwr ++) distr_size *= 4;
    MemWrapper <int> scores (distr_size);
    for (int ii = 0; ii < distr_size; ii ++)
        scores [ii] = 100;

    // calculate the average ktuple score and mismatch score for gip/gep evaluation
    mismatch_score = (int) ((100 / ksize) * sim_level) / 100;

    return scores.release ();

}

unsigned Search_helper_files::load_queries_nn (const char* fname, std::vector <NN_SEQ>& dest, unsigned min_len, unsigned max_len, unsigned* total, unsigned beg, unsigned end)
{
    if (!query_.open (fname))
        ers << "Unable to open file " << fname << " for reading" << Throw;

    // load sequence data
    unsigned seqno = 0;
    if (total) *total = 0;
    while (query_.next ())
    {
        if (seqno >= end) break;
        if (seqno >= beg)
        {
            if (query_.cur_seq_len () >= min_len && query_.cur_seq_len () <= max_len)
            {
                NN_SEQ n;
                n.len = query_.cur_seq_len ();
                n.seq = new char [(n.len + 3) >> 2];
                n_ascii2binary (n.seq, int ((n.len + 3) >> 2), query_.cur_seq (), 0, n.len);
                n.rev = 0;
                n.uid = query_.cur_recstart ();
                n.owner = true;
                dest.push_back (n);
            }
            if (total) (*total)++;
        }
        seqno ++;
    }
    return dest.size ();
}
unsigned Search_helper_files::load_queries_aa (const char* fname, std::vector <AA_SEQ>& dest, unsigned min_len, unsigned max_len, unsigned* total, unsigned beg, unsigned end)
{
    if (!query_.open (fname))
        ers << "Unable to open file " << fname << " for reading" << Throw;

    // load sequence data
    unsigned seqno = 0;
    if (total) *total = 0;
    while (query_.next ())
    {
        if (seqno >= end) break;
        if (seqno >= beg)
        {
            if (query_.cur_seq_len () >= min_len && query_.cur_seq_len () <= max_len)
            {
                AA_SEQ n;
                n.len = query_.cur_seq_len ();
                n.seq = new char [n.len];
                a_ascii2binary (n.seq, n.len, query_.cur_seq (), 0, n.len);
                n.rev = 0;
                n.uid = query_.cur_recstart ();
                n.owner = true;
                dest.push_back (n);
            }
            if (total) (*total)++;
        }
        seqno ++;
    }
    return dest.size ();
}

unsigned Search_helper_files::load_queries_na (const char* fname, std::vector <NA_SEQ>& dest, unsigned min_len, unsigned max_len, unsigned* total, unsigned beg, unsigned end)
{
    if (!query_.open (fname))
        ers << "Unable to open file " << fname << " for reading" << Throw;

    // load sequence data
    unsigned seqno = 0;
    if (total) *total = 0;
    while (query_.next ())
    {
        if (seqno >= end) break;
        if (seqno >= beg)
        {
            if (query_.cur_seq_len () >= min_len && query_.cur_seq_len () <= max_len)
            {
                NA_SEQ n;
                n.len = query_.cur_seq_len ();
                n.seq = new char [(n.len + 3) >> 2];
                n_ascii2binary (n.seq, int ((n.len + 3) >> 2), query_.cur_seq (), 0, n.len);
                n.rev = 0;
                n.uid = query_.cur_recstart ();
                n.owner = true;
                dest.push_back (n);
            }
            if (total) (*total)++;
        }
        seqno ++;
    }
    return dest.size ();
}

void Search_helper_files::init_nn_search (const char* fname, unsigned beg, unsigned end)
{
    skipped_ = 0;
    beg_ = beg, end_ = end;
    target_file_name_ = fname;
    if (!target_.open (fname))
        ers << "Unable to open file " << fname << ThrowEx (OSRerror);
    iter_init_ = true;
}

NN_SEQ*      Search_helper_files::next_nn_seq (unsigned min_len, unsigned max_len)
{
    if (!iter_init_ || iter_done_)
        return false;
    while (target_.next ())
    {
        if (target_no_ >= end_) break;
        if (target_no_ >= beg_)
        {
            if (target_.cur_seq_len () >= min_len && target_.cur_seq_len () <= max_len)
            {
                n_xseq_.len = target_.cur_seq_len ();
                n_xseq_.seq = new char [(n_xseq_.len + 3) >> 2];
                n_ascii2binary (n_xseq_.seq, int ((n_xseq_.len + 3) >> 2), target_.cur_seq (), 0, n_xseq_.len);
                n_xseq_.rev = 0;
                n_xseq_.uid = target_.cur_recstart ();
                target_no_ ++;
                return &n_xseq_;
            }
            else
                skipped_ ++;
        }
        target_no_ ++;
    }
    iter_done_ = true;
    return NULL;
}

void Search_helper_files::init_aa_search (const char* fname, unsigned beg, unsigned end)
{
    init_nn_search (fname, beg, end);
}
AA_SEQ* Search_helper_files::next_aa_seq (unsigned min_len, unsigned max_len)
{
    if (!iter_init_ || iter_done_)
        return false;
    while (target_.next ())
    {
        if (target_no_ >= end_) break;
        if (target_no_ >= beg_)
        {
            if (target_.cur_seq_len () >= min_len && target_.cur_seq_len () <= max_len)
            {
                p_xseq_.len = target_.cur_seq_len ();
                p_xseq_.seq = new char [p_xseq_.len];
                a_ascii2binary (p_xseq_.seq, p_xseq_.len, target_.cur_seq (), 0, p_xseq_.len);
                p_xseq_.rev = 0;
                p_xseq_.uid = target_.cur_recstart ();
                target_no_ ++;
                return &p_xseq_;
            }
            else
                skipped_ ++;
        }
        target_no_ ++;
    }
    iter_done_ = true;
    return NULL;
}

void Search_helper_files::init_na_search (const char* fname, unsigned beg, unsigned end)
{
    init_nn_search (fname, beg, end);
}

NA_SEQ* Search_helper_files::next_na_seq (unsigned min_len, unsigned max_len)
{
    if (!iter_init_ || iter_done_)
        return false;
    while (target_.next ())
    {
        if (target_no_ >= end_) break;
        if (target_no_ >= beg_)
        {
            if (target_.cur_seq_len () >= min_len && target_.cur_seq_len () <= max_len)
            {
                na_xseq_.len = target_.cur_seq_len ();
                na_xseq_.seq = new char [(na_xseq_.len + 3) >> 2];
                n_ascii2binary (na_xseq_.seq, int ((na_xseq_.len + 3) >> 2), target_.cur_seq (), 0, na_xseq_.len);
                na_xseq_.rev = 0;
                na_xseq_.uid = target_.cur_recstart ();
                target_no_ ++;
                return &na_xseq_;
            }
            else
                skipped_ ++;
        }
        target_no_ ++;
    }
    iter_done_ = true;
    return NULL;
}


bool Search_helper_files::output_results_nn (AlignResultStorage& resrec, unsigned res_per_query, unsigned alignments_no, std::vector<NN_SEQ>& f_qry, std::vector<NN_SEQ>& r_qry, WMatrix* w, std::ostream& out_file, unsigned max_x_len)
{
    MemWrapper <char> bin_buf ((max_x_len + 3) >> 2);
    out_file << std::endl;
    for (unsigned query_no = 0; query_no < f_qry.size (); query_no ++)
    {
        NN_SEQ& fwd_qry = f_qry [query_no];
        unsigned sim_found = resrec.resPerQuery (fwd_qry.uid);

        query_.seek (fwd_qry.uid);
        if (!query_.next ())
            ers << "Unable to read query sequence at offset " << fwd_qry.uid << Throw;


        out_file << "[Query " << query_no / (rev_?1:2) + 1 << "] " << sim_found << " similarities found." << std::endl;
        out_file << query_.cur_name () << " (" << fwd_qry.len << " bp), " << query_.cur_hdr ();
        out_file << std::endl;

        if (sim_found)
        {
            // get results order for current seq
            unsigned resno = 0;
            // make results order
            ResultQueue::ElemVec& query_results = *resrec.getQueryResults (fwd_qry.uid);
            for (ResultQueue::ElemVec::iterator cur_res = query_results.begin (); cur_res != query_results.end (); cur_res ++, resno ++)
            {
                NN_SEQ cur_search;
                cur_search.uid = cur_res->uid_;
                target_.seek (cur_search.uid);
                if (!target_.next ())
                    ers << "Unable to read target sequence at offset " << cur_search.uid << Throw;
                cur_search.len = target_.cur_seq_len ();

                unsigned tot_sim_len = 0;
                for (unsigned batch_idx = 0; batch_idx < cur_res->batch_no_; batch_idx ++)
                    tot_sim_len += cur_res->batches_[batch_idx].len;

                // print batch parameters
                out_file << std::endl << "     [Result " << resno+1 << "]";
                out_file << std::endl << "     " << target_.cur_name () << " (" << cur_search.len << " bp), " << target_.cur_hdr ();
                out_file << std::endl << "     " << ((cur_res->reverse_)?("REVERSE") : ("FORWARD")) << ", SCORE " << cur_res->al_score_;
                if (cur_res->evalue_ != 0)   out_file << ", EVALUE " << cur_res->evalue_;
                if (cur_res->chi2_ != 0)     out_file << ", CHI-2 " << cur_res->chi2_;
                if (cur_res->bitscore_ != 0) out_file << ", BITSCORE " << cur_res->bitscore_;
                if (cur_res->score_ != 0)    out_file << ", MATCH_SCORE " << cur_res->score_;
                if (cur_res->q_auto_score_ || cur_res->t_auto_score_)
                {
                    out_file << std::endl << "     ";
                    if (cur_res->q_auto_score_) out_file << "QUERY AUTO SCORE " << cur_res->q_auto_score_;
                    if (cur_res->q_auto_score_ && cur_res->t_auto_score_) out_file << ", ";
                    if (cur_res->q_auto_score_) out_file << "TARGET AUTO SCORE " << cur_res->t_auto_score_;
                }

                out_file << std::endl << "     QUERY START " << cur_res->batches_ [0].xpos << ", " << query_.cur_name () << " START " << cur_res->batches_ [0].ypos << ", LENGTH " << tot_sim_len;

                // print batch
                if (resno < alignments_no)
                {
                    n_ascii2binary (bin_buf, ((max_x_len + 3) >> 2), target_.cur_seq (), 0, cur_search.len);
                    cur_search.seq = bin_buf;
                    print_batches (((cur_res->reverse_) ? r_qry [query_no] : fwd_qry), cur_search, cur_res->batches_, cur_res->batch_no_, w, out_file);
                }
                out_file << std::endl;
            }
        }
    }
    return true;
}

bool Search_helper_files::output_results_aa (AlignResultStorage& resrec, unsigned res_per_query, unsigned alignments_no, std::vector<AA_SEQ>& qry, WMatrix* w, std::ostream& out_file, unsigned max_x_len)
{
    MemWrapper <char> bin_buf (max_x_len);
    out_file << std::endl;

    for (unsigned query_no = 0; query_no < qry.size (); query_no ++)
    {
        AA_SEQ& q = qry [query_no];

        unsigned sim_found = resrec.resPerQuery (q.uid);

        query_.seek (q.uid);
        if (!query_.next ())
            ers << "Unable to read query sequence at offset " << q.uid << Throw;


        out_file << "[Query " << query_no + 1 << "] " << sim_found << " similarities found." << std::endl;
        out_file << query_.cur_name () << " (" << q.len << " bp), " << query_.cur_hdr ();
        out_file << std::endl;

        if (sim_found)
        {
            // get results order for current seq
            unsigned resno = 0;
            // make results order
            ResultQueue::ElemVec& query_results = *resrec.getQueryResults (q.uid);
            for (ResultQueue::ElemVec::iterator cur_res = query_results.begin (); cur_res != query_results.end (); cur_res ++, resno ++)
            {
                AA_SEQ cur_search;
                cur_search.uid = cur_res->uid_;
                target_.seek (cur_search.uid);
                if (!target_.next ())
                    ers << "Unable to read target sequence at offset " << cur_search.uid << Throw;
                cur_search.len = target_.cur_seq_len ();

                unsigned tot_sim_len = 0;
                for (unsigned batch_idx = 0; batch_idx < cur_res->batch_no_; batch_idx ++)
                    tot_sim_len += cur_res->batches_[batch_idx].len;

                // print batch parameters
                out_file << std::endl << "     [Result " << resno+1 << "]";
                out_file << std::endl << "     " << target_.cur_name () << " (" << cur_search.len << " bp), " << target_.cur_hdr ();
                out_file << std::endl << "     SCORE " << cur_res->al_score_;
                if (cur_res->evalue_ != 0)   out_file << ", EVALUE " << cur_res->evalue_;
                if (cur_res->chi2_ != 0)     out_file << ", CHI-2 " << cur_res->chi2_;
                if (cur_res->bitscore_ != 0) out_file << ", BITSCORE " << cur_res->bitscore_;
                if (cur_res->score_ != 0)    out_file << ", MATCH_SCORE " << cur_res->score_;
                if (cur_res->q_auto_score_ || cur_res->t_auto_score_)
                {
                    out_file << std::endl << "     ";
                    if (cur_res->q_auto_score_) out_file << "QUERY AUTO SCORE " << cur_res->q_auto_score_;
                    if (cur_res->q_auto_score_ && cur_res->t_auto_score_) out_file << ", ";
                    if (cur_res->q_auto_score_) out_file << "TARGET AUTO SCORE " << cur_res->t_auto_score_;
                }
                out_file << std::endl << "     QUERY START " << cur_res->batches_ [0].xpos << ", " << query_.cur_name () << " START " << cur_res->batches_ [0].ypos << ", LENGTH " << tot_sim_len;

                // print batch
                if (resno < alignments_no)
                {
                    a_ascii2binary (bin_buf, max_x_len, target_.cur_seq (), 0, cur_search.len);
                    cur_search.seq = bin_buf;
                    print_batches (q, cur_search, cur_res->batches_, cur_res->batch_no_, w, out_file);
                }
                out_file << std::endl;
            }
        }
    }
    return true;
}

bool Search_helper_files::output_results_an (AlignResultStorage& resrec, unsigned res_per_query, unsigned alignments_no, std::vector<AA_SEQ>& qry, WMatrix* w, std::ostream& out_file, unsigned max_x_len)
{
    MemWrapper <char> bin_buf ((max_x_len+3)>>2);
    MemWrapper <char> rev_buf ((max_x_len+3)>>2);
    out_file << std::endl;

    for (unsigned query_no = 0; query_no < qry.size (); query_no ++)
    {
        AA_SEQ& q = qry [query_no];
        query_.seek (q.uid);
        if (!query_.next ())
            ers << "Unable to read query sequence at offset " << q.uid << Throw;

        unsigned sim_found = resrec.resPerQuery (q.uid);

        out_file << "[Query " << query_no + 1 << "] " << sim_found << " similarities found." << std::endl;
        out_file << query_.cur_name () << " (" << q.len << " bp), " << query_.cur_hdr ();
        out_file << std::endl << std::endl;

        if (sim_found)
        {
            // get results order for current seq
            unsigned resno = 0;
            // make results order
            ResultQueue::ElemVec& query_results = *resrec.getQueryResults (q.uid);
            for (ResultQueue::ElemVec::iterator cur_res = query_results.begin (); cur_res != query_results.end (); cur_res ++, resno ++)
            {
                NA_SEQ cur_search;

                cur_search.uid = cur_res->uid_;
                target_.seek (cur_search.uid);
                if (!target_.next ())
                    ers << "Unable to read target sequence at offset " << cur_search.uid << Throw;
                cur_search.len = target_.cur_seq_len ();

                unsigned tot_sim_len = 0;
                for (unsigned batch_idx = 0; batch_idx < cur_res->batch_no_; batch_idx ++)
                    tot_sim_len += cur_res->batches_ [batch_idx].len;

                // print batch parameters
                out_file << std::endl << "     [Result " << resno+1 << "]";
                out_file << std::endl << "     " << target_.cur_name () << " (" << cur_search.len << " bp), " << target_.cur_hdr ();
                if (cur_res->reverse_) out_file << std::endl << "     Reverse chain ";
                else                   out_file << std::endl << "     Forward chain ";
                out_file << std::endl << "     SCORE " << cur_res->al_score_;
                if (cur_res->evalue_ != 0)   out_file << ", EVALUE " << cur_res->evalue_;
                if (cur_res->chi2_ != 0)     out_file << ", CHI-2 " << cur_res->chi2_;
                if (cur_res->bitscore_ != 0) out_file << ", BITSCORE " << cur_res->bitscore_;
                if (cur_res->score_ != 0)    out_file << ", MATCH_SCORE " << cur_res->score_;
                if (cur_res->q_auto_score_ || cur_res->t_auto_score_)
                {
                    out_file << std::endl << "     ";
                    if (cur_res->q_auto_score_) out_file << "QUERY AUTO SCORE " << cur_res->q_auto_score_;
                    if (cur_res->q_auto_score_ && cur_res->t_auto_score_) out_file << ", ";
                    if (cur_res->q_auto_score_) out_file << "TARGET AUTO SCORE " << cur_res->t_auto_score_;
                }
                out_file << std::endl << "     QUERY START " << cur_res->batches_ [0].xpos << ", " << query_.cur_name () << " START " << cur_res->batches_ [0].ypos << ", LENGTH " << tot_sim_len;

                // print batch
                if (resno < alignments_no)
                {
                    n_ascii2binary (bin_buf, ((max_x_len + 3) >> 2), target_.cur_seq (), 0, cur_search.len);
                    cur_search.seq = bin_buf;
                    if (cur_res->reverse_)
                    {
                        n_revert_seq (rev_buf, cur_search.seq, cur_search.len);
                        cur_search.seq = rev_buf;
                    }
                    print_batches (q, cur_search, cur_res->batches_, cur_res->batch_no_, w, out_file);
                }
                out_file << std::endl;
            }
        }
    }
    return true;
}

bool Search_helper_files::output_results_na (AlignResultStorage& resrec, unsigned res_per_query, unsigned alignments_no, std::vector<NA_SEQ>& f_qry, std::vector<NA_SEQ>& r_qry, WMatrix* w, std::ostream& out_file, unsigned max_x_len)
{
    MemWrapper <char> bin_buf (max_x_len);
    out_file << std::endl;

    for (unsigned query_no = 0; query_no < f_qry.size (); query_no ++)
    {
        NA_SEQ& fwd_qry = f_qry [query_no];
        NA_SEQ* rev_qry_p = r_qry.size () ? &(r_qry [query_no]) : NULL;
        // sanity check
        if (rev_qry_p && rev_qry_p->uid != fwd_qry.uid) ERR(ERR_Internal);

        unsigned sim_found = resrec.resPerQuery (fwd_qry.uid);

        query_.seek (fwd_qry.uid);
        if (!query_.next ())
            ers << "Unable to read query sequence at offset " << fwd_qry.uid << Throw;


        out_file << "[Query " << query_no + 1 << "] " << sim_found << " similarities found." << std::endl;
        out_file << query_.cur_name () << " (" << fwd_qry.len << " bp), " << query_.cur_hdr ();
        out_file << std::endl;

        if (sim_found)
        {
            // get results order for current seq
            unsigned resno = 0;
            // make results order
            ResultQueue::ElemVec& query_results = *resrec.getQueryResults (fwd_qry.uid);
            for (ResultQueue::ElemVec::iterator cur_res = query_results.begin (); cur_res != query_results.end (); cur_res ++, resno ++)
            {
                AA_SEQ cur_search;
                cur_search.uid = cur_res->uid_;
                target_.seek (cur_search.uid);
                if (!target_.next ())
                    ers << "Unable to read target sequence at offset " << cur_search.uid << Throw;
                cur_search.len = target_.cur_seq_len ();

                unsigned tot_sim_len = 0;
                for (unsigned batch_idx = 0; batch_idx < cur_res->batch_no_; batch_idx ++)
                    tot_sim_len += cur_res->batches_ [batch_idx].len;

                // print batch parameters
                out_file << std::endl << "     [Result " << resno+1 << "]";
                out_file << std::endl << "     " << target_.cur_name () << "   " << cur_search.len << "aa, " << target_.cur_hdr ();
                if (!cur_res->reverse_) out_file << std::endl << "     Matches REVERSE chain of query";
                else                    out_file << std::endl << "     Matches FORWARD chain of query";
                out_file << std::endl << "     SCORE " << cur_res->al_score_;
                if (cur_res->evalue_ != 0)   out_file << ", EVALUE " << cur_res->evalue_;
                if (cur_res->chi2_ != 0)     out_file << ", CHI-2 " << cur_res->chi2_;
                if (cur_res->bitscore_ != 0) out_file << ", BITSCORE " << cur_res->bitscore_;
                if (cur_res->score_ != 0)    out_file << ", MATCH_SCORE " << cur_res->score_;
                if (cur_res->q_auto_score_ || cur_res->t_auto_score_)
                {
                    out_file << std::endl << "     ";
                    if (cur_res->q_auto_score_) out_file << "QUERY AUTO SCORE " << cur_res->q_auto_score_;
                    if (cur_res->q_auto_score_ && cur_res->t_auto_score_) out_file << ", ";
                    if (cur_res->q_auto_score_) out_file << "TARGET AUTO SCORE " << cur_res->t_auto_score_;
                }

                out_file << std::endl << "     QUERY START " << cur_res->batches_ [0].xpos << ", " << query_.cur_name () << " START " << cur_res->batches_ [0].ypos << ", LENGTH " << tot_sim_len;

                // print batch
                if (resno < alignments_no)
                {
                    a_ascii2binary (bin_buf, max_x_len, target_.cur_seq (), 0, cur_search.len);
                    cur_search.seq = bin_buf;
                    if (cur_res->reverse_ && !rev_qry_p) ERR (ERR_Internal);
                    print_batches (((cur_res->reverse_) ? *rev_qry_p : fwd_qry), cur_search, cur_res->batches_, cur_res->batch_no_, w, out_file);
                }
                out_file << std::endl;
            }
        }
    }
    return true;
}


void eval_align (SEQ& xseq, SEQ& yseq, WMatrix* w, BATCH* batches, unsigned batch_no, double* p_identity, unsigned* mismatches, unsigned* al_length, unsigned* gap_openings, unsigned* gap_length)
{

    *p_identity = 0, *mismatches = 0, *al_length = 0, *gap_openings = 0;
    unsigned xp, yp;

    for (unsigned i = 0; i < batch_no; i ++)
    {
        if (i)
        {
            (*gap_openings) ++;
            (*gap_length) += (batches [i].xpos - xp) + (batches [i].ypos - yp);
        }

        xp = batches [i].xpos;
        yp = batches [i].ypos;
        for (unsigned p = 0; p < batches [i].len; p ++)
        {
            if (xseq.get_code (xp) != yseq.get_code (yp))
                (*mismatches) ++;
            xp ++;
            yp ++;
        }
        (*al_length) += batches [i].len;
    }
    *p_identity = (double (*al_length) - double (*mismatches)) * 100.0 / double (*al_length);
}

static char const* TAB_STR = "\t";

bool Search_helper_files::output_results_m8_aa (AlignResultStorage& resrec, unsigned res_per_query, unsigned alignments_no, std::vector<AA_SEQ>& qry, WMatrix* w, std::ostream& o, const char* hdr)
{

    const char* Query_id = "";
    const char* Subject_id = "";
    double p_identity = 0;
    unsigned alignment_length = 0;
    unsigned mismatches = 0;
    unsigned gap_openings = 0;
    unsigned gap_length = 0;
    unsigned q_start = 0;
    unsigned q_end = 0;
    unsigned s_start = 0;
    unsigned s_end = 0;
    double e_value = 0;
    double bit_score = 0;

    for (unsigned query_no = 0; query_no < qry.size (); query_no ++)
    {
        AA_SEQ& cur_query = qry [query_no];

        unsigned sim_found = resrec.resPerQuery (cur_query.uid);

        if (sim_found)
        {
            query_.seek (cur_query.uid);
            if (!query_.next ())
                ers << "Unable to read query sequence at offset " << cur_query.uid << Throw;
            Query_id = query_.cur_name ();

            if (hdr)
            {
                // # program
                o << "# " << hdr << std::endl;
                // # Query
                o << "# Query: " << Query_id << " " << query_.cur_hdr () << std::endl;
                // # Database
                o << "# Database: " << target_file_name_ << std::endl;
                // # Fields
                o << "# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score" << std::endl;
            }

            // get results order for current sequence
            unsigned resno = 0;
            // make results order
            ResultQueue::ElemVec& query_results = *resrec.getQueryResults (cur_query.uid);
            for (ResultQueue::ElemVec::iterator cur_res = query_results.begin (); cur_res != query_results.end (); cur_res ++, resno ++)
            {
                AA_SEQ cur_search;
                cur_search.uid = cur_res->uid_;

                target_.seek (cur_search.uid);
                if (!target_.next ())
                    ers << "Unable to read target sequence at offset " << cur_search.uid << Throw;
                cur_search.len = strlen (target_.cur_seq ());
                cur_search.rev = 0;
                cur_search.seq = const_cast <char*> (target_.cur_seq ());
                // inplace conversion
                a_ascii2binary (cur_search.seq, cur_search.len, cur_search.seq, 0, cur_search.len);

                // Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
                Subject_id = target_.cur_name ();

                BATCH* batches = cur_res->batches_;
                unsigned batch_no = cur_res->batch_no_;

                eval_align (cur_query, cur_search, w, batches, batch_no, &p_identity, &mismatches, &alignment_length, &gap_openings, &gap_length);

                // here query IS x!
                q_start = batches [0].xpos + 1;
                q_end = batches [batch_no-1].xpos + batches [batch_no-1].len;
                s_start = batches [0].ypos + 1;
                s_end = batches [batch_no-1].ypos + batches [batch_no-1].len;
                e_value = cur_res->evalue_;
                bit_score = cur_res->bitscore_;

                o <<  Query_id << TAB_STR <<
                      Subject_id << TAB_STR <<
                      std::setprecision (2) << std::fixed << std::noshowpoint << p_identity << TAB_STR <<
                      alignment_length << TAB_STR <<
                      mismatches << TAB_STR <<
                      gap_openings << TAB_STR <<
                      q_start+1 << TAB_STR <<
                      q_end << TAB_STR <<
                      s_start+1 << TAB_STR <<
                      s_end << TAB_STR;
                o.unsetf (std::ios::fixed | std::ios::scientific);
                o <<  std::noshowpoint << std::setprecision (2) << e_value << TAB_STR <<
                      std::setprecision (2) << std::fixed << std::noshowpoint << bit_score << std::endl;
            }
        }
    }
    return true;
}

bool Search_helper_files::output_results_m8_nn (AlignResultStorage& resrec, unsigned res_per_query, unsigned alignments_no, std::vector<NN_SEQ>& f_qry, std::vector<NN_SEQ>& r_qry, WMatrix* w, std::ostream& o, const char* hdr)
{
    const char* Query_id = "";
    const char* Subject_id = "";
    double p_identity = 0;
    unsigned alignment_length = 0;
    unsigned mismatches = 0;
    unsigned gap_openings = 0;
    unsigned gap_length = 0;
    unsigned q_start = 0;
    unsigned q_end = 0;
    unsigned s_start = 0;
    unsigned s_end = 0;
    double e_value = 0;
    double bit_score = 0;

    for (unsigned query_no = 0; query_no < f_qry.size (); query_no ++)
    {
        NN_SEQ& fwd_qry = f_qry [query_no];
        NN_SEQ* rev_qry_p = r_qry.size () ? &(r_qry [query_no]) : NULL;
        // sanity check
        if (rev_qry_p && rev_qry_p->uid != fwd_qry.uid) ERR(ERR_Internal);

        unsigned sim_found = resrec.resPerQuery (fwd_qry.uid);

        if (sim_found)
        {
            query_.seek (fwd_qry.uid);
            if (!query_.next ())
                ers << "Unable to read query sequence at offset " << fwd_qry.uid << Throw;
            Query_id = query_.cur_name ();

            if (hdr)
            {
                // # program
                o << "# " << hdr << std::endl;
                // # Query
                o << "# Query: " << Query_id << " " << query_.cur_hdr () << std::endl;
                // # Database
                o << "# Database: " << target_file_name_ << std::endl;
                // # Fields
                o << "# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score" << std::endl;
            }

            // get results order for current sequence
            unsigned resno = 0;
            // make results order
            ResultQueue::ElemVec& query_results = *resrec.getQueryResults (fwd_qry.uid);
            for (ResultQueue::ElemVec::iterator cur_res = query_results.begin (); cur_res != query_results.end (); cur_res ++, resno ++)
            {
                NN_SEQ cur_search;
                cur_search.uid = cur_res->uid_;

                target_.seek (cur_search.uid);
                if (!target_.next ())
                    ers << "Unable to read target sequence at offset " << cur_search.uid << Throw;
                cur_search.len = strlen (target_.cur_seq ());
                cur_search.rev = 0;
                cur_search.seq = const_cast <char*> (target_.cur_seq ());
                // inplace conversion
                n_ascii2binary (cur_search.seq, cur_search.len, cur_search.seq, 0, cur_search.len);

                // Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
                Subject_id = target_.cur_name ();

                BATCH* batches = cur_res->batches_;
                unsigned batch_no = cur_res->batch_no_;

                if (cur_res->reverse_ && !rev_qry_p) ERR (ERR_Internal);
                eval_align (cur_res->reverse_ ? *rev_qry_p : fwd_qry, cur_search, w, batches, batch_no, &p_identity, &mismatches, &alignment_length, &gap_openings, &gap_length);

                // here query IS x!
                q_start = cur_res->reverse_ ? fwd_qry.len - (batches [batch_no-1].xpos + batches [batch_no-1].len) : batches [0].xpos + 1;
                q_end = cur_res->reverse_ ? fwd_qry.len - batches [0].xpos : batches [batch_no-1].xpos + batches [batch_no-1].len;
                s_start = batches [0].ypos + 1;
                s_end = batches [batch_no-1].ypos + batches [batch_no-1].len;
                e_value = cur_res->chi2_;
                bit_score = cur_res->al_score_;

                o <<  Query_id << TAB_STR <<
                      Subject_id << TAB_STR <<
                      std::setprecision (2) << std::fixed << std::noshowpoint << p_identity << TAB_STR <<
                      alignment_length << TAB_STR <<
                      mismatches << TAB_STR <<
                      gap_openings << TAB_STR <<
                      q_start+1 << TAB_STR <<
                      q_end << TAB_STR <<
                      s_start+1 << TAB_STR <<
                      s_end << TAB_STR;
                o.unsetf (std::ios::fixed | std::ios::scientific);
                o <<  std::setprecision (2) << std::noshowpoint << e_value << TAB_STR <<
                      std::setprecision (2) << std::fixed << std::noshowpoint << bit_score << std::endl;
            }
        }
    }
    return true;
}

bool Search_helper_files::output_results_tab_aa (AlignResultStorage& resrec, unsigned res_per_query, unsigned alignments_no, std::vector<AA_SEQ>& qry, WMatrix* w, std::ostream& o, bool hdr)
{

    if (hdr) o << "#Q_id\tS_id\tp_inden\tal_len\tmism\tgaps\tgap_len\tqry_beg\tqry_end\tqry_len\ttrg_beg\ttrg_end\ttrg_len\tevalue\tsw_score\tqry_auto\ttrg_auto" << std::endl;

    const char* Query_id = "";
    const char* Subject_id = "";
    double p_identity = 0;
    unsigned alignment_length = 0;
    unsigned mismatches = 0;
    unsigned gap_openings = 0;
    unsigned gap_length = 0;
    unsigned q_start = 0;
    unsigned q_end = 0;
    unsigned s_start = 0;
    unsigned s_end = 0;
    double e_value = 0;
    double sw_score = 0;

    for (unsigned query_no = 0; query_no < qry.size (); query_no ++)
    {
        AA_SEQ& cur_query = qry [query_no];

        unsigned sim_found = resrec.resPerQuery (cur_query.uid);

        if (sim_found)
        {
            query_.seek (cur_query.uid);
            if (!query_.next ())
                ers << "Unable to read query sequence at offset " << cur_query.uid << Throw;
            Query_id = query_.cur_name ();

            // get results order for current sequence
            unsigned resno = 0;
            // make results order
            ResultQueue::ElemVec& query_results = *resrec.getQueryResults (cur_query.uid);
            for (ResultQueue::ElemVec::iterator cur_res = query_results.begin (); cur_res != query_results.end (); cur_res ++, resno ++)
            {
                AA_SEQ cur_search;
                cur_search.uid = cur_res->uid_;

                target_.seek (cur_search.uid);
                if (!target_.next ())
                    ers << "Unable to read target sequence at offset " << cur_search.uid << Throw;
                cur_search.len = strlen (target_.cur_seq ());
                cur_search.rev = 0;
                cur_search.seq = const_cast <char*> (target_.cur_seq ());
                // inplace conversion
                a_ascii2binary (cur_search.seq, cur_search.len, cur_search.seq, 0, cur_search.len);

                // Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
                Subject_id = target_.cur_name ();

                BATCH* batches = cur_res->batches_;
                unsigned batch_no = cur_res->batch_no_;

                eval_align (cur_query, cur_search, w, batches, batch_no, &p_identity, &mismatches, &alignment_length, &gap_openings, &gap_length);

                // here query IS x!
                q_start = batches [0].xpos;
                q_end = batches [batch_no-1].xpos + batches [batch_no-1].len;
                s_start = batches [0].ypos;
                s_end = batches [batch_no-1].ypos + batches [batch_no-1].len;
                e_value = cur_res->evalue_;
                sw_score = cur_res->al_score_;

                // QueryId SubjectId P_indentity AlignmentLength Mismatches GapOpenings GapLength QueryStart QueryEnd QueryLen TargetStart TargetEnd TargetLen EValue SW_score AutoQuery AutoTarget

                o <<
                    Query_id << TAB_STR <<
                    Subject_id << TAB_STR <<
                    std::setprecision (2) << std::fixed << std::noshowpoint << p_identity << TAB_STR <<
                    alignment_length << TAB_STR <<
                    mismatches << TAB_STR <<
                    gap_openings << TAB_STR <<
                    gap_length << TAB_STR <<
                    q_start << TAB_STR <<
                    q_end << TAB_STR <<
                    cur_query.len << TAB_STR <<
                    s_start << TAB_STR <<
                    s_end << TAB_STR <<
                    cur_search.len << TAB_STR;
                o.unsetf (std::ios::fixed | std::ios::scientific);
                o<< std::setprecision (2) << std::noshowpoint << e_value << TAB_STR <<
                    std::setprecision (2) << std::noshowpoint << std::fixed << sw_score << TAB_STR <<
                    cur_res->q_auto_score_ << TAB_STR <<
                    cur_res->t_auto_score_ << std::endl;
            }
        }
    }
    return true;
}

bool Search_helper_files::output_results_tab_nn (AlignResultStorage& resrec, unsigned res_per_query, unsigned alignments_no, std::vector<NN_SEQ>& f_qry, std::vector<NN_SEQ>& r_qry, WMatrix* w, std::ostream& o, bool hdr)
{

    if (hdr) o << "#Q_id\tS_id\tp_inden\tal_len\tmism\tgaps\tgap_len\tqry_beg\tqry_end\tqry_len\ttrg_beg\ttrg_end\ttrg_len\tevalue\tsw_score\tqry_auto\ttrg_auto" << std::endl;

    const char* Query_id = "";
    const char* Subject_id = "";
    double p_identity = 0;
    unsigned alignment_length = 0;
    unsigned mismatches = 0;
    unsigned gap_openings = 0;
    unsigned gap_length = 0;
    unsigned q_start = 0;
    unsigned q_end = 0;
    unsigned s_start = 0;
    unsigned s_end = 0;
    double e_value = 0;
    double sw_score = 0;

    for (unsigned query_no = 0; query_no < f_qry.size (); query_no ++)
    {
        NN_SEQ& fwd_qry = f_qry [query_no];
        NN_SEQ* rev_qry_p = r_qry.size () ? &(r_qry [query_no]) : NULL;
        // sanity check
        if (rev_qry_p && rev_qry_p->uid != fwd_qry.uid) ERR(ERR_Internal);

        unsigned sim_found = resrec.resPerQuery (fwd_qry.uid);

        if (sim_found)
        {
            query_.seek (fwd_qry.uid);
            if (!query_.next ())
                ers << "Unable to read query sequence at offset " << fwd_qry.uid << Throw;
            Query_id = query_.cur_name ();

            // get results order for current sequence
            unsigned resno = 0;
            // make results order
            ResultQueue::ElemVec& query_results = *resrec.getQueryResults (fwd_qry.uid);
            for (ResultQueue::ElemVec::iterator cur_res = query_results.begin (); cur_res != query_results.end (); cur_res ++, resno ++)
            {
                NN_SEQ cur_search;
                cur_search.uid = cur_res->uid_;

                target_.seek (cur_search.uid);
                if (!target_.next ())
                    ers << "Unable to read target sequence at offset " << cur_search.uid << Throw;
                cur_search.len = strlen (target_.cur_seq ());
                cur_search.rev = 0;
                cur_search.seq = const_cast <char*> (target_.cur_seq ());
                // inplace conversion
                n_ascii2binary (cur_search.seq, cur_search.len, cur_search.seq, 0, cur_search.len);

                // Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
                Subject_id = target_.cur_name ();

                BATCH* batches = cur_res->batches_;
                unsigned batch_no = cur_res->batch_no_;

                if (cur_res->reverse_ && !rev_qry_p) ERR (ERR_Internal);
                eval_align (cur_res->reverse_ ? *rev_qry_p : fwd_qry, cur_search, w, batches, batch_no, &p_identity, &mismatches, &alignment_length, &gap_openings, &gap_length);

                // here query IS x!
                q_start = cur_res->reverse_ ? fwd_qry.len - (batches [batch_no-1].xpos + batches [batch_no-1].len) : batches [0].xpos + 1;
                q_end = cur_res->reverse_ ? fwd_qry.len - batches [0].xpos : batches [batch_no-1].xpos + batches [batch_no-1].len;
                s_start = batches [0].ypos + 1;
                s_end = batches [batch_no-1].ypos + batches [batch_no-1].len;
                e_value = cur_res->chi2_;
                sw_score = cur_res->al_score_;

                // QueryId SubjectId P_indentity AlignmentLength Mismatches GapOpenings GapLength QueryStart QueryEnd QueryLen TargetStart TargetEnd TargetLen EValue SW_score AutoQuery AutoTarget

                o <<
                    Query_id << TAB_STR <<
                    Subject_id << TAB_STR <<
                    std::setprecision (2) << std::fixed << std::noshowpoint << p_identity << TAB_STR <<
                    alignment_length << TAB_STR <<
                    mismatches << TAB_STR <<
                    gap_openings << TAB_STR <<
                    gap_length << TAB_STR <<
                    q_start << TAB_STR <<
                    q_end << TAB_STR <<
                    fwd_qry.len << TAB_STR <<
                    s_start << TAB_STR <<
                    s_end << TAB_STR <<
                    cur_search.len << TAB_STR;
                o.unsetf (std::ios::fixed | std::ios::scientific);
                o<< std::setprecision (2) << std::noshowpoint << e_value << TAB_STR <<
                    std::setprecision (2) << std::fixed << std::noshowpoint << sw_score << TAB_STR <<
                    cur_res->q_auto_score_ << TAB_STR <<
                    cur_res->t_auto_score_ << std::endl;
            }
        }
    }
    return true;
}


ulonglong Search_helper_files::current () const
{
    return target_.cur_pos ();
}

ulonglong Search_helper_files::total () const
{
    return target_.tot_len ();
}

unsigned Search_helper_files::skipped () const
{
    return skipped_;
}
