
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

#include <rerror.h>
#include "scan_results.h"
#include "align.h"
#include "filters.h"

ScanResults :: ScanResults (ALIGN* align, int keep_per_query, SEQ_SCAN_TYPE st, WEIGHTS<int, 24>* wm, double score_thresh, int rep_len, double rep_perc)
:
AlignResultStorage (keep_per_query),
_aligner        (align),
_score_thresh   (int (score_thresh)),
_rep_percent    (int (rep_perc)),
_rep_len        (rep_len),
_total_found    (0),
_passed_repeats (0),
_scan_type      (st),
_wm             (wm)
{
    _rep_buf = NULL;
    if (rep_len)
    {
        try { _rep_buf = new int [rep_len];
        } catch (std::bad_alloc&) {}
        if (!_rep_buf) ERR(NOEMEM);
    }
}

ScanResults :: ~ScanResults ()
{
    if (_rep_buf)
        delete [] _rep_buf;
}

// this is actually the dispatcher
void ScanResults :: match_found (SEQ* query_seq,  SEQ* search_seq, int score)
{
    switch (_scan_type)
    {
    case SST_NN: process_match_nn ((NN_SEQ*) query_seq, (NN_SEQ*) search_seq, score); break;
    case SST_NN_REV: process_match_nn ((NN_SEQ*) search_seq, (NN_SEQ*) query_seq, score); break;
    case SST_AA: process_match_aa ((AA_SEQ*) query_seq, (AA_SEQ*) search_seq, score); break;
    case SST_AA_REV: process_match_aa ((AA_SEQ*) search_seq, (AA_SEQ*) query_seq, score); break;
    case SST_AN: process_match_an ((AA_SEQ*) query_seq, (NA_SEQ*) search_seq, score); break;
    case SST_NA: process_match_na ((NA_SEQ*) query_seq, (AA_SEQ*) search_seq, score); break;
    case SST_NA_REV: process_match_na_rev ((AA_SEQ*) query_seq, (NA_SEQ*) search_seq, score); break;
    }
}

void ScanResults :: process_match_nn (NN_SEQ* query_seq,  NN_SEQ* search_seq, int score)
{
    _total_found ++;

    if (!_aligner)
        return;

    // perform the alignment
    int al_score = _aligner -> align (*query_seq, *search_seq);
    if (al_score < _score_thresh) return;
    int batch_no = _aligner -> backtrace (_batches, MAX_BATCHES);


    // recalculate the score
    int matches = count_matches (query_seq->seq, search_seq->seq, _batches, batch_no);
    int tot_blen = len_batches (_batches, batch_no);
    float chi2score = nu_all_chi2_b (query_seq->seq, search_seq->seq, _batches, batch_no);

    // filter results by repeats
    if (_rep_len)
    {
        int max_idx = nu_offs_score_b (query_seq->seq, search_seq->seq, _batches, batch_no, _rep_len, _rep_buf);
        if (!matches || (_rep_buf [max_idx] * 100) / matches > _rep_percent)
            return;
    }

    _passed_repeats ++;

    // save the result
    add_result (query_seq->uid, search_seq->uid, query_seq->rev?true:false, al_score, score, chi2score, 0, 0, tot_blen, tot_blen, batch_no, _batches, search_seq->seq, search_seq->len);
}


void ScanResults :: process_match_aa (AA_SEQ* query_seq, AA_SEQ* search_seq, int score)
{
    _total_found ++;

    if (!_aligner)
        return;

    // perform the alignment
    int al_score = _aligner -> align (*query_seq, *search_seq);
    if (al_score < _score_thresh) return;
    int batch_no = _aligner -> backtrace (_batches, MAX_BATCHES);

    int q_score;
    prot_score (query_seq->seq, search_seq->seq, _batches, batch_no, _wm, &q_score);

    _passed_repeats ++;

    // save the result
    add_result (query_seq->uid, search_seq->uid, query_seq->rev?true:false, al_score, score, al_score, 0, 0, q_score, 0, batch_no, _batches);
}

void ScanResults :: process_match_an (AA_SEQ* query_seq, NA_SEQ* search_seq, int score)
{
    _total_found ++;

    if (!_aligner)
        return;

    // perform the alignment
    int al_score = _aligner -> align_na (*query_seq, *search_seq);
    if (al_score < _score_thresh) return;
    int batch_no = _aligner -> backtrace (_batches, MAX_BATCHES);

    _passed_repeats ++;

    int q_score;
    int sc = an_score (query_seq->seq, search_seq->seq, _batches, batch_no, false, _wm, &q_score);

    // printf ("AN: bno = %d, al_score = %d, sc = %d, q_score = %d\n", batch_no, al_score, sc, q_score);

    // save the result
    add_result (query_seq->uid, search_seq->uid, query_seq->rev?true:false, al_score, score, al_score, 0, 0, q_score, 0, batch_no, _batches);
}
void ScanResults :: process_match_na (NA_SEQ* query_seq, AA_SEQ* search_seq, int score)
{
    _total_found ++;

    if (!_aligner)
        return;

    // perform the alignment
    int al_score = _aligner -> align_na (*search_seq, *query_seq);
    if (al_score < _score_thresh) return;
    int batch_no = _aligner -> backtrace (_batches, MAX_BATCHES);

    _passed_repeats ++;

    int q_score, s_score;
    int sc = an_score (search_seq->seq, query_seq->seq, _batches, batch_no, false, _wm, &s_score, &q_score);

    // printf ("AN: bno = %d, al_score = %d, sc = %d, q_score = %d\n", batch_no, al_score, sc, q_score);

    // save the result
    add_result (query_seq->uid, search_seq->uid, query_seq->rev?true:false, al_score, score, al_score, q_score, 0, 0, 0, batch_no, _batches);
}

void ScanResults :: process_match_na_rev (AA_SEQ* search_seq, NA_SEQ* query_seq, int score)
{
    _total_found ++;

    if (!_aligner)
        return;

    // perform the alignment
    int al_score = _aligner -> align_na (*search_seq, *query_seq);
    if (al_score < _score_thresh) return;
    int batch_no = _aligner -> backtrace (_batches, MAX_BATCHES);

    _passed_repeats ++;

    int q_score, s_score;
    int sc = an_score (search_seq->seq, query_seq->seq, _batches, batch_no, false, _wm, &s_score, &q_score);

    // printf ("\nAN: bno = %d, al_score = %d, sc = %d, q_score = %d", batch_no, al_score, sc, q_score);

    // save the result
    add_result (query_seq->uid, search_seq->uid, query_seq->rev?true:false, al_score, score, al_score, 0, 0, q_score, 0, batch_no, _batches);
}


