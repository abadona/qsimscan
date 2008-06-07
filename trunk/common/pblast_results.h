
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

#ifndef __pblast_results_h__
#define __pblast_results_h__

#pragma warning (disable : 4786)
#include <platform.h>
#include "result_reciever_pblast.h"
#include "merging_result_storage.h"
#include "gaplibshim.h"


class PblastResults : public ResultReciever_pblast, public MergingResultStorage
{
    static SimMergerBase null_merger_;
    int total_found_;
    int min_len_;
    double min_score_;
    WMatrixType *w_;
	bool eval_eval_;
    // Cached K, lambda, H and alpha for target sequence
    // we usually called several times with same target sequence,
    // so to save on calculation of base statistics for the same
    // sequence we cache its uid and calculated parameters here
    longlong uid_cached_;
    double K0mean_, lambda0mean_, Hmean_, alpha_mean_;
    ResFreqType trg_freq_;
public:
    PblastResults (double min_score, int min_len, WMatrixType *w, bool eval_eval = true, unsigned res_per_qry = 0, SimMergerBase& merger = null_merger_);
    bool match_found  (SEQ& query_seq, SEQ& target_seq, BATCH* batches, int batch_no, double score, double q_auto_score, double t_auto_score);
    int totalFound ();
};


#endif
