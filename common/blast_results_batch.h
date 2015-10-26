
//////////////////////////////////////////////////////////////////////////////
// This software module is developed by SCIDM team in 1999-2015.
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

#ifndef __blast_results_batch_h__
#define __blast_results_batch_h__

#include "sequtil.h"
#include "biosequence.h"
#include "result_reciever_blast_batch.h"
#include "merging_result_storage.h"

class Search_helper_frag;

class BlastResultsBatch : public ResultReciever_blast_batch, public MergingResultStorage
{
    static SimMergerBase null_merger_;

    int total_found_;
    int passed_repeats_;

    int  rep_percent_;
    int  rep_len_;
    int* rep_buf_;

    bool triangle_only_;

public:
    BlastResultsBatch (int keep_per_query, int rep_len = 0, int rep_perc = 0, bool triangle_only = false, SimMergerBase& merger = null_merger_);
    virtual ~BlastResultsBatch () {}
    virtual bool match_found  (NN_SEQ& xseq, NN_SEQ& yseq, BATCH* batches, int batch_no, int matches);
    int totalFound () const {return total_found_;}
    int passedRepeats () const {return passed_repeats_;}
};


#endif // __blast_results_frag_batch_h__
