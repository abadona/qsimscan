
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

#ifndef __merging_result_storage_h__
#define __merging_result_storage_h__

#include <vector>
#include <map>
#include "align_result_storage.h"
#include "sim_merger_base.h"


class MergingResultStorage : public AlignResultStorage
{
    SimMergerBase& merger_;

    typedef std::vector < AlignResult > ARVect;
    typedef std::map < longlong, ARVect > QryResults;

    QryResults accum_;
    unsigned accum_no_;
    longlong cur_sid_;

    void clear_accum ();
public:
    MergingResultStorage (SimMergerBase& merger, unsigned capacity);
    bool add_result (longlong qid, longlong sid, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, BATCH* batches, const char* binsubj = NULL, QWORD subjid = 0, DWORD subjlen = 0);
    void flush ();
    bool reset ();
};

#endif
