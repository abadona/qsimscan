
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

#ifndef __align_result_storage_h__
#define __align_result_storage_h__

#include <map>
#include <platform.h>
#include <pqueue.h>
#include "align_result.h"

typedef PQueue<AlignResult> ResultQueue;
typedef std::map<longlong, ResultQueue> ResultSet;

class AlignResultStorage
{
private:

    ResultSet heads_;
    unsigned total_stored_;
    unsigned max_per_query_;

public:

    AlignResultStorage (unsigned capacity);
    virtual ~AlignResultStorage ();

    bool add_result (longlong qid, longlong sid, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, BATCH* batches);

    unsigned totalStored () {return total_stored_;}
    unsigned resPerQuery (longlong id);
    ResultQueue::ElemVec* getQueryResults (longlong id);
    ResultSet& heads () {return heads_;}

    bool reset ();
};


#endif // __align_result_storage_h__
