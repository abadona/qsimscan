
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
// For any questions please contact SciDM team by email at scidmteam@yahoo.com
//////////////////////////////////////////////////////////////////////////////

#pragma warning (disable : 4786)
#include "sequtil.h"
#include "align_result_storage.h"

AlignResultStorage :: AlignResultStorage (unsigned capacity)
{
    max_per_query_ = capacity;
    total_stored_ = 0;
}
AlignResultStorage :: ~AlignResultStorage ()
{
}

bool AlignResultStorage :: add_result (longlong query_id, longlong search_id, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, BATCH* batches, const char* binsubj, QWORD subjid, DWORD subjlen)
{
    // if this is a first result for given y, expand the queue to desired size
    ResultQueue* q;

    ResultSet::iterator itr = heads_.find (query_id);
    if (itr == heads_.end ())
    {
        q = &(heads_ [query_id]);
        q->setCapacity (max_per_query_);
    }
    else
    {
        // check for the need to add (to avoid extra allocation of batches buffer)
        q = &(itr->second);
        if (q->size () == max_per_query_)
        {
            AlignResult d (search_id, reverse, al_score, score, chi2, evalue, bitscore, 0, 0, 0, NULL, NULL); // no batch array allocation / copying here
            if (d < q->top ())
                return false;  // result Ok but no space left
        }
    }

    // create the Result object
    AlignResult r (search_id, reverse, al_score, score, chi2, evalue, bitscore, q_auto_score, t_auto_score, batch_no, batches, binsubj, subjid, subjlen);

    // add to heap.
    q->push (r);

    total_stored_ ++;
    return true;
}

bool AlignResultStorage :: add_result (longlong query_id, AlignResult& r)
{
    // if this is a first result for given y, expand the queue to desired size
    ResultQueue* q;

    ResultSet::iterator itr = heads_.find (query_id);
    if (itr == heads_.end ())
    {
        q = &(heads_ [query_id]);
        q->setCapacity (max_per_query_);
    }
    else
        q = &(itr->second);

    // add to heap. This will call default assignment operator that will transfer batches array ownership to the object in queue.
    if (q->push (r))
        total_stored_ ++;
    return true;
}

unsigned AlignResultStorage :: resPerQuery (longlong query_id)
{
    ResultSet::iterator itr = heads_.find (query_id);
    if (itr == heads_.end ())
        return 0;
    else
        return itr->second.size ();
}

/*
AlignResult* AlignResultStorage :: getFirstResult (longlong query_id)
{
    ResultSet::iterator itr = heads_.find (query_id);
    if (itr == heads_.end ())
        return NULL;

    if (itr->second.size ())
        return &(itr->second.top ());
    else
        return NULL;
}

AlignResult* AlignResultStorage :: getNextResult (longlong query_id)
{
    ResultSet::iterator itr = heads_.find (query_id);
    if (itr == heads_.end ())
        return NULL;

    if (itr->second.size ())
    {
        itr->second.pop ();
        if (itr->second.size ())
            return &(itr->second.top ());
    }
    return NULL;
}
*/

ResultQueue::ElemVec* AlignResultStorage :: getQueryResults  (longlong query_id)
{
    ResultSet::iterator itr = heads_.find (query_id);
    if (itr == heads_.end ())
        return NULL;
    else
        return &(*itr).second.sort ();
}

bool AlignResultStorage :: reset ()
{
    heads_.clear ();
    total_stored_ = 0;
    return true;
}

