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

#ifdef _MSC_VER
#pragma warning (disable: 4786)
#endif

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

bool AlignResultStorage :: add_result (longlong query_id, const AlignResult& r)
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

