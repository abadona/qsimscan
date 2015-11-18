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

#include "merging_result_storage.h"

MergingResultStorage::MergingResultStorage (SimMergerBase& merger, unsigned capacity, unsigned res_per_target)
:
AlignResultStorage (capacity),
merger_ (merger),
res_per_target_ (res_per_target),
accum_no_ (0),
cur_sid_ (-1)
{
}

bool MergingResultStorage::add_result (longlong qid, longlong sid, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, BATCH* batches, const char* binsubj, QWORD subjid, DWORD subjlen)
{
    if (cur_sid_ != sid)
    {
        // flush (binsubj);
        cur_sid_ = sid;
    }
    AlignResult r (sid, reverse, al_score, score, chi2, evalue, bitscore, q_auto_score, t_auto_score, batch_no, batches, binsubj, subjid, subjlen);
    if (reverse)
        accum_rev_ [qid].push_back (r);
    else
        accum_fwd_ [qid].push_back (r);
    accum_no_ ++;
    return true;
}

struct AlignResultWrapper
{
    const AlignResult* align_;
    longlong query_;

    AlignResultWrapper ()
    :
    align_ (NULL),
    query_ (0)
    {
    }
    AlignResultWrapper (const AlignResult& align, longlong query)
    :
    align_ (&align),
    query_ (query)
    {
    }
    AlignResultWrapper (const AlignResultWrapper& oth)
    :
    align_ (oth.align_),
    query_ (oth.query_)
    {
    }
    AlignResultWrapper& operator = (const AlignResultWrapper& oth)
    {
        align_ = oth.align_;
        query_ = oth.query_;
    }
    bool operator < (const AlignResultWrapper& other) const
    {
        return *align_ < *(other.align_);
    }
    bool operator == (const AlignResultWrapper& other) const
    {
        return *align_ == *(other.align_);
    }
};
typedef PQueue <AlignResultWrapper> ARWQueue;


void MergingResultStorage::flush (const char* tseq)
{
    if (!accum_no_)
        return;
    QryResults::iterator qi;
    ARWQueue besthits (res_per_target_);
    for (qi = accum_fwd_.begin (); qi != accum_fwd_.end (); qi ++)
    {
        ARVect& alignments = (*qi).second;
        merger_.merge (alignments, tseq);
        for (ARVect::iterator ri = alignments.begin (); ri != alignments.end (); ri ++)
        {
            if (res_per_target_)
                besthits.push (AlignResultWrapper (*ri, (*qi).first));
            else
                AlignResultStorage::add_result ((*qi).first, *ri);
        }
    }
    for (qi = accum_rev_.begin (); qi != accum_rev_.end (); qi ++)
    {
        ARVect& alignments = (*qi).second;
        merger_.merge (alignments, tseq);
        for (ARVect::iterator ri = alignments.begin (); ri != alignments.end (); ri ++)
        {
            if (res_per_target_)
                besthits.push (AlignResultWrapper (*ri, (*qi).first));
            else
                AlignResultStorage::add_result ((*qi).first, *ri);
        }
    }
    if (res_per_target_)
    {
        const ARWQueue::ElemVec& wrappers = besthits.data ();
        for (unsigned pp = 0, sent = besthits.size (); pp != sent; ++pp)
        {
            const AlignResultWrapper& wrapper = wrappers [pp];
            AlignResultStorage::add_result (wrapper.query_, *wrapper.align_);
        }
    }
    clear_accum ();
}

void MergingResultStorage::clear_accum ()
{
    QryResults::iterator ii;
    for (ii = accum_fwd_.begin (); ii != accum_fwd_.end (); ii ++)
        (*ii).second.clear ();
    for (ii = accum_rev_.begin (); ii != accum_rev_.end (); ii ++)
        (*ii).second.clear ();
    accum_no_ = 0;
}

bool MergingResultStorage::reset ()
{
    clear_accum ();
    return AlignResultStorage::reset ();
}
