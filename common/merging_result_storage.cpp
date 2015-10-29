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

MergingResultStorage::MergingResultStorage (SimMergerBase& merger, unsigned capacity)
:
AlignResultStorage (capacity),
merger_ (merger),
cur_sid_ (-1),
accum_no_ (0)
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

void MergingResultStorage::flush (const char* tseq)
{
    if (!accum_no_)
        return;
    QryResults::iterator qi;
    for (qi = accum_fwd_.begin (); qi != accum_fwd_.end (); qi ++)
    {
        ARVect& alignments = (*qi).second;
        merger_.merge (alignments, tseq);
        for (ARVect::iterator ri = alignments.begin (); ri != alignments.end (); ri ++)
            AlignResultStorage::add_result ((*qi).first, *ri);
    }
    for (qi = accum_rev_.begin (); qi != accum_rev_.end (); qi ++)
    {
        ARVect& alignments = (*qi).second;
        merger_.merge (alignments, tseq);
        for (ARVect::iterator ri = alignments.begin (); ri != alignments.end (); ri ++)
            AlignResultStorage::add_result ((*qi).first, *ri);
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
