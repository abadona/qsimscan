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

protected:
    bool add_result (longlong query_id, const AlignResult& r);

public:

    AlignResultStorage (unsigned capacity);
    virtual ~AlignResultStorage ();

    virtual bool add_result (longlong qid, longlong sid, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, BATCH* batches, const char* binsubj = NULL, QWORD subjid = 0, DWORD subjlen = 0);
    unsigned totalStored () {return total_stored_;}
    unsigned resPerQuery (longlong id);
    ResultQueue::ElemVec* getQueryResults (longlong id);
    ResultSet& heads () {return heads_;}
    virtual void flush (const char* tseq) {}
    virtual bool reset ();
};


#endif // __align_result_storage_h__
