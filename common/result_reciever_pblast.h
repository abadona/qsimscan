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

#ifndef __result_reciever_pblast_h__
#define __result_reciever_pblast_h__
#pragma warning (disable:4786)

#include "biosequence.h"

class ResultReciever_pblast
{
public:
    virtual bool match_found  (SEQ& xseq, SEQ& yseq, BATCH* batches, int batch_no, double score, double q_auto_score, double t_auto_score) = 0;
};


#endif
