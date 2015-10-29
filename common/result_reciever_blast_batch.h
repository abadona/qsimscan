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

#ifndef __result_reciever_blast_batch_h__
#define __result_reciever_blast_batch_h__

class ResultReciever_blast_batch
{
public:
    virtual bool match_found  (NN_SEQ& xseq, NN_SEQ& yseq, BATCH* batches, int batch_no, int matches) = 0;
};


#endif // __result_reciever_blast_batch_h__
