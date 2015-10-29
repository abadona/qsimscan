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

#ifndef __blast_results_h__
#define __blast_results_h__

#include <pqueue.h>
#include <set>
#include "scan_results.h"
#include "result_reciever_blast.h"
#include "blast_results_req.h"


class blast_results	: public ResultReciever_blast, public ScanResults
{
    req curx;
    std::set <req> seen_y;

public:
	blast_results (ALIGN* align, int keep_per_query, double score_thresh = 0, int rep_len = 0, double rep_perc = 0);

	virtual bool match_found  (NN_SEQ& xseq, int xpos, NN_SEQ& yseq, int ypos, int len, int matches);
};

#endif // __blast_results_h__
