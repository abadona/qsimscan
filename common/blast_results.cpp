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

#pragma warning (disable : 4786)
#include "blast_results.h"

blast_results :: blast_results (ALIGN* align, int keep_per_query, double score_thresh, int rep_len, double rep_perc)
:
ScanResults (align, keep_per_query, SST_NN, NULL, score_thresh, rep_len, rep_perc)
{
}

bool blast_results :: match_found  (NN_SEQ& xseq, int xpos, NN_SEQ& yseq, int ypos, int len, int matches)
{
    if (!(curx == xseq))
    {
        curx = xseq;
        seen_y.erase (seen_y.begin (), seen_y.end ());
    }
    else
    {
        req t (yseq);
        if (!seen_y.insert (t).second) return true;
    }

    ScanResults::match_found (&yseq, &xseq, matches);
    return true;
}
