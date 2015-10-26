
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
