
//////////////////////////////////////////////////////////////////////////////
// This software module is developed by SCIDM team in 1999-2008.
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
