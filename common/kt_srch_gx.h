
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

#ifndef __kt_srch_h__
#define __kt_srch_h__

#include <platform.h>
#include "biosequence.h"
#include "result_reciever_blast.h"
#include "result_reciever_blast_batch.h"

#define MAX_KTUP 12
#define DEF_K_WT 100

#pragma pack (push, 1)

//align at 32 byte boundary
struct DIAG_ENTRY
{
	longlong xuid;
    int xrev;
	int ynum;
	int xpos;
	int score;
	int b_end;
	int pad1;
};

struct KTUP_ENTRY
{
	int ynum;
	int ypos;
};

//align at 16 byte boundary
struct KTUP_INFO
{
	int cnt;
	int wt;
	KTUP_ENTRY*	ptr;
	int pad1;
};

struct SEQ_INFO
{
	int     start;
    int     len;
	longlong uid;
    int     rev;
	char*   seq;
};

#pragma pack (pop)

class KT_SEARCH
{
	// output mode
    int batches_out;

    //K-tuple info
	int kt_size;
	unsigned kt_mask;

	//K-tuple info table
	int ki_len;
	KTUP_INFO* ki;

	//K-tuple position array
	int ke_len;
	KTUP_ENTRY* ke;

	//diag array
	int di_len;
	DIAG_ENTRY* di;

	//K-tuple history array
	int kt_last;
	unsigned prev_ktup[8];

	//Y sequences info
	int y_tnum;				//total number of Y sequences loaded
	int y_tlen;				//total length of Y sequences
	SEQ_INFO* yi;

	//search limits
	int max_x_len;
	int max_y_num;
	int max_y_len;

    //status
    int armed;

	//runinig sequence variables
	char *xseq, *yseq;
	int xlen, ylen;
	int xpos, ypos;
	int ynum;
    longlong xuid;
    int xrev, yrev;
	int diag;

	//pass 2 - delimit parameters
	int offs_lev;			//set to b_thresh - 5
	int stop_lev;			//set to 3 * offs_lev
	int penalty;			//set to (l_thresh * (b_thresh_s - b_thresh_l)) / 100;

	void scan_l1 ();
	void scan_l1_fast ();
	void scan_l2 ();
	void scan_l2_g ();
	void scan_l2_gb_1 ();
	void scan_l2_gb_2 ();
	void ki_count (char* seq, int len);
	void ke_add (char* seq, int len, int ynum, int ystart);

	//result reciever
	ResultReciever_blast* batch;
    ResultReciever_blast_batch* batch_processor;

    void init_vars (int max_ynum, int *wts);

public:
	//SEARCH PARAMETERS
	//No of adjacent diagonal to check
	int max_offs;

	//K-tuple repeat elimination
	int rep_del;			//0-no filtering, 8-look up 8 prev K-tuples

	//pass 1 parameters
	int fast_mode;			//fast, high similarity level >80% search mode (K-tup xstep = 4)
	int apprx_match;		//allow 1 base mismatch in K-tuple
	int k_thresh;			//diagonal score rejection threshold
	int k_gip;				//pass 1 gap initialisation penalty
	int k_gep;				//pass 1 gap extention penalty

	//pass 2 parameters
	int b_thresh_s;			//batch similarity level (%) at minimum length (l_thresh)
	int b_thresh_l;			//batch similarity level (%) at infinite length
	int l_thresh;			//minimum batch length
	int g_period;			//mimimum gap period (in nucleotides)

	//misc globals
	int k_matches;			//internal pass 1 result count
	int b_matches;			//internal pass 2 result count

	//interface functions
		KT_SEARCH 	(int kt_sz,		//k-tuple size
					 int *wts,		//k-tuple weights array
					 int max_ynum,	//maximum number of Y sequences
					 int max_total_ylen, //maximum sum of y lengths
					 int max_xlen,	//maximum X seq lenght
					 ResultReciever_blast* res_reciever);

        KT_SEARCH 	(int kt_sz,		//k-tuple size
					 int *wts,		//k-tuple weights array
					 int max_ynum,	//maximum number of Y sequences
					 int max_total_ylen, //maximum sum of y lengths
					 int max_xlen,	//maximum X seq lenght
					 ResultReciever_blast_batch* res_reciever);

		~KT_SEARCH 	();

	int	add_yseq 	(NN_SEQ& yseq);
	int	init	 	();
	int  search	    (NN_SEQ& xseq);
	void reset_y	();


	// debug helpers
	int _diags_found;
	int _ktups_found;
	int _nuc_scanned;

};

#endif // __kt_srch_h__
