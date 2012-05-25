
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

#include "biosequence.h"
#include "rerror.h"

const char* const aa2str[AANUM] = {\
"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",\
"His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser",\
"Thr", "Trp", "Tyr", "Val", "BBB", "ZZZ", "XXX", "STP"};

// standard genetic code - kept as const data to avoid loading from database

const char* const gcode_c_rev = "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF";

const char* const gcode_c     = "KEQ*RGR*TAPSIVLLKEQ*RGRWTAPSMVLLNDHYSGRCTAPSIVLFNDHYSGRCTAPSIVLF";

const int gcode_b_rev[64] = {
11, 11, 2, 2, 1, 1, 15, 15, 16, 16, 16, 16, 9, 12, 9, 9, 6, 6, 3, 3, 7, 7, 7, 7, 0, 0, 0, 0,
19, 19, 19, 19, 5, 5, 8, 8, 1, 1, 1, 1, 14, 14, 14, 14, 10, 10, 10, 10, 23, 23, 18, 18, 23,
17, 4, 4, 15, 15, 15, 15, 10, 10, 13, 13};

const int gcode_b[64] = {
11, 6, 5, 23, 1, 7, 1, 23, 16, 0, 14, 15, 9, 19, 10, 10, 11, 6, 5, 23, 1, 7, 1, 17, 16, 0,
14, 15, 12, 19, 10, 10, 2, 3, 8, 18, 15, 7, 1, 4, 16, 0, 14, 15, 9, 19, 10, 13, 2, 3, 8, 18,
15, 7, 1, 4, 16, 0, 14, 15, 9, 19, 10, 13};


void reverse_seqs (std::vector <NN_SEQ>& src, std::vector <NN_SEQ>& dest)
{
    for (int ii = 0; ii < src.size (); ii ++)
    {
        NN_SEQ& cur_src = src [ii];
        NN_SEQ rs;
        rs.uid = cur_src.uid;
        rs.len = cur_src.len;
        rs.rev = 1;
        rs.seq = NULL;
        try
        {
            rs.seq = new char [(rs.len+3)>>2];
        }
        catch (std::bad_alloc&)
        {
        }
        if (!rs.seq)
            ERR(NOEMEM);
        n_revert_seq (rs.seq, cur_src.seq, rs.len);

        dest.push_back (rs);
    }
}

void reverse_seqs (std::vector <NA_SEQ>& src, std::vector <NA_SEQ>& dest)
{
    for (int ii = 0; ii < src.size (); ii ++)
    {
        NA_SEQ rs;
        rs.uid = src [ii].uid;
        rs.len = src [ii].len;
        rs.rev = 1;
        rs.seq = NULL;
        try
        {
            rs.seq = new char [(rs.len+3)>>2];
        }
        catch (std::bad_alloc&)
        {
        }
        if (!rs.seq)
            ERR(NOEMEM);
        n_revert_seq (rs.seq, src [ii].seq, rs.len);

        dest.push_back (rs);
    }
}

