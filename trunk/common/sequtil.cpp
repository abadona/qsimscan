
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

#define __sequtil_cpp__

#include "sequtil.h"
#include <string.h>
#include <iomanip>
#include <ios>

void n_revert_seq (char* dest, const char* src, int len, int beg)
{
    for (int ii = 0; ii < len; ii ++)
        put_base (dest, len - ii - 1, ((~get_base (src, ii+beg)) & 0x3));
}

bool n_ascii2binary (char* dest, int destlen, const char* source, int start, int len)
{
    if (destlen < ((len + 3) >> 2))
        return false;

    for (int i = start; i < start + len; i ++)
        put_base (dest, i - start, char2base (source [i]));

    return true;
}

bool n_binary2ascii (char* dest, int destlen, const char* source, int start, int len)
{
    if (destlen < len)
        return false;

    for (int i = start; i < start + len; i ++)
        dest [i - start] = base2char (get_base (source, i));

    return true;
}

void n_countnucs (int& a, int& g, int& c, int& t, const char* src, int start, int len)
{
    a = g = t = c = 0;
    if (len == -1) len = strlen (src) - start;

    for (int i = start; i < start + len; i ++)
        switch (tolower (src [i]))
        {
        case 't':   t ++; break;
        case 'a':   a ++; break;
        case 'g':   g ++; break;
        case 'c':   c ++; break;
        }
}

std::ostream& output_seq (std::ostream& o, const char* buffer, int buflen, int chars_in_line, bool counts, bool decades, int first_off)
{
    int cnt = 0;

    while (cnt < buflen)
    {
        if (cnt % chars_in_line == 0)
        {
            if (cnt) o << std::endl;
            if (counts) o << std::setw (first_off - 1) << std::right << cnt + 1 << " " << std::setw (0);
            else o << std::setw (first_off) << /*std::right <<*/ "" << std::setw (0);
        }
        else if (decades && (cnt % 10 == 0))
        {
            o << " ";
        }
        o << buffer [cnt];
        cnt ++;
    }
    o << std::endl;
    return o;
}

std::ostream& fasta_output_seq (std::ostream& o, const char* buffer, int buflen, int chars_per_line)
{
    return output_seq (o, buffer, buflen, chars_per_line, false, false, 0);
}

std::ostream& gb_output_seq (std::ostream& o, const char* buffer, int buflen)
{
    return output_seq (o, buffer, buflen, 60, true, true, 10);
}

                        // A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z   *
const char aa2number [] = {0,  23, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23, 14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23 };

                        // 0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23
const char number2aa [] = {'A','R','N','D','C','Q','E','G','H','I','L','K', 'M','F','P','S','T','W','Y','V',' ','Z','X','*' };



bool a_ascii2binary (char* dest, int destlen, const char* source, int start, int len)
{
    if (destlen < len)
        return false;

    for (int i = start; i < start + len; i ++)
        dest [i - start] = char2aa (source [i]);

    return true;
}

bool a_binary2ascii (char* dest, int destlen, const char* source, int start, int len)
{
    if (destlen < len)
        return false;

    for (int i = start; i < start + len; i ++)
        dest [i - start] = aa2char (source [i]);

    return true;
}


/*
int char2aa (char c)
{
    switch (toupper(c))
    {
        case 'A':
            return 0;
        case 'R':
            return 1;
        case 'N':
            return 2;
        case 'D':
            return 3;
        case 'C':
            return 4;
        case 'Q':
            return 5;
        case 'E':
            return 6;
        case 'G':
            return 7;
        case 'H':
            return 8;
        case 'I':
            return 9;
        case 'L':
            return 10;
        case 'K':
            return 11;
        case 'M':
            return 12;
        case 'F':
            return 13;
        case 'P':
            return 14;
        case 'S':
            return 15;
        case 'T':
            return 16;
        case 'W':
            return 17;
        case 'Y':
            return 18;
        case 'V':
            return 19;
        case 'B':
            return 20;
        case 'Z':
            return 21;
        case 'X':
            return 22;
        default:
            return 23;
    }
}
*/

int swap_batches (BATCH* batches, int batch_no)
{
    int totlen = 0;
    int xpos;
    BATCH* b = batches;
    for (int idx = 0; idx < batch_no; idx ++)
    {
        xpos = b->xpos;
        b->xpos = b->ypos;
        b->ypos = xpos;
        totlen += b->len;
        b ++;
    }
    return totlen;
}

int len_batches (const BATCH* batches, int batch_no)
{
    int totlen = 0;
    for (int idx = 0; idx < batch_no; idx ++)
        totlen += batches [idx].len;
    return totlen;
}

int count_matches (const char* xseq, const char* yseq, const BATCH* batches, int batch_no)
{
    int matches = 0;
    const BATCH* curb;
    int xp, yp;
    int pos;
    for (int idx = 0; idx < batch_no; idx ++)
    {
        curb = batches + idx;
        xp = curb->xpos;
        yp = curb->ypos;
        for (pos = 0; pos < curb->len; pos++, xp++, yp++)
        {
            if (get_base (xseq, xp) == get_base (yseq, yp))
                matches ++;
        }
    }
    return matches;
}
int count_gaps (const BATCH* batches, int batch_no)
{
    int gaplen = 0;
    int lastx, lasty;
    const BATCH* curb;
    for (int idx = 0; idx < batch_no; idx ++)
    {
        curb = batches + idx;
        if (idx)
            gaplen += (curb->xpos - lastx) + (curb->ypos - lasty);
        lastx = curb->xpos + curb->len;
        lasty = curb->ypos + curb->len;
    }
    return gaplen;
}
