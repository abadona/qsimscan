
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

#include "print_batches.h"
#include <limits.h>
#include <iomanip>
#include <sciminmax.h>
#include <rerror.h>
#include "weights.h"
#include <cassert>

using namespace std;

/*
expanded format -

NA 1x insert:

LysLys--LysLys
...:::  ...|||
LysLys--LysLys
agtgtgattccgta
^f1		^f2

NA 3x insert:

LysLys---LysLys
...:::   ...|||
LysLysLysLysLys
agtgtgatctccgta
^f1		^f2

AA 1x insert:

LysLysLysLysLys
...:::   ...|||
LysLys---LysLys
agtgtg---attccg
^f1		 ^f1

compact format -

ASDFGHJKL
.: |..:|
ASDFGHJKL

alternative ???

ASD-FGHJKL
.: |..:|
ASD-FGHJKL
aagttatgta
ctc	gcgctt
tga	gtacag
1	2
*/


//create text batch AAxAA, NNxNN
void print_batches (SEQ& xseq, SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, ostream& stream, int margin, int width)
{
    int slen = 0;
    int blen = 0;
    int x = b_ptr->xpos;
    int y = b_ptr->ypos;
    int xstart = x;
    int ystart = y;
    int xc, yc;
    char s[3][256];

    stream << endl << setiosflags (ios::left);
    while (b_cnt > 0)
    {
        xc = xseq.get_code (x);
        yc = yseq.get_code (y);

        //special handling may be needed for (x < b_ptr->xpos && y < b_ptr->xpos)
        //X insert
        if (x < b_ptr->xpos)
        {
            s[0][slen] = xseq.code2char (xc);
            s[1][slen] = ' ';
            s[2][slen] = '-';
            x++, slen++;
        }
        //Y insert
        else if (y < b_ptr->ypos)
        {
            s[0][slen] = '-';
            s[1][slen] = ' ';
            s[2][slen] = yseq.code2char (yc);
            y++, slen++;
        }
        //emit text batch
        else if (blen < b_ptr->len)
        {
            s[0][slen] = xseq.code2char (xc);
            s[1][slen] = w->get_char (xc, yc);
            s[2][slen] = yseq.code2char (yc);
            x++, y++, slen++, blen++;
        }
        else
            blen = 0, b_cnt--, b_ptr++;

        //print accumulated lines
        if ((slen > width - 7) || b_cnt <= 0)
        {
            //null terminate all strings
            for (int i = 0; i < 3; i++)
                s[i][slen] = 0;

            int xdisp = xseq.rev ? xseq.len - xstart - 1 : xstart;
            int ydisp = yseq.rev ? yseq.len - ystart - 1 : ystart;
            stream << endl << setw (margin) << "" << setw (6) << xdisp  << setw (0) << s[0];
            stream << endl << setw (margin) << "" << setw (6) << " "    << setw (0) << s[1];
            stream << endl << setw (margin) << "" << setw (6) << ydisp  << setw (0) << s[2];
            stream << endl;

            xstart = x, ystart = y, slen = 0;
        }
    }
    stream << flush;
}

void print_batches_ascii_subj (SEQ& xseq, const char* ascy, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin, int width)
{
    int slen = 0;
    int blen = 0;
    int x = b_ptr->xpos;
    int y = b_ptr->ypos;
    int xstart = x;
    int ystart = y;
    int y0 = y;
    char xc, yc;
    char s[3][256];

    stream << endl << setiosflags (ios::left);
    while (b_cnt > 0)
    {
        xc = xseq.code2char (xseq.get_code (x));
        yc = ascy [y - y0];

        //special handling may be needed for (x < b_ptr->xpos && y < b_ptr->xpos)
        //X insert
        if (x < b_ptr->xpos)
        {
            s[0][slen] = xc;
            s[1][slen] = ' ';
            s[2][slen] = '-';
            x++, slen++;
        }
        //Y insert
        else if (y < b_ptr->ypos)
        {
            s[0][slen] = '-';
            s[1][slen] = ' ';
            s[2][slen] = yc;
            y++, slen++;
        }
        //emit text batch
        else if (blen < b_ptr->len)
        {
            s[0][slen] = xc;
            s[1][slen] = (xc == yc)? '*' : ' ';
            s[2][slen] = yc;
            x++, y++, slen++, blen++;
        }
        else
            blen = 0, b_cnt--, b_ptr++;
        
        //print accumulated lines
        if ((slen > width - 7) || b_cnt <= 0)
        {
            //null terminate all strings
            for (int i = 0; i < 3; i++)
                s[i][slen] = 0;
            
            int xdisp = xseq.rev ? xseq.len - xstart - 1 : xstart;
            int ydisp = ystart;
            stream << endl << setw (margin) << "" << setw (6) << xdisp  << setw (0) << s[0];
            stream << endl << setw (margin) << "" << setw (6) << " "    << setw (0) << s[1];
            stream << endl << setw (margin) << "" << setw (6) << ydisp  << setw (0) << s[2];
            stream << endl;
            
            xstart = x, ystart = y, slen = 0;
        }
    }
        stream << flush;
}

//create text batch AAxAA
void print_batches_3 (AA_SEQ& xseq, AA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, ostream& stream, int margin, int width)
{
    int slen = 0;
    int blen = 0;
    int x = b_ptr->xpos;
    int y = b_ptr->xpos;
    int xstart = x;
    int ystart = y;
    int xc, yc;
    char s[3][256];

    stream << endl << setiosflags (ios::left);
    while (b_cnt > 0)
    {
        xc = xseq.get_code (x);
        yc = yseq.get_code (y);

        //special handling may be needed for (x < b_ptr->xpos && y < b_ptr->xpos)
        //X insert
        if (x < b_ptr->xpos)
        {
            sprintf (s[0] + slen, "%s", aa2str [xc]);
            sprintf (s[1] + slen, "   ");
            sprintf (s[2] + slen, "---");
            x++, slen += 3;
        }
        //Y insert
        else if (y < b_ptr->xpos)
        {
            sprintf (s[0] + slen, "---");
            sprintf (s[1] + slen, "   ");
            sprintf (s[2] + slen, "%s", aa2str [yc]);
            y++, slen += 3;
        }
        //emit text batch
        else if (blen < b_ptr->len)
        {
            sprintf (s[0] + slen, "%s", aa2str [xc]);
            sprintf (s[1] + slen, "%2d", w->mx[xc][yc]);
/*
            char wc = w->get_char (xc, yc);
            for (int i = 0; i < 3; i++)
                s[1][slen + i] = wc;
*/
            sprintf (s[2] + slen, "%s", aa2str [yc]);
            x++, y++, slen += 3, blen++;
        }
        else
            blen = 0, b_cnt--, b_ptr++;

        //print accumulated lines
        if ((slen > width - 7) || b_cnt <= 0)
        {
            //null terminate all strings
            for (int i = 0; i < 3; i++)
                s[i][slen] = 0;

            stream << endl << setw (margin) << "" << setw (6) << xstart << setw (0) << s [0];
            stream << endl << setw (margin) << "" << setw (6) << " "    << setw (0) << s [1];
            stream << endl << setw (margin) << "" << setw (6) << ystart << setw (0) << s [2];
            stream << endl;

            xstart = x, ystart = y, slen = 0;
        }
    }
    stream << flush;
}


//create text batch AAxNA
void print_batches (AA_SEQ& xseq, NA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, ostream& stream, int margin, int width)
{
    print_batches_an (xseq, yseq, b_ptr, b_cnt, w, stream, margin, width);
}


void print_batches_an (AA_SEQ& xseq, NA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, ostream& stream, int margin, int width)
{
    int slen = 0;
    int blen = 0;
    int x = b_ptr->ypos;
    int y = b_ptr->xpos;
    int xstart = x;
    int ystart = y;
    int xc, yc;
    char s[5][256];

    stream << endl << setiosflags (ios::left);
    while (b_cnt > 0)
    {
        xc = xseq.get_code (x);
        yc = yseq.get_code (y);

        //special handling may be needed for (x < b_ptr->ypos && y < b_ptr->xpos)
        //X (AA) insert
        if (x < b_ptr->ypos)
        {
            sprintf (s[0] + slen, "%s", aa2str [xc]);
            sprintf (s[1] + slen, "   ");
            sprintf (s[2] + slen, "---");
            sprintf (s[3] + slen, "---");
            sprintf (s[4] + slen, "   ");
            x++, slen += 3;
        }
        //Y (NA) insert 3x
        else if (y < b_ptr->xpos - 3)
        {
            sprintf (s[0] + slen, "---");
            sprintf (s[1] + slen, "   ");
            sprintf (s[2] + slen, "%s", aa2str [yc]);
            for (int i = 0; i < 3; i++)
                s[3][slen + i] = nn2char[yseq.get_base (y + i)];
            sprintf (s[4] + slen, "   ");
            y += 3, slen += 3;
        }
        //Y (NA) insert 1x
        else if (y < b_ptr->xpos)
        {
            s[0][slen] = '-';
            s[1][slen] = ' ';
            s[2][slen] = '-';
            s[3][slen] = nn2char[yseq.get_base (y)];
            s[4][slen] = ' ';
            y++, slen++;
        }
        //emit text batch
        else if (blen < b_ptr->len)
        {
            sprintf (s[0] + slen, "%s", aa2str [xc]);


            sprintf (s[1] + slen, "%3d", w->mx[xc][yc]);

            char wc = w->get_char (xc, yc);
            for (int i = 0; i < 3; i++)
            {
                // s[1][slen + i] = wc;
                s[3][slen + i] = nn2char[yseq.get_base (y + i)];
            }
            sprintf (s[2] + slen, "%s", aa2str [yc]);

            if (blen == 0)
                sprintf (s[4] + slen, "^f%1d", y % 3);
            else
                sprintf (s[4] + slen, "   ");

            blen++, slen += 3, x++, y += 3;
        }
        //advance to next batch
        else
            blen = 0, b_cnt--, b_ptr++;

        //print accumulated lines
        if ((slen > width - 7) || b_cnt <= 0)
        {
            //null terminate all strings
            for (int i = 0; i < 5; i++)
                s[i][slen] = 0;

            int xdisp = xseq.rev ? xseq.len - xstart - 1 : xstart;
            int ydisp = yseq.rev ? yseq.len - ystart - 1 : ystart;

            stream << endl << setw (margin) << "" << setw (6) << xdisp  << setw (0) << s [0];
            stream << endl << setw (margin) << "" << setw (6) << " "    << setw (0) << s [1];
            stream << endl << setw (margin) << "" << setw (6) << " "    << setw (0) << s [2];
            stream << endl << setw (margin) << "" << setw (6) << ydisp  << setw (0) << s [3];
            stream << endl << setw (margin) << "" << setw (6) << " "    << setw (0) << s [4];
            stream << endl;

            xstart = x, ystart = y, slen = 0;
        }
    }
    stream << flush;
}

unsigned decout (unsigned num, char* dest, unsigned destlen) // returns 0 on failure, number of chars on success
{
    unsigned pp = num;
    unsigned decposno = 1;
    while ((pp /= 10)) decposno ++;
    unsigned rval = decposno;
    if (decposno >= destlen) return 0;
    do
    {
        dest [--decposno] = char ('0' + num % 10);
    }
    while ((num /= 10));
    return rval;
}

int print_batches_cigar (const BATCH *b_ptr, int b_cnt, char* dest, unsigned destlen)
{
    unsigned dpos = 0, pl;
    int curb = 0;
    const BATCH *pb = NULL;
    assert (destlen > 1);
    while (curb <= b_cnt)
    {
        if (pb)
        {
            pl = decout (pb->len, dest + dpos, destlen - dpos);
            if (!pl)
                break;
            dpos += pl;
            if (dpos == destlen - 1)
                break;
            dest [dpos++] = 'M';
            if (curb < b_cnt)
            {
                if (pb->xpos + pb->len < b_ptr->xpos) // skip on x (subject) == gap on y (query)
                {
                    pl = decout (b_ptr->xpos - (pb->xpos + pb->len), dest + dpos, destlen - dpos);
                    if (!pl)
                        break;
                    dpos += pl;
                    if (dpos == destlen - 1)
                        break;
                    //dest [dpos++] = b_ptr->neutral_gap ? 'N' : 'I';
                    dest [dpos++] = 'I';
                }
                if (pb->ypos + pb->len < b_ptr->ypos) // skip on y (query) == gap on x (subject)
                {
                    pl = decout (b_ptr->ypos - (pb->ypos + pb->len), dest + dpos, destlen - dpos);
                    if (!pl)
                        break;
                    dpos += pl;
                    if (dpos == destlen - 1)
                        break;
                    dest [dpos++] = 'D';
                }
            }
        }
        pb = b_ptr;
        b_ptr ++;
        curb ++;
    }
    dest [dpos] = 0;
    return dpos;
}
