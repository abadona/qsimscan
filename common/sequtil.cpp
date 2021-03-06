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

#define __sequtil_cpp__

#include <cstring>
#include <iomanip>
#include <ios>
#include "sequtil.h"
#include "common_str.h"

void n_revert_seq (char* dest, const char* src, unsigned len, unsigned beg)
{
    for (int ii = 0; ii < len; ii ++)
        put_base (dest, len - ii - 1, ((~get_base (src, ii+beg)) & 0x3));
}

bool n_ascii2binary (char* dest, unsigned destlen, const char* source, unsigned start, unsigned len)
{
    if (destlen < ((len + 3) >> 2))
        return false;

    for (int i = start; i < start + len; i ++)
        put_base (dest, i - start, char2base (source [i]));

    return true;
}

bool n_binary2ascii (char* dest, unsigned destlen, const char* source, unsigned start, unsigned len)
{
    if (destlen < len)
        return false;

    for (unsigned i = start; i < start + len; i ++)
        dest [i - start] = base2char (get_base (source, i));
    // only zero-terminate if there is enough space in dest
    if (destlen > len)
        dest [len] = 0;

    return true;
}

void n_countnucs (int& a, unsigned & g, unsigned & c, unsigned & t, const char* src, unsigned start, unsigned len)
{
    a = g = t = c = 0;
    if (len == -1) len = strlen (src) - start;

    for (unsigned i = start; i < start + len; i ++)
        switch (tolower (src [i]))
        {
        case 't':   t ++; break;
        case 'a':   a ++; break;
        case 'g':   g ++; break;
        case 'c':   c ++; break;
        }
}

std::ostream& output_seq (std::ostream& o, const char* buffer, unsigned buflen, unsigned chars_in_line, bool counts, bool decades, unsigned first_off)
{
    unsigned cnt = 0;

    while (cnt < buflen)
    {
        if (cnt % chars_in_line == 0)
        {
            if (cnt) o << std::endl;
            if (counts) o << std::setw (first_off - 1) << std::right << cnt + 1 << " " << std::setw (0);
            else o << std::setw (first_off) << /*std::right <<*/ EMPTY_STR << std::setw (0);
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

std::ostream& fasta_output_seq (std::ostream& o, const char* buffer, unsigned buflen, unsigned chars_per_line)
{
    return output_seq (o, buffer, buflen, chars_per_line, false, false, 0);
}

std::ostream& gb_output_seq (std::ostream& o, const char* buffer, unsigned buflen)
{
    return output_seq (o, buffer, buflen, 60, true, true, 10);
}

                         // A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z   *
const char aa2number []  = {0,  23, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23, 14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23 };

                         // 0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23
const char number2aa [] = {'A','R','N','D','C','Q','E','G','H','I','L','K', 'M','F','P','S','T','W','Y','V',' ','Z','X','*' };

const char number2base [] = {'a', 'g', 'c' , 't'};


bool a_ascii2binary (char* dest, unsigned destlen, const char* source, unsigned start, unsigned len)
{
    if (destlen < len)
        return false;

    for (unsigned i = start; i < start + len; i ++)
        dest [i - start] = char2aa (source [i]);

    return true;
}

bool a_binary2ascii (char* dest, unsigned destlen, const char* source, unsigned start, unsigned len)
{
    if (destlen < len)
        return false;

    for (unsigned i = start; i < start + len; i ++)
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

unsigned swap_batches (BATCH* batches, unsigned batch_no)
{
    unsigned totlen = 0;
    unsigned xpos;
    BATCH* b = batches;
    for (unsigned idx = 0; idx < batch_no; idx ++)
    {
        xpos = b->xpos;
        b->xpos = b->ypos;
        b->ypos = xpos;
        totlen += b->len;
        b ++;
    }
    return totlen;
}

unsigned len_batches (const BATCH* batches, unsigned batch_no)
{
    int totlen = 0;
    for (int idx = 0; idx < batch_no; idx ++)
        totlen += batches [idx].len;
    return totlen;
}

unsigned count_matches (const char* xseq, const char* yseq, const BATCH* batches, unsigned batch_no)
{
    unsigned matches = 0;
    const BATCH* curb;
    unsigned xp, yp;
    unsigned pos;
    for (unsigned idx = 0; idx < batch_no; idx ++)
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
unsigned count_gaps (const BATCH* batches, unsigned batch_no)
{
    unsigned gaplen = 0;
    unsigned lastx, lasty;
    const BATCH* curb;
    for (unsigned idx = 0; idx < batch_no; idx ++)
    {
        curb = batches + idx;
        if (idx)
            gaplen += (curb->xpos - lastx) + (curb->ypos - lasty);
        lastx = curb->xpos + curb->len;
        lasty = curb->ypos + curb->len;
    }
    return gaplen;
}

QWORD LTUPLE_MASKS [] = {
    0x0000000000000000LL,
    0x0000000000000003LL,
    0x000000000000000FLL,
    0x000000000000003FLL,
    0x00000000000000FFLL,
    0x00000000000003FFLL,
    0x0000000000000FFFLL,
    0x0000000000003FFFLL,
    0x000000000000FFFFLL,
    0x000000000003FFFFLL,
    0x00000000000FFFFFLL,
    0x00000000003FFFFFLL,
    0x0000000000FFFFFFLL,
    0x0000000003FFFFFFLL,
    0x000000000FFFFFFFLL,
    0x000000003FFFFFFFLL,
    0x00000000FFFFFFFFLL,
    0x00000003FFFFFFFFLL,
    0x0000000FFFFFFFFFLL,
    0x0000003FFFFFFFFFLL,
    0x000000FFFFFFFFFFLL,
    0x000003FFFFFFFFFFLL,
    0x00000FFFFFFFFFFFLL,
    0x00003FFFFFFFFFFFLL,
    0x0000FFFFFFFFFFFFLL,
    0x0003FFFFFFFFFFFFLL,
    0x000FFFFFFFFFFFFFLL,
    0x003FFFFFFFFFFFFFLL,
    0x00FFFFFFFFFFFFFFLL,
    0x03FFFFFFFFFFFFFFLL,
    0x0FFFFFFFFFFFFFFFLL,
    0x3FFFFFFFFFFFFFFFLL,
    0xFFFFFFFFFFFFFFFFLL
};

DWORD TUPLE_MASKS [] = {
    0x00000000,
    0x00000003,
    0x0000000F,
    0x0000003F,
    0x000000FF,
    0x000003FF,
    0x00000FFF,
    0x00003FFF,
    0x0000FFFF,
    0x0003FFFF,
    0x000FFFFF,
    0x003FFFFF,
    0x00FFFFFF,
    0x03FFFFFF,
    0x0FFFFFFF,
    0x3FFFFFFF,
    0xFFFFFFFF,
};

BYTE LTUPLE_SHIFTS [] = {
    64,
    62,
    60,
    58,
    56,
    54,
    52,
    50,
    48,
    46,
    44,
    42,
    40,
    38,
    36,
    34,
    32,
    30,
    28,
    26,
    24,
    22,
    20,
    18,
    16,
    14,
    12,
    10,
    8,
    6,
    4,
    2,
    0
};

const BYTE* const TUPLE_SHIFTS = LTUPLE_SHIFTS + 16;


