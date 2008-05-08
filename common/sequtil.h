
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

#ifndef __sequtil_h__
#define __sequtil_h__

#include <ctype.h>
#include <ostream>
#include <platform.h>
#include "align_batch.h"


void n_revert_seq   (char* dest, const char* source, int len, int beg=0);
bool n_ascii2binary (char* dest, int destlen, const char* source, int start, int len);
bool n_binary2ascii (char* dest, int destlen, const char* source, int start, int len);
void n_countnucs    (int& a, int& g, int& c, int& t, const char* src, int start = 0, int len = -1);
std::ostream& output_seq (std::ostream& o, const char* buffer, int buflen, int chars_in_line, bool counts, bool decades, int first_off);
std::ostream& fasta_output_seq (std::ostream& o, const char* buffer, int buflen, int chars_per_line);
std::ostream& gb_output_seq (std::ostream& o, const char* buffer, int buflen);
int swap_batches (BATCH* batches, int batch_no); // returns total batch length
int len_batches (const BATCH* batches, int batch_no);
int count_gaps (const BATCH* batches, int batch_no);
int count_matches (const char* xeq, const char* yseq, const BATCH* batches, int batch_no);
bool a_ascii2binary (char* dest, int destlen, const char* source, int start, int len);
bool a_binary2ascii (char* dest, int destlen, const char* source, int start, int len);

//get single base
inline unsigned get_base (const char *seq, int pos)
{
    return (seq[pos >> 2] >> ((pos & 3) << 1)) & 3;
}

//put single base
inline void put_base (char *seq, int pos, int base)
{
    // assert (!(base & ~0x3));
    seq [pos >> 2] &= ~(0x3 << ((pos & 3) << 1)); // zero the two bits
    seq [pos >> 2] |= (base << ((pos & 3) << 1)); // set the two bits
}

//get k-tuple
inline unsigned get_ktup (char *seq, int pos, unsigned kt_mask)
{
    return (GET32_U (seq + (pos >> 2)) >> ((pos & 3) << 1)) & kt_mask;
}

//get 12 bases, upper bits undefined
inline unsigned get_12_bases (char *seq, int pos)
{
    return GET32_U (seq + (pos >> 2)) >> ((pos & 3) << 1);
}


//special case of bit counting algorithm
inline int count_12 (unsigned r)
{
        r += r >> 2;		//16 x 2-bit vector addition
        r &= 0x333333;
        r += r >> 4;  		// 8 x 4-bit vector addition
        r += r >> 8;  		// 4 x 4-bit vector addition
        r += r >> 16;		// 2 x 4-bit vector addition

        return (r & 0xf);
}

inline char base2char (int b)
{
    switch (b)
    {
        case 1:	return 'g';
        case 2:	return 'c';
        case 3:	return 't';
        default: return 'a';
    }
}

inline int char2base (char c)
{
    switch (toupper(c))
    {
        case 'R':
        case 'W':
        case 'G':
            return 1;

        case 'Y':
        case 'C':
            return 2;

        case 'K':
        case 'U':
        case 'T':
            return 3;
        default:
            return 0;
    }
}

extern const char aa2number [];
extern const char number2aa [];

inline int char2aa (char c)
{
    if (c >= 'a' && c <= 'z') c -= ('a' - 'A');
    if (c <'A' || c > 'Z') return 23;
    return aa2number [c - 'A'];
}

inline char aa2char (char aanum)
{
    if (aanum < 0 || aanum > 23) return 'X';
    return number2aa [aanum];
}


#endif // __sequtil_h__
