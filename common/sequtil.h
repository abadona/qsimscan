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

#ifndef __sequtil_h__
#define __sequtil_h__

#include <ctype.h>
#include <ostream>
#include <platform.h>
#include <tracer.h>
#include <bitops.h>
#include <myassert.h>
#include "align_batch.h"


// definitions related to binary sequence encoding
#define NNNUM 4
#define AANUM 24

#define BITS_PER_BASE 2
#define BITS_PER_BASE_SHIFT 1
#define BASES_PER_BYTE (BITS_PER_BYTE/BITS_PER_BASE)
#define BASES_PER_BYTE_SHIFT 2 // how much the length of sequence in bases should be shifted right to get number of bytes the packed sequence will take

void n_revert_seq   (char* dest, const char* source, unsigned len, unsigned beg=0);
bool n_ascii2binary (char* dest, unsigned destlen, const char* source, unsigned start, unsigned len);
bool n_binary2ascii (char* dest, unsigned destlen, const char* source, unsigned start, unsigned len);
void n_countnucs    (unsigned& a, unsigned& g, unsigned& c, unsigned& t, const char* src, unsigned start = 0, unsigned len = -1);
std::ostream& output_seq (std::ostream& o, const char* buffer, unsigned buflen, unsigned chars_in_line, bool counts, bool decades, unsigned first_off);
std::ostream& fasta_output_seq (std::ostream& o, const char* buffer, unsigned buflen, unsigned chars_per_line);
std::ostream& gb_output_seq (std::ostream& o, const char* buffer, unsigned buflen);
unsigned swap_batches (BATCH* batches, unsigned batch_no); // returns total batch length
unsigned len_batches (const BATCH* batches, unsigned batch_no);
unsigned count_gaps (const BATCH* batches, unsigned batch_no);
unsigned count_matches (const char* xeq, const char* yseq, const BATCH* batches, unsigned batch_no);
bool a_ascii2binary (char* dest, unsigned destlen, const char* source, unsigned start, unsigned len);
bool a_binary2ascii (char* dest, unsigned destlen, const char* source, unsigned start, unsigned len);

// get single base
inline unsigned get_base (const char *seq, unsigned pos)
{
    return (seq[pos >> 2] >> ((pos & 3) << 1)) & 3;
}

// put single base
inline void put_base (char *seq, unsigned pos, unsigned base)
{
    // assert (!(base & ~0x3));
    seq [pos >> 2] &= ~(0x3 << ((pos & 3) << 1)); // zero the two bits
    seq [pos >> 2] |= (base << ((pos & 3) << 1)); // set the two bits
}

// get k-tuple, 'short' version. Extraction of tuples not longer then 13 bp
// in k-tuple, the earlier bits appear as less significant
inline DWORD get_ktup (const char *seq, unsigned pos, DWORD kt_mask)
{
    return (GET32_U (seq + (pos >> 2)) >> ((pos & 3) << 1)) & kt_mask;
}

// get k-tuple, 'long' version. Extraction of tuples not longer then 29 bp
// in k-tuple, the earlier bits appear as less significant
inline QWORD get_lktup (const char *seq, unsigned pos, QWORD lkt_mask)
{
    return (GET64_U (seq + (pos >> 2)) >> ((pos & 3) << 1)) & lkt_mask;
}

extern BYTE LTUPLE_SHIFTS [];
extern const BYTE* const TUPLE_SHIFTS;

// this could be more efficient if templetized and computed at compile time. 
// However this will mean miltiplication (bloating) of same code for all possible ktuples

inline BYTE tuple_shift (BYTE len)
{
    // return (((sizeof (DWORD) << BASES_PER_BYTE_SHIFT) - len) << BITS_PER_BASE_SHIFT);
    return TUPLE_SHIFTS [len];
}

inline BYTE ltuple_shift (BYTE len)
{
    // return (((sizeof (QWORD) << BASES_PER_BYTE_SHIFT) - len) << BITS_PER_BASE_SHIFT);
    return LTUPLE_SHIFTS [len];
}

extern QWORD LTUPLE_MASKS [];
extern DWORD TUPLE_MASKS [];


// compute mask for a tuple of given length, 'long' version - tuples below 32 bases
inline QWORD ltuple_mask (BYTE len)
{
    // return (~((QWORD) 0)) >> (sizeof (QWORD)*BITS_PER_BYTE - len * BITS_PER_BASE);
    return LTUPLE_MASKS [len];
}

// compute mask for a tuple of given length, 'long' version - tuples below 32 bases
inline DWORD tuple_mask (BYTE len)
{
    // return (~((QWORD) 0)) >> (sizeof (QWORD)*BITS_PER_BYTE - len * BITS_PER_BASE);
    return TUPLE_MASKS [len];
}

inline BYTE swapbases_hard (BYTE in)
{
    BYTE out = 0;
    out |= (in & 0x03) << 6;
    out |= (in & 0x0C) << 2;
    out |= (in & 0x30) >> 2;
    out |= (in & 0xC0) >> 6;
    return out;
}

static const unsigned SWAP_TABLE_SIZE = 0x100;
static BYTE SWAP_TABLE [SWAP_TABLE_SIZE];
static bool make_swap_table ()
{
    for (unsigned a = 0; a != SWAP_TABLE_SIZE; ++a)
    {
        SWAP_TABLE [a] = swapbases_hard (a);
    }
    return true;
}
static bool table_made = make_swap_table ();

inline BYTE swapbases (BYTE in)
{
    return SWAP_TABLE [in];
}

#ifdef SCIDM_LITTLE_ENDIAN

#define GET16BASES(ptr) ( (((DWORD) swapbases (((const BYTE*) ptr) [0])) << 24) | \
                          (((DWORD) swapbases (((const BYTE*) ptr) [1])) << 16) | \
                          (((DWORD) swapbases (((const BYTE*) ptr) [2])) << 8 ) | \
                          (((DWORD) swapbases (((const BYTE*) ptr) [3]))      ) )

#define GET32BASES(ptr) ( (((QWORD) swapbases (((const BYTE*) ptr) [0])) << 56) | \
                          (((QWORD) swapbases (((const BYTE*) ptr) [1])) << 48) | \
                          (((QWORD) swapbases (((const BYTE*) ptr) [2])) << 40) | \
                          (((QWORD) swapbases (((const BYTE*) ptr) [3])) << 32) | \
                          (((QWORD) swapbases (((const BYTE*) ptr) [4])) << 24) | \
                          (((QWORD) swapbases (((const BYTE*) ptr) [5])) << 16) | \
                          (((QWORD) swapbases (((const BYTE*) ptr) [6])) << 8 ) | \
                          (((QWORD) swapbases (((const BYTE*) ptr) [7]))      ) )

#else

#define GET16BASES(ptr) ( (((DWORD) ((const BYTE*) ptr) [3]) << 24) | \
                          (((DWORD) ((const BYTE*) ptr) [2]) << 16) | \
                          (((DWORD) ((const BYTE*) ptr) [1]) << 8 ) | \
                          (((DWORD) ((const BYTE*) ptr) [0])      ) )

#define GET32BASES(ptr) ( (((QWORD) ((const BYTE*) ptr) [7]) << 56) | \
                          (((QWORD) ((const BYTE*) ptr) [6]) << 48) | \
                          (((QWORD) ((const BYTE*) ptr) [5]) << 40) | \
                          (((QWORD) ((const BYTE*) ptr) [4]) << 32) | \
                          (((QWORD) ((const BYTE*) ptr) [3]) << 24) | \
                          (((QWORD) ((const BYTE*) ptr) [2]) << 16) | \
                          (((QWORD) ((const BYTE*) ptr) [1]) << 8 ) | \
                          (((QWORD) ((const BYTE*) ptr) [0])      ) )


#endif


// get tuple, 'short' version. Extraction of tuples not longer then 13 bp
// in tuple, the earlier bits appear as more significant
inline DWORD get_tup (const char *seq, unsigned pos, BYTE len)
{
    return ((GET16BASES (seq + (pos >> 2)) << ((pos & 3) << 1)) >> tuple_shift (len)) & tuple_mask (len);
}

// get tuple, 'long' version. Extraction of tuples not longer then 29 bp
// in tuple, the earlier bits appear as more significant
inline QWORD get_ltup (const char *seq, unsigned pos, BYTE len)
{
/*
    QWORD t = GET32BASES (seq + (pos >> 2));
    // t is bitpair-inverse
    BYTE sh = ((pos & 3) << 1);
    t <<= sh;
    BYTE tsh = ltuple_shift (len);
    t >>= tsh;
    QWORD mask = ltuple_mask (len);
    t &= mask;
    return t;
*/
    return ((GET32BASES (seq + (pos >> 2)) << ((pos & 3) << 1)) >> ltuple_shift (len)) & ltuple_mask (len);
}

// get 12 bases, upper bits undefined
inline unsigned get_12_bases (const char *seq, int pos) // DO NOT CHANGE SECOND ARG TO unsigned type! This will break the nsimscan code!
{
    // HACK / IMPROPER! this code uses features of negative number shifting/
    return GET32_U (seq + (pos >> 2)) >> ((pos & 3) << 1);
}

// special case of bit counting algorithm
inline unsigned  count_12 (unsigned r)
{
    r += r >> 2;        //16 x 2-bit vector addition
    r &= 0x333333;
    r += r >> 4;          // 8 x 4-bit vector addition
    r += r >> 8;          // 4 x 4-bit vector addition
    r += r >> 16;        // 2 x 4-bit vector addition
    return (r & 0xf);
}

extern const char aa2number [];
extern const char number2aa [];
extern const char number2base [];

inline char base2char (unsigned b)
{
    if (b > 4)
        b = 0;
    return number2base [b];
}

inline unsigned char2base (char c)
{
    switch (tolower(c))
    {
        case 'r':
        case 'w':
        case 'g':
            return 1;

        case 'y':
        case 'c':
            return 2;

        case 'k':
        case 'u':
        case 't':
            return 3;
        default:
            return 0;
    }
}

inline unsigned char2aa (char c)
{
    if (c >= 'a' && c <= 'z') c -= ('a' - 'A');
    if (c <'A' || c > 'Z') return 23;
    return aa2number [c - 'A'];
}

inline char aa2char (unsigned aanum)
{
    if (aanum < 0 || aanum > 23) return 'X';
    return number2aa [aanum];
}

inline QWORD get_tuple (const char* seq, unsigned off, unsigned len, DWORD* mask)
{
    register QWORD rv = 0;
    register unsigned char c;
    if (mask) *mask = 0UL;
    seq += off;
    myassert (len <= 32);
    while (len)
    {
        if (mask) *mask <<= 1;
        c = (unsigned char) char2base (*seq);
        if (c == 0xff)
        {
            if (mask) *mask |= 1;
            c = 0;
        }
        rv <<= 2;
        rv |= c;
        seq ++;
        len --;
    }
    return rv;
}

inline void tuple2ascii (QWORD tuple, unsigned len, char* dest) // omitted  destlen !!!
{
    dest [len] = 0;
    while (len)
    {
        dest [len - 1] = base2char (tuple & 3);
        tuple >>= 2;
        len --;
    }
}

// helper class for printing out tuples
struct TUPLE
{
    public:
        QWORD tuple_;
        BYTE l_;
        TUPLE (QWORD t, BYTE l)
        :
        tuple_ (t),
        l_ (l)
        {
        }
        TUPLE ()
        :
        l_ (0)
        {
        }
};

inline std::ostream& operator << (std::ostream& ostr, const TUPLE& tuple)
{
    char buf [33];
    tuple2ascii (tuple.tuple_, tuple.l_, buf);
    ostr << buf;
    return ostr;
}

inline Logger& operator << (Logger& logger, const TUPLE& tuple)
{
    if (logger.enabled ())
        logger.o_ << tuple;
    return logger;
}



#endif // __sequtil_h__
