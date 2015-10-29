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

#ifndef __biosequence_h__
#define __biosequence_h__

#include <ctype.h>
#include <vector>
#include <platform.h>
#include "sequtil.h"

//common sequence translation helpers


extern const char* const aa2str [AANUM];
extern const char* const gcode_c;
extern const char* const gcode_c_rev;
extern const int gcode_b [64];
extern const int gcode_b_rev [64];



//---------------------------------------------------------------------------------
//sequence base class
struct SEQ
{
    char*           seq;
    int             len: 31;
    int             rev: 1;
    longlong        uid;
    unsigned char   kt_size:7; //kt_size variable holds ktuple mask for NN sequences
    mutable bool    owner: 1;

    SEQ ()
    :
    seq (NULL),
    len (0),
    rev (0),
    uid (0),
    kt_size (0),
    owner (false)
    {
    }
    SEQ (const SEQ& tc)
    :
    seq (tc.seq),
    len (tc.len),
    rev (tc.rev),
    uid (tc.uid),
    kt_size (tc.kt_size),
    owner (false)
    {
        if (tc.owner)
        {
            tc.owner = false;
            owner = true;
        }
    }
    virtual ~SEQ ()
    {
        clean ();
    }
    SEQ& operator = (const SEQ& other)
    {
       seq = other.seq;
       len = other.len;
       rev = other.rev;
       uid = other.uid;
       kt_size = other.kt_size;
       if (owner = other.owner) other.owner = false;
       return *this;
    }
    virtual void clean ()
    {
        if (owner && seq)
        {
            delete seq;
            seq = NULL;
        }
    }
    virtual int      bases_per_int () = 0;
    virtual char     code2char (int code) = 0;

    virtual unsigned get_bases (int pos) = 0;
    virtual unsigned get_base (int pos) = 0;
    virtual void     set_base (int pos, int base) = 0;

    virtual unsigned set_kt_size (int size) = 0;
    virtual unsigned get_ktup (int pos) = 0;

    virtual unsigned get_code (int pos) = 0;
};


//---------------------------------------------------------------------------------
//bit packed binary nucleotide sequence
//byte packing - MSB 33221100 LSB
struct NN_SEQ:SEQ
{
    NN_SEQ& operator = (const NN_SEQ& other)
    {
        return (NN_SEQ&) SEQ::operator = (other);
    }
    int         bases_per_int () {return 12;}
    char        code2char (int code) {return number2base [code];}

    unsigned    get_bases (int pos);
    unsigned    get_base (int pos);
    void        set_base (int pos, int base);

    unsigned    set_kt_size (int size) {kt_size = ~(0xffffffff << (size << 1)); return (NNNUM ^ size);}
    unsigned    get_ktup (int pos);

    unsigned    get_code (int pos) {return get_base (pos);}
};

//get 12 bases, upper bits undefined
inline unsigned NN_SEQ::get_bases (int pos)
{
    return GET32_U (seq + (pos >> 2)) >> ((pos & 3) << 1);
}

//get single base
inline unsigned NN_SEQ::get_base (int pos)
{
    return (seq[pos >> 2] >> ((pos & 3) << 1)) & 3;
}

//put single base
inline void NN_SEQ::set_base (int pos, int base)
{
    seq [pos >> 2] &= ~(0x3 << ((pos & 3) << 1)); // zero the two bits
    seq [pos >> 2] |= (base << ((pos & 3) << 1)); // set the two bits
}

//get k-tuple
inline unsigned NN_SEQ::get_ktup (int pos)
{
    return (GET32_U (seq + (pos >> 2)) >> ((pos & 3) << 1)) & kt_size;  //kt_size variable holds ktuple mask
}

class SEQ_compare_to_id
{
public:
    static bool cmp (const SEQ& c1, longlong uid)
    {
        return (c1.uid < uid);
    }
};

class SEQ_compare
{
public:
    static bool cmp (const SEQ& c1, const SEQ& c2)
    {
        if (c1.uid < c2.uid) return true;
        else if (c1.uid > c2.uid) return false;
        else return (c1.rev <= c2.rev);
    }
};

//----------------------------------------------------------------------------------

//byte packed binary aminoacid sequence
//0 < *seq < 24
struct AA_SEQ:SEQ
{
    AA_SEQ& operator = (const AA_SEQ& other)
    {
        return (AA_SEQ&) SEQ::operator = (other);
    }
    int         bases_per_int () {return 4;}
    char        code2char (int code) {return aa2char (code);}

    unsigned    get_bases (int pos);
    unsigned    get_base (int pos);
    void        set_base (int pos, int base);

    unsigned    set_kt_size (int size) {kt_size = size; return (AANUM ^ size);}
    unsigned    get_ktup (int pos);

    unsigned    get_code (int pos) {return get_base (pos);}
};


//get 4 bases
inline unsigned AA_SEQ::get_bases (int pos)
{
    return GET32_U (seq + pos);
}

//get single base
inline unsigned AA_SEQ::get_base (int pos)
{
    return (unsigned) seq[pos];
}

//put single base
inline void AA_SEQ::set_base (int pos, int base)
{
    seq[pos] = base;
}

//get k-tuple
inline unsigned AA_SEQ::get_ktup (int pos)
{
    unsigned r = seq[pos];
    for (int i = 1; i < (int)kt_size; i++)
        r *= (unsigned) seq[++pos];

    return r;
}

//----------------------------------------------------------------------------------
//aminoacid sequence represented by
//bit packed binary nucleotide sequence
//byte packing - MSB 33221100 LSB
struct NA_SEQ:SEQ
{
    NA_SEQ& operator = (const AA_SEQ& other)
    {
        return (NA_SEQ&) SEQ::operator = (other);
    }
    int         bases_per_int () {return 12;}
    char        code2char (int code) {return number2base [code];}

    unsigned    get_bases (int pos);
    unsigned    get_base (int pos);
    void        set_base (int pos, int base);

    unsigned    set_kt_size (int size) {kt_size = size; return (AANUM ^ size);}
    unsigned    get_ktup (int pos);

    unsigned    get_code (int pos) {return gcode_b [get_bases (pos) & 0x3f];}
};

//get 12 bases, upper bits undefined
inline unsigned NA_SEQ::get_bases (int pos)
{
    return GET32_U (seq + (pos >> 2)) >> ((pos & 3) << 1);
}

//get single base
inline unsigned NA_SEQ::get_base (int pos)
{
    return (seq[pos >> 2] >> ((pos & 3) << 1)) & 3;
}

//put single base
inline void NA_SEQ::set_base (int pos, int base)
{
    // assert (!(base & ~0x3));
    seq [pos >> 2] &= ~(0x3 << ((pos & 3) << 1)); // zero the two bits
    seq [pos >> 2] |= (base << ((pos & 3) << 1)); // set the two bits
}

//get k-tuple
inline unsigned NA_SEQ::get_ktup (int pos)
{
    unsigned r = get_code (pos);
    for (int i = 1; i < (int)kt_size; i++)
        pos += 3, r *= get_code(pos);

    return r;
}

void reverse_seqs (std::vector <NN_SEQ>& src, std::vector <NN_SEQ>& dest);
void reverse_seqs (std::vector <NA_SEQ>& src, std::vector <NA_SEQ>& dest);

#endif // __biosequence_h__
