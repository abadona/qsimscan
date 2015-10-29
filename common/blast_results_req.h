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

#ifndef __blast_results_req_h__
#define __blast_results_req_h__

#include <platform.h>
#include "biosequence.h"

struct req
{
    longlong uid;
    char rev;
    req ()
    {
        uid = -1;
        rev = 0;
    }
    req (const NN_SEQ& seq)
    {
        uid = seq.uid;
        rev = seq.rev;
    }
    bool operator < (const req& other) const
    {
        if (uid < other.uid) return true;
        else if (uid > other.uid) return false;
        else if (rev < other.rev) return true;
        else return false;
    }
    bool operator == (const req& other) const
    {
        return ((uid == other.uid) && (rev == other.rev));
    }
    bool operator == (const NN_SEQ& seq) const
    {
        return ((uid == seq.uid) && (rev == seq.rev));
    }
    void operator = (NN_SEQ& seq)
    {
        uid = seq.uid;
        rev = seq.rev;
    }
};

#endif // __blast_results_req_h__
