
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
// For any questions please contact SciDM team by email at team@scidm.org
//////////////////////////////////////////////////////////////////////////////

#ifndef __fasta_h__
#define __fasta_h__

#include <cstdio>
#include "platform.h"

#define MAX_HDR_LEN 8000
#define MAX_NAME_LEN 160
#define MAX_LINE_LEN 10000
// #define INIT_SEQ_LEN 300000000ULL // 300M
#define INIT_SEQ_LEN 1000000ULL // 1M

class FastaFile
{
    FILE* f_;
    ulonglong tot_len_;
    ulonglong cur_pos_;
    ulonglong prev_pos_;
    ulonglong cur_recstart_;
    unsigned cur_reclen_;
    unsigned seq_no_;
    char *seqbuf_;
    unsigned seq_buf_sz_;
    unsigned seqlen_;
    char hdrbuf_  [MAX_HDR_LEN+1];
    char namebuf_ [MAX_NAME_LEN+1];
    char linebuf_ [MAX_LINE_LEN];
    char* l_;

    void reset ();
    void parse_hdr ();
    void add_seq ();

public:
    FastaFile (ulonglong init_sz = INIT_SEQ_LEN);
    FastaFile (const char* name, ulonglong init_sz = INIT_SEQ_LEN);
    ~FastaFile ();

    bool        open (const char* name);
    bool        close ();
    bool        next ();
    bool        seek (ulonglong off);
    bool        is_open () const {return f_ != NULL;}
    void        fetch_hdr () { parse_hdr (); }

    const char* cur_name () const;
    const char* cur_hdr  () const;
    const char* cur_seq  () const;
    char*       cur_seq_buf ();
    unsigned    cur_seq_len () const { return seqlen_; }
    unsigned    cur_no  () const;
    ulonglong   tot_len () const;
    ulonglong   cur_pos () const;
    unsigned    cur_reclen () const;
    ulonglong   cur_recstart () const;
};

#endif
