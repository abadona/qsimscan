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

#include <cstring>
#include <ctype.h>
#include <cstdlib>
#include <cerrno>

#include "fasta.h"
#include <rerror.h>
#include <tracer.h>
#include <myassert.h>
#include <fileutils.h>
#include <portability.h>

static const char file_not_open [] = "INTERNAL ERROR: Operating on unopen fasta file";

void FastaFile::reset ()
{
    cur_pos_ = 0;
    prev_pos_ = 0;
    *linebuf_ = 0;
    *seqbuf_ = 0;
    *namebuf_ = 0;
    l_ = linebuf_;
}


void FastaFile::parse_hdr ()
{
    if (!*l_) return;
    myassert (*l_ == '>');

    char* p = l_ + 1;

    // skip to first non-space
    while (*p && isspace (*p)) p ++;
    // if at the end of the line, err
    if (!*p) ERR ("Empty FASTA header");

    char* p1 = p;
    // skip to first space
    while (*p1 && !isspace (*p1) && p1 - p < MAX_NAME_LEN) p1 ++;
    // copy name into namebuf
    memcpy (namebuf_, p, p1 - p);
    namebuf_ [p1 - p] = 0;

    // skip to non-space
    p = p1;
    while (*p && isspace (*p)) p ++;
    // if not at the end, copy into hdrbuf
    if (*p)
    {
        strncpy (hdrbuf_, p, MAX_HDR_LEN-1);
        hdrbuf_ [MAX_HDR_LEN-1] = 0;
        // trim end spaces / newlines
        p = hdrbuf_ + strlen (hdrbuf_) - 1;
        while (p >= hdrbuf_ && isspace (*p))
        {
            *p = 0;
            p --;
        }
    }
    else
        *hdrbuf_ = 0;
}

void FastaFile::add_seq ()
{
    char* p = l_;
    while (*p)
    {
        if (!isspace (*p))
        {
            if (seqlen_ == seq_buf_sz_)
            {
                warnlog << "Reallocating FastaReader buffer!" << std::endl;
                // reallocate seqbuf_: increment space twice
                unsigned new_sz = seq_buf_sz_ * 2;
                char* newbuf;
                try
                {
                    newbuf = new char [new_sz + 1];
                }
                catch (std::bad_alloc&)
                {
                    newbuf = NULL;
                }
                if (!newbuf)
                    ers << "ERROR: Not enough memory to hold the sequence" << namebuf_ << ": error allocating " << new_sz + 1 << " bytes" <<  Throw;
                memcpy (newbuf, seqbuf_, seq_buf_sz_);
                delete [] seqbuf_;
                seqbuf_ = newbuf;
                seq_buf_sz_ = new_sz;
            }
            seqbuf_ [seqlen_] = tolower (*p);
            seqlen_ ++;
        }
        p ++;
    }
}

FastaFile::FastaFile (ulonglong init_sz)
:
f_ (NULL)
{
    try
    {
        seqbuf_ = new char [init_sz+1];
    }
    catch (std::bad_alloc&)
    {
        seqbuf_ = NULL;
    }
    if (!seqbuf_)
        ers << "Memory allocation error: unable to allocate " << init_sz + 1 << " bytes" << Throw;

    seq_buf_sz_ = init_sz;
    reset ();
}

FastaFile::FastaFile (const char* name, ulonglong init_sz)
:
f_ (NULL)
{
    try
    {
        seqbuf_ = new char [init_sz+1];
    }
    catch (std::bad_alloc&)
    {
        seqbuf_ = NULL;
    }
    if (!seqbuf_)
        ers << "Memory allocation error: unable to allocate " << init_sz + 1 << " bytes"  << Throw;

    seq_buf_sz_ = init_sz;
    reset ();
    if (!open (name))
    {
        ers << "ERROR: Could not open fasta file" << name << Throw;
    }
}

FastaFile::~FastaFile ()
{
    close ();
    delete [] seqbuf_;
    seq_buf_sz_ = 0;
}

bool FastaFile::open (const char* name)
{
    if (f_) fclose (f_);
    reset ();

    l_ = linebuf_;
    *l_ = 0;
    f_ = fopen (name, "rb");
    if (f_ != NULL)
    {
        fseek (f_, 0, SEEK_END);
        tot_len_ = ftell (f_);
        fseek (f_, 0, SEEK_SET);
        while (*l_ != '>')
        {
            prev_pos_ = cur_pos_;
            if (!(l_ = fgets (linebuf_, MAX_LINE_LEN, f_)))
            {
                if (ferror (f_))
                {
                    ers << "Error reading fasta file: " << strerror (errno) << Throw;
                }
                else
                {
                    fclose (f_);
                    f_ = NULL;
                    *linebuf_ = 0;
                    ers << "Fasta file: " << name << " has no valid records" << Throw;
                }
            }
            cur_pos_ = ftell (f_);
        }
        seq_no_ = unsigned (-1);
        return true;
    }
    return false;
}

bool FastaFile::close ()
{
    if (f_)
    {
        fclose  (f_);
        f_ = NULL;
        return true;
    }
    else
        return false;
}

bool FastaFile::next ()
{
    if (cur_pos_ == tot_len_) return false;
    cur_recstart_ = prev_pos_;
    parse_hdr ();
    seqlen_ = 0;
    l_ = fgets (linebuf_, MAX_LINE_LEN, f_);
    if (l_)
    {
        do
        {
            prev_pos_ = cur_pos_;
            cur_pos_ = ftell (f_);
            if (*l_ == '>')
                break;
            add_seq ();
        }
        while ((l_ = fgets (linebuf_, MAX_LINE_LEN, f_)));
    }
    else
    {
        cur_pos_ = tot_len_;
    }
    seq_no_ ++;
    cur_reclen_ = prev_pos_ - cur_recstart_;
    seqbuf_ [seqlen_] = 0;
    return true;
}

bool FastaFile::seek (ulonglong off)
{
    if (!f_) return false;
    if (off >= tot_len_) return false;
    fseek (f_, off, SEEK_SET);
    prev_pos_ = off;
    l_ = fgets (linebuf_, MAX_LINE_LEN, f_);
    if (!l_ || *l_ != '>') ERR ("Invalid offset in fasta file: not at record start");
    cur_pos_ = ftell (f_);
    return true;
}

const char* FastaFile::cur_name ()  const
{
    return namebuf_;
}

const char* FastaFile::cur_hdr  () const
{
    return hdrbuf_;
}

const char* FastaFile::cur_seq  () const
{
    return seqbuf_;
}

char* FastaFile::cur_seq_buf  ()
{
    return seqbuf_;
}

unsigned FastaFile::cur_no  () const
{
    return seq_no_;
}

ulonglong FastaFile::tot_len () const
{
    return tot_len_;
}

ulonglong FastaFile::cur_pos () const
{
    return cur_pos_;
}

ulonglong FastaFile::cur_recstart () const
{
    return cur_recstart_;
}

unsigned FastaFile::cur_reclen () const
{
    return cur_reclen_;
}


#ifdef FASTA_TEST
int main (int argc, char* argv [])
{
    try
    {
        ulonglong opa [100];
        int c = 0;
        if (argc < 2) ERR ("No argument");
        FastaFile f (argv [1]);
        while (f.next () && c < 100)
        {
            ((int*) (opa + c)) [0] = f.cur_recstart ();
            ((int*) (opa + c)) [1] = f.cur_reclen ();
            c ++;
        }

        for (c = 0; c < 100; c ++)
        {
            f.seek (((int*) (opa + c)) [0]);
            f.next ();
            printf (">%s %s\n%s\n", f.cur_name (), f.cur_hdr (), f.cur_seq ());
        }

    }
    catch (Rerror& e)
    {
        printf ("%s\n", (const char*) e);
        return 1;
    }
    return 0;
}
#endif
