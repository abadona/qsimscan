
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

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "fasta.h"
#include "rerror.h"
#include "fileutils.h"
#include "portability.h"

static const char file_not_open [] = "INTERNAL ERROR: Operating on unopen fasta file";
static const int INIT_SEQ_LEN = 1000000; // 1M

void FastaFile::reset ()
{
    cur_pos_ = 0;
    prev_pos_ = 0;
    *linebuf_ = 0;
    seqbuf_ = NULL;
    seq_buf_sz_ = 0;
    *namebuf_ = 0;
    *hdrbuf_  = 0;
    l_ = linebuf_;
    nameptr_ = namebuf_;
    hdrptr_ = hdrbuf_;
    namelen_ = 0;
    hdrlen_ = 0;
    for (SeqCache::iterator ii = seq_cache_.begin (), sent = seq_cache_.end (); ii != sent; ii ++)
    {
        Seq& s = (*ii).second;
        delete [] s.name;
        delete [] s.hdr;
        delete [] s.seq;
    }
    seq_cache_.clear ();
}


void FastaFile::parse_hdr ()
{
    if (!*l_) return;
    assert (*l_ == '>');

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
    namelen_ = p1 - p;

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
        hdrlen_ = p - hdrbuf_;
    }
    else
    {
        *hdrbuf_ = 0;
        hdrlen_ = 0;
    }
    hdrptr_ = hdrbuf_;
    nameptr_ = namebuf_;
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
                // reallocate seqbuf_: increment space twice
                unsigned new_sz = seq_buf_sz_ * 2;
                char* newbuf = new char [new_sz + 1];
                if (!newbuf)
                {
                    ers << "ERROR: Not enough memory to hold the sequence" << namebuf_ << Throw;
                }
                memcpy (newbuf, seqbuf_, seq_buf_sz_);
                delete [] seqbuf_;
                seqbuf_ = newbuf;
                seq_buf_sz_ = new_sz;
            }
            seqbuf_ [seqlen_] = *p;
            seqlen_ ++;
        }
        p ++;
    }
}

FastaFile::FastaFile ()
:
f_ (NULL)
{
    reset ();
}

FastaFile::FastaFile (const char* name)
:
f_ (NULL)
{
    reset ();
    if (!open (name))
    {
        ers << "ERROR: Could not open fasta file" << name << Throw;
    }
}

FastaFile::~FastaFile ()
{
    close ();
}

bool FastaFile::open (const char* name)
{
    if (f_) 
        close ();
    else
        reset ();

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
        seq_no_ = -1;
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
        reset ();
        return true;
    }
    else
        return false;
}

bool FastaFile::next ()
{
    if (cur_pos_ == tot_len_) return false;
    cur_recstart_ = prev_pos_;
    
    // lookup if this record is already in cache
    SeqCache::iterator cur = seq_cache_.find (cur_recstart_);
    if (cur != seq_cache_.end ())
    {
        cur_reclen_ = (*cur).second.reclen;
        seqbuf_ = (*cur).second.seq;
        seqlen_ = (*cur).second.seqlen;
        hdrptr_ = (*cur).second.hdr;
        nameptr_ = (*cur).second.name;
        cur_pos_ += prev_pos_ + cur_reclen_;
    }
    else
    {
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
            while (l_ = fgets (linebuf_, MAX_LINE_LEN, f_));
    
        }
        else
        {
            cur_pos_ = tot_len_;
        }
        seq_no_ ++;
        cur_reclen_ = prev_pos_ - cur_recstart_;
        seqbuf_ [seqlen_] = 0;
        
        // place into cache
        Seq& newentry = seq_cache_ [cur_recstart_];
        newentry.reclen = cur_reclen_;
        newentry.seq = new char [seqlen_ + 1];
        memcpy (newentry.seq, seqbuf_, seqlen_+1);
        newentry.hdr = new char [hdrlen_ + 1];
        memcpy (newentry.hdr, hdrbuf_, hdrlen_+1);
        newentry.name = new char [namelen_ + 1];
        memcpy (newentry.name, namebuf_, namelen_+1);
    }
    return true;
}

bool FastaFile::seek (ulonglong off)
{
    if (!f_) return false;
    if (off >= tot_len_) return false;
    SeqCache::iterator cur = seq_cache_.find (off);
    prev_pos_ = off;
    if (cur == seq_cache_.end ())
    {
        fseek (f_, off, SEEK_SET);
        l_ = fgets (linebuf_, MAX_LINE_LEN, f_);
        if (!l_ || *l_ != '>') ERR ("Invalid offset in fasta file: not at record start");
        cur_pos_ = ftell (f_);
    }
    return true;
}

const char* FastaFile::cur_name ()  const
{
    return nameptr_;
}

const char* FastaFile::cur_hdr  () const
{
    return hdrptr_;
}

const char* FastaFile::cur_seq  () const
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
