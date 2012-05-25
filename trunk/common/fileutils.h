
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

#ifndef __fileutils_h__
#define __fileutils_h__

#pragma warning (disable:4786)
#include "structures.h"
#include "platform.h"

void expand_wildcards (const char* expr, std::vector <std::string>& files, bool verbose = true);
void checkCreateOutputFile (const char* output_name, bool overwrite, bool append);

StrVec      listdir (const char *path);    // returns directory listing
StrVec      listdir (const std::string& path);
StrVec      split_path (const char* path);  // splits path into components
StrVec      split_path (const std::string& path);
std::string join_path (StrVec components);  // merges components into path
std::string join_path (const char* comp0, ...); // last component must be NULL
bool        file_exists (const char* name); // checks for file existence
bool        file_exists (const std::string& name);
bool        is_file (const char* fname);
bool        is_file (const std::string& fname);
bool        is_dir (const char* fname);
bool        is_dir (const std::string& fname);
time_t      file_time (const char* fname);
time_t      file_time (const std::string& fname);
ulonglong   file_size (const char* fname);
ulonglong   file_size (const std::string& fname);
ulonglong   get_open_file_size (int fhandle);

// warning: this function is HEAVY (reads directory, uses heap-allocated data etc.; also returns std::string)
// warning: this function creates race condition if used simultanously with same prefix on same directory
std::string temp_dir (const char* tmpdir = NULL);
std::string make_temp_fname (const char* tmpdir = NULL, const char* prefix = NULL);
int make_linked_temp_file (std::string& dest, const char* tmpdir = NULL, const char* prefix = NULL);
int make_temp_file ();


class LineReader
{
private:
    int fhandle_;
    char *buffer_;
    int cur_line_beg_;
    int cur_line_end_;
    int buf_end_;
    char prev_char_;
    longlong cur_pos_;

    void nextChunk ();
public:
    LineReader (const char* fname);
    ~LineReader ();
    char* nextLine ();
    bool isOpen () const;
    longlong curPos () const;
    void close ();
};

#undef putc
#define WRITE_BUF_SZ 1024*1024*4 // 4Mb
class BufferedWriter
{
    int fhandle;
    char *buffer;
    int bufst;
    int cpos;
public:
    BufferedWriter (const char* fname = NULL, bool append = false);
    ~BufferedWriter ();
    bool open (const char* fname, bool append = false);
    void close ();
    void flush ();
    bool is_open ()
    {
        return (fhandle != -1 && buffer);
    }
    void write (char* buf, unsigned blen)
    {
        while (blen--) putc (*buf++);
    }
    int puts (const char* s)
    {
        int toR = 0;
        while (*s) putc (*s++), toR ++;
        return toR;
    }
    void putc (char c)
    {
        if (cpos + 1 >= WRITE_BUF_SZ) flush ();
        buffer [cpos++] = c;
    }
    int tell ()
    {
        return bufst + cpos;
    }
};


StrVec read_lists (const StrVec& fnames);
StrVec read_list (const char* fname);
StrVec read_list (std::string fname);

#endif

