
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

#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <fstream>
#include <ios>

#include "portability.h"

// file access mode flags
#if defined (_MSC_VER)
    #include <io.h>
    #define WRFLAGS O_CREAT|O_WRONLY|O_BINARY
#else
    #include <unistd.h>
    #define WRFLAGS O_CREAT|O_WRONLY|O_LARGEFILE
#endif

#include "portability.h"
#include "fileutils.h"
#include "rerror.h"
#include "common_errors.h"

// directory listing
#if defined (_MSC_VER)
    const char PATH_SEPARATOR = '\\';
    #include <windows.h>
#else
    #include <dirent.h>
    #define NAMLEN(dirent) strlen((dirent)->d_name)
    const char PATH_SEPARATOR = '/';
#endif


// globbing (wildcard expansion)
#if defined (_MSC_VER)

void expand_wildcards (const char* expr, std::vector <std::string>& files, bool verbose)
{
    struct _finddata_t gbfile;
    long hFile;
    if ((hFile = _findfirst (expr, &gbfile)) != -1L)
    {
        char drive [_MAX_DRIVE]; char dir [_MAX_DIR]; char name [_MAX_FNAME]; char ext [_MAX_EXT]; char fullp [_MAX_PATH];
        _splitpath (expr, drive, dir, name, ext);
        do
        {
            _makepath (fullp, drive, dir, gbfile.name, "");
            files.push_back (fullp);
        }
        while (_findnext (hFile, &gbfile) == 0);
        _findclose (hFile);
    }
}
#else

#include <glob.h>

static std::string glob_errpath;
static int glob_errno = 0;

int glob_errf (const char *epath, int eerrno)
{
    glob_errpath = epath;
    glob_errno = eerrno;
}

//int glob_errf (const char *epath, int eerrno)
//{
//    ers << "Glob error on: " << epath << "; error " << eerrno << ": " << strerror (eerrno) << Throw;
//}

void expand_wildcards (const char* expr, std::vector <std::string>& files, bool verbose)
{
    glob_t gdata;
    gdata.gl_offs = 0;
    int glob_res = glob (expr, 0, glob_errf, &gdata);
    if (glob_errno)
        ers << "Glob error on: " << glob_errpath.c_str () << "; error " << glob_errno << ": " << strerror (glob_errno) << Throw;
    switch (glob_res)
    {
        case 0:
                            int mm;
                            for (mm = 0; mm < gdata.gl_pathc; mm ++)
                                files.push_back ((char*) gdata.gl_pathv [mm]);
        case GLOB_NOMATCH:
        case GLOB_NOSPACE:
        case GLOB_ABORTED:  break;
    }
    globfree (&gdata);
}

#endif

void checkCreateOutputFile (const char* output_name, bool overwrite, bool append)
{
    // walk to the parent dir, creating directories as needed
    std::ifstream tifile;
    tifile.open (output_name);
    if (tifile.is_open ())
    {
        tifile.close ();
        if (!overwrite && !append)
            ers << "Output file allready exist : " << output_name << " and neither overwite nor append are enabled" << Throw;
        else if (!append)
        {
            if (unlink (output_name) != 0)
                ers << "Unable to delete existing object: " << output_name << ThrowEx (OSRerror);
            std::ofstream tofile (output_name, std::ios::out);
            if (!tofile)
                ers << "Error creating new output file: " << output_name << ThrowEx (OSRerror);
            tofile.close ();
            if (unlink (output_name) != 0)
                ers << "Output file " << output_name << " is created but cannot be unlinked" << ThrowEx (OSRerror);
        }
        else
        {
            std::ofstream tofile (output_name, std::ios::app);
            if (!tofile)
                ers << "Error opening output file for append: " << output_name << ThrowEx (OSRerror);
            tofile.close ();
        }
    }
    else
    {
            std::ofstream tofile (output_name, std::ios::out);
            if (!tofile)
                ers << "Error creating new output file: " << output_name << ThrowEx (OSRerror);
            tofile.close ();
            if (unlink (output_name) != 0)
                ers << "Output file " << output_name << " is created but cannot be unlinked" << ThrowEx (OSRerror);
    }
}



StrVec listdir (const char *name)
{
    StrVec file_list;
#ifdef _MSC_VER
    HANDLE hFindFile;
    WIN32_FIND_DATA FileData;
    char namebuf [MAX_PATH*2+5];

    strncpy (namebuf, name, sizeof (namebuf));

    int len = strlen (namebuf);
    if (len > 0)
    {
        char ch = namebuf [len-1];
        if (ch != PATH_SEPARATOR && ch != ':')
            namebuf [len++] = PATH_SEPARATOR;
    }
    strncpy (namebuf + len, "*", sizeof (namebuf) - len);
    hFindFile = FindFirstFile (namebuf, &FileData);
    if (hFindFile == INVALID_HANDLE_VALUE) return file_list;
    do {
        if (FileData.cFileName[0] == '.' &&
            (FileData.cFileName[1] == '\0' ||
            FileData.cFileName[1] == '.' &&
            FileData.cFileName[2] == '\0'))
            continue;
        file_list.push_back (FileData.cFileName);
    } while (FindNextFile (hFindFile, &FileData) == TRUE);
    FindClose (hFindFile);
#else
    DIR *dirp;
    struct dirent *ep;
    if ((dirp = opendir (name)) == NULL) {
        return file_list;
    }
    while ((ep = readdir (dirp)) != NULL) {
        if (ep->d_name[0] == '.' &&
            (NAMLEN (ep) == 1 ||
            (ep->d_name[1] == '.' && NAMLEN(ep) == 2)))
            continue;
        file_list.push_back (std::string(ep->d_name, NAMLEN (ep)));
    }
    closedir(dirp);
#endif
    return file_list;
}

StrVec listdir (const std::string& path)
{
    return listdir (path.c_str ());
}

StrVec split_path (const char* path)
{
    StrVec rv;
    std::string pstr = path;
    std::string::size_type pos, ppos = 0;
    while ((pos = pstr.find (PATH_SEPARATOR, ppos)) != std::string::npos)
    {
        if (pos != ppos || ppos == 0)
            rv.push_back (pstr.substr (ppos, pos - ppos));
        ppos = pos + 1;
    }
    if (ppos < pstr.length ())
        rv.push_back (pstr.substr (ppos, pstr.length () - ppos));

    return rv;
}

StrVec split_path (const std::string& path)
{
    return split_path (path.c_str ());
}

std::string join_path (StrVec components)
{
    printf ("\njoin_path\n");
    std::string rv = "";
    StrVec::const_iterator itr = components.begin ();
    for (;itr != components.end (); itr ++)
    {
        rv += *itr;
        int rvlen = rv.length ();
        if ((itr == components.begin () || components.end () - itr > 1) && (rvlen || itr == components.begin ()) && (rv [rvlen - 1] != PATH_SEPARATOR))
            rv += PATH_SEPARATOR;
    }
    return rv;
}

std::string join_path (const char* comp0, ...)
{
    std::string rv = "";
    va_list argl;
    va_start (argl, comp0);
    const char* comp = comp0;
    int cnt = 0;
    while (comp)
    {
        if (comp != comp0)
        {
            int rvlen = rv.length ();
            if (rvlen && rv [rvlen-1] != PATH_SEPARATOR)
                rv += PATH_SEPARATOR;
        }
        rv += comp;
        comp = va_arg (argl, const char*);
        cnt ++;
    }
    va_end (argl);
    if (cnt == 1 && comp0 && !*comp0) rv += PATH_SEPARATOR;
    return rv;
}

#if !defined (S_ISDIR)
#define S_ISDIR(x) ((x)&S_IFDIR)
#endif
#if !defined (S_ISREG)
#define S_ISREG(x) ((x)&S_IFREG)
#endif

bool file_exists (const char* fname)
{
    return is_file (fname);
}
bool file_exists (const std::string& name)
{
    return file_exists (name.c_str ());
}

bool is_file (const char* fname)
{
    struct sci_stat_struc stat_buf;
    if (-1 == sci_stat (fname, &stat_buf)) return false;
    return S_ISREG (stat_buf.st_mode) != 0;
}
bool is_file (const std::string& name)
{
    return is_file (name.c_str ());
}

bool is_dir (const char* fname)
{
    struct sci_stat_struc stat_buf;
    if (-1 == sci_stat (fname, &stat_buf)) return false;
    return S_ISDIR (stat_buf.st_mode) != 0;
}
bool is_dir (const std::string& name)
{
    return is_dir (name.c_str ());
}

time_t file_time (const char* fname)
{
    struct sci_stat_struc stat_buf;
    if (-1 == sci_stat (fname, &stat_buf)) return 0;
    return stat_buf.st_mtime;
}
time_t file_time (const std::string& name)
{
    return file_time (name.c_str ());
}

ulonglong file_size (const char* fname)
{
    struct sci_stat_struc stat_buf;
    if (-1 == sci_stat (fname, &stat_buf)) return 0;
    return stat_buf.st_size;
}
ulonglong file_size (const std::string& name)
{
    return file_size (name.c_str ());
}

#define LR_BUF_SIZE 1024*256

#ifndef O_BINARY
#define O_BINARY 0
#endif

#ifndef O_SEQUENTIAL
#define O_SEQUENTIAL 0
#endif

LineReader::LineReader (const char* fname)
:
fhandle_ (-1),
buffer_ (NULL),
cur_pos_ (-1L)
{
    fhandle_ = sci_open (fname, O_BINARY|O_RDONLY|O_SEQUENTIAL);
    if (fhandle_ == -1) ers << "Unable to open file " << fname << " for reading." << std::endl << strerror (errno) << Throw;
    try
    {
        buffer_ = new char [LR_BUF_SIZE+1]; // +1 is needed because there may be a need to zero-terminate the full buffer
    }
    catch (std::bad_alloc&)
    {
    }
    if (!buffer_) Error (MemoryRerror);
    cur_line_beg_ = 0;
    cur_line_end_ = 0;
    buf_end_ = 0;
    prev_char_ = 0;
}
LineReader::~LineReader ()
{
    close ();
}

void LineReader::nextChunk ()
{
    int reminder = buf_end_ - cur_line_beg_;
    if (reminder)
    {
        if (!cur_line_beg_)
            ers << "line in file being read too ling to fit in buffer" << Throw;
        memmove (buffer_, buffer_ + cur_line_beg_, reminder);
    }
    cur_line_end_ -= cur_line_beg_;
    cur_line_beg_ = 0;
    buf_end_ = reminder;
    int rd = sci_read (fhandle_, buffer_ + buf_end_, LR_BUF_SIZE - buf_end_);
    buf_end_ += rd;
}

char* LineReader::nextLine ()
{
    char* rv;
    char cur_char;
    while (1)
    {
        cur_pos_ ++;
        if (cur_line_end_ == buf_end_)
        {
            nextChunk ();
            if (cur_line_end_ == buf_end_) return NULL;
        }
        cur_char = buffer_ [cur_line_end_];
        if (cur_char == '\r' || cur_char == '\n')
        {
            rv = buffer_ + cur_line_beg_;
            buffer_ [cur_line_end_] = 0;
            cur_line_beg_ = cur_line_end_ + 1;
            cur_line_end_ = cur_line_beg_;

            if (prev_char_ == '\r' && cur_char == '\n')
            {
                prev_char_ = cur_char;
                continue;
            }
            else
            {
                prev_char_ = cur_char;
                return rv;
            }
        }
        cur_line_end_ ++;
    }
}

longlong LineReader::curPos () const
{
    return cur_pos_;
}

bool LineReader::isOpen () const
{
    return fhandle_ != -1;
}
void LineReader::close ()
{
    if (fhandle_ != -1)
    {
        ::sci_close (fhandle_);
        fhandle_ = -1;
    }
    if (buffer_)
    {
        delete [] buffer_;
        buffer_ = NULL;
    }
}

StrVec read_lists (const StrVec& fnames)
{
    StrVec rv;
    StrVec::const_iterator itr = fnames.begin ();
    while (itr != fnames.end ())
    {
        StrVec vals = read_list (*itr);
        std::copy (vals.begin (), vals.end (), std::back_inserter (rv));
        itr ++;
    }
    return rv;
}
StrVec read_list (const char* fname)
{
    StrVec rv;
    LineReader rd (fname);
    char *l, *lc;
    while (l = rd.nextLine ())
    {
        // find first token and add to rv
        while (*l && isspace (*l)) l ++;
        if (!*l) continue;
        lc = l + 1;
        while (*lc && !isspace (*lc)) lc ++;
        *lc = 0;
        rv.push_back (l);
    }
    rd.close ();
    return rv;
}
StrVec read_list (const std::string fname)
{
    return read_list (fname.c_str ());
}


BufferedWriter::BufferedWriter (const char* fname, bool append)
:
fhandle (-1),
buffer (NULL),
cpos (0),
bufst (0)
{
    if (fname) sci_open (fname, append);
}
BufferedWriter::~BufferedWriter ()
{
    close ();
}
bool BufferedWriter::open (const char* fname, bool append)
{
    if (is_open ()) ers << "writer allready open" << Throw;
    int flags = WRFLAGS;
    if (append) flags |= O_APPEND;
    else flags |= O_TRUNC;
    fhandle = ::sci_open (fname, flags, S_IREAD | S_IWRITE);
    if (fhandle == -1)
        return false;
    if (!buffer)
    {
        try
        {
            buffer = new char [WRITE_BUF_SZ];
        }
        catch (std::bad_alloc&)
        {
        }
        if (!buffer) Error (MemoryRerror);
    }
    return true;
}
void BufferedWriter::close ()
{
    flush ();
    if (fhandle != -1)
    {
        ::sci_close (fhandle);
        fhandle = -1;
    }
    if (buffer)
    {
        delete [] buffer;
        buffer = NULL;
    }
}
void BufferedWriter::flush ()
{
    if (fhandle != -1 && buffer && cpos)
    {
        ::sci_write (fhandle, buffer, cpos);
        bufst += cpos;
        cpos = 0;
    }
}

