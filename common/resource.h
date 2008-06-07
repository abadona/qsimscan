
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

#ifndef __RESOURCE_H__
#define __RESOURCE_H__

#include <stdio.h>
#include "rerror.h"

class FileWrapper
{
    FILE* f_;
    mutable bool controlled_;
public:

    FileWrapper (const char* fname, const char* mode, bool control = true)
    {
        f_ = fopen (fname, mode);
        controlled_ = control;
    }
    FileWrapper (FILE* f, bool control = true)
    :
    f_ (f),
    controlled_ (control)
    {
    }
    FileWrapper ()
    :
    f_ (NULL),
    controlled_ (true)
    {
    }
    FileWrapper (const FileWrapper& from)
    :
    f_ (from.f_),
    controlled_ (from.controlled_)
    {
        from.controlled_ = false;
    }
    ~FileWrapper ()
    {
        free ();
    }
    FileWrapper& operator = (const FileWrapper& from)
    {
        free ();
        f_ = from.f_;
        controlled_ = from.controlled_;
        from.controlled_ = false;
        return *this;
    }
    FileWrapper& operator = (FILE* nf)
    {
        free ();
        f_ = nf;
        controlled_ = true;
        return *this;
    }
    bool operator ! () const
    {
        return (f_ == NULL);
    }
    void free ()
    {
        if (controlled_ && f_)
        {
            fclose (f_);
            f_ = NULL;
        }
    }
    FILE* release ()
    {
        controlled_ = false;
        return f_;
    }
    FILE* operator* ()
    {
        return f_;
    }
};

template <class ValueType> class MemWrapper
{
    ValueType* ptr_;
    mutable bool controlled_;
public:
    MemWrapper (unsigned size, bool control = true)
    {
        if (size)
        {
            ptr_ = new ValueType [size];
            if (!ptr_) Error (MemoryRerror);
        }
        else
            ptr_ = NULL;
        controlled_ = control;
    }
    MemWrapper (ValueType* ptr, bool control = true)
    :
    ptr_ (ptr),
    controlled_ (control)
    {
    }
    MemWrapper ()
    :
    ptr_ (NULL),
    controlled_ (false)
    {
    }
    MemWrapper (const MemWrapper& from)
    :
    ptr_ (from.ptr_),
    controlled_ (from.controlled_)
    {
        from.controlled_ = false;
    }
    ~MemWrapper ()
    {
        free ();
    }
    MemWrapper& operator = (const MemWrapper& from)
    {
        free ();
        ptr_ = from.ptr_;
        controlled_ = from.controlled_;
        from.controlled_ = false;
        return *this;
    }
    MemWrapper& operator = (ValueType* np)
    {
        free ();
        ptr_ = np;
        controlled_ = true;
        return *this;
    }
    void free ()
    {
        if (controlled_ && ptr_)
        {
            delete [] ptr_;
            ptr_ = NULL;
        }
    }
    bool operator! ()
    {
        return ptr_ == NULL;
    }
    ValueType* release ()
    {
        controlled_ = false;
        return ptr_;
    }
    ValueType& operator* ()
    {
        return *ptr_;
    }
    ValueType* operator-> ()
    {
        return ptr_;
    }
//    ValueType& operator [] (unsigned idx)
//    {
//        return ptr_ [idx];
//    }
//    ValueType& operator [] (int idx)
//    {
//        return ptr_ [(unsigned) idx];
//    }
//    ValueType operator [] (unsigned idx) const
//    {
//        return ptr_ [idx];
//    }
//    ValueType operator [] (int idx) const
//    {
//        return ptr_ [(unsigned) idx];
//    }
    operator ValueType* ()
    {
        return ptr_;
    }
    operator const ValueType* () const
    {
        return ptr_;
    }
};

template <class ValueType> class ConstMemWrapper
{
    const ValueType* ptr_;
    mutable bool controlled_;
public:
    ConstMemWrapper (unsigned size, bool control = true)
    {
        if (size)
        {
            ptr_ = new ValueType [size];
            if (!ptr_) Error (MemoryRerror);
        }
        else
            ptr_ = NULL;
        controlled_ = control;
    }
    ConstMemWrapper (const ValueType* ptr, bool control = true)
    :
    ptr_ (ptr),
    controlled_ (control)
    {
    }
    ConstMemWrapper ()
    :
    ptr_ (NULL),
    controlled_ (false)
    {
    }
    ConstMemWrapper (const ConstMemWrapper& from)
    :
    ptr_ (from.ptr_),
    controlled_ (from.controlled_)
    {
        from.controlled_ = false;
    }
    ~ConstMemWrapper ()
    {
        free ();
    }
    ConstMemWrapper& operator = (const ConstMemWrapper& from)
    {
        free ();
        ptr_ = from.ptr_;
        controlled_ = from.controlled_;
        from.controlled_ = false;
        return *this;
    }
    ConstMemWrapper& operator = (const ValueType* np)
    {
        free ();
        ptr_ = np;
        controlled_ = true;
        return *this;
    }
    void free ()
    {
        if (controlled_ && ptr_)
        {
            delete [] (ValueType*) ptr_;
            ptr_ = NULL;
        }
    }
    bool operator! ()
    {
        return ptr_ == NULL;
    }
    const ValueType* release ()
    {
        controlled_ = false;
        return ptr_;
    }
    const ValueType& operator* ()
    {
        return *ptr_;
    }
    const ValueType* operator-> ()
    {
        return ptr_;
    }
    const ValueType& operator [] (unsigned idx)
    {
        return ptr_ [idx];
    }
    const ValueType& operator [] (int idx)
    {
        return ptr_ [(unsigned) idx];
    }
    ValueType operator [] (unsigned idx) const
    {
        return ptr_ [idx];
    }
    ValueType operator [] (int idx) const
    {
        return ptr_ [(unsigned) idx];
    }
    operator const ValueType* ()
    {
        return ptr_;
    }
    operator const ValueType* () const
    {
        return ptr_;
    }
};

template <class ValueType> class ObjWrapper
{
    ValueType* ptr_;
    mutable bool controlled_;
public:
    ObjWrapper (ValueType* ptr, bool control = true)
    :
    ptr_ (ptr),
    controlled_ (control)
    {
    }
    ObjWrapper ()
    :
    ptr_ (NULL),
    controlled_ (false)
    {
    }
    ObjWrapper (const ObjWrapper& from)
    :
    ptr_ (from.ptr_),
    controlled_ (from.controlled_)
    {
        from.controlled_ = false;
    }
    ~ObjWrapper ()
    {
        free ();
    }
    ObjWrapper &operator = (const ObjWrapper& from)
    {
        free ();
        ptr_ = from.ptr_;
        controlled_ = from.controlled_;
        from.controlled_ = false;
        return *this;
    }
    ObjWrapper &operator = (ValueType* np)
    {
        free ();
        ptr_ = np;
        controlled_ = true;
        return *this;
    }
    void free ()
    {
        if (controlled_ && ptr_)
        {
            delete ptr_;
            ptr_ = NULL;
        }
    }
    bool operator! ()
    {
        return ptr_ == NULL;
    }
    ValueType* release ()
    {
        controlled_ = false;
        return ptr_;
    }
    ValueType& operator* ()
    {
        return *ptr_;
    }
    ValueType* operator-> ()
    {
        return ptr_;
    }
    operator ValueType* ()
    {
        return ptr_;
    }
};


#endif
