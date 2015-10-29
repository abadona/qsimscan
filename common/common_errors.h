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

#ifndef __common_errors_h__
#define __common_errors_h__

#ifndef __common_errors_cpp__

extern const char* ERR_NoMemory;
extern const char* ERR_Internal;
extern const char* ERR_FileNotFound;
extern const char* ERR_OSError;
extern const char* ERR_OutOfBounds;


// synonyms
#define NOEMEM ERR_NoMemory
#define INTERNAL ERR_Internal

#define MAKE_ERROR_TYPE(C,N) class C : public Rerror\
{\
public:\
    C (const char* s = "")\
    : Rerror (N)\
    {\
        if (*s)\
        {\
            msg_ += ": ";\
            msg_ += s;\
        }\
    }\
};

MAKE_ERROR_TYPE (MemoryRerror, ERR_NoMemory);
MAKE_ERROR_TYPE (InternalRerror, ERR_Internal);
MAKE_ERROR_TYPE (FileNotFoundRerror, ERR_FileNotFound);
MAKE_ERROR_TYPE (OutOfBoundsRerror, ERR_OutOfBounds);

#endif

class OSRerror : public Rerror
{
    const char* get_err_str () const;
    const char* get_errno_str () const;
public:
    OSRerror (const char* s = "")
    : Rerror (ERR_OSError)
    {
        msg_ += ": errno ";
        msg_ += get_errno_str ();
        msg_ += ": ";
        msg_ += get_err_str ();
        if (*s)
        {
            msg_ += ": ";
            msg_ += s;
        }
    }
};


#endif // __common_errors_h__
