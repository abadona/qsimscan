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

#include "misc_utils.h"
#include "rerror.h"
#include "fileutils.h"
#ifndef _MSC_VER
#include <unistd.h>
#endif


void checkCreateFile (bool overwrite, const char* output_name, std::ofstream& out_file)
{
    if (file_exists (output_name))
    {
        if (!overwrite) // ask weather to overwrite existing user object
            ers << "Output file with the name " << output_name << " allready exists" << Throw;
        if (unlink (output_name)) ers << "Unable to overwrite output file" << Throw;
    }
    out_file.open (output_name);
    if (!out_file.is_open ()) ers << "Unable to open output file" << Throw;
}
