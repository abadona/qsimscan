
//////////////////////////////////////////////////////////////////////////////
// This software module is developed by SCIDM team in 1999-2015.
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
