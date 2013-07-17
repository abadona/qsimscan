
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

#define __ncbi_weight_matrix_cpp__

#include "resource.h"
#include <string.h>
#include "ncbi_weight_matrix.h"

const char* ERR_BadMatrix = "Bad weight matrix format";
#ifndef _MSC_VER
const int STD_PROT_MATRIX_SIZE = 24;
#endif

const unsigned BUFSZ = 2048;

unsigned readNcbiMatrix (const char* fname, unsigned max_alpha_size, char* alphaBuf, int* valueBuf)
{
    char buf [BUFSZ];
    int tmp;
    unsigned col = 0, row = 0;
    unsigned alphaSize = 0;

    FileWrapper fi (fname, "rt");
    if (!fi) 
        ers << fname << ". This or other aminoacid substitution weight matrices can be downloaded from NCBI (ftp://ftp.ncbi.nlm.nih.gov/blast/matrices)" << ThrowEx(FileNotFoundRerror);

    while (fgets (buf, BUFSZ, *fi))
    {
        if (buf [0] == '#') continue;
        if (buf [0] == ' ')
        {
            // read alphabet
            char* tok = strtok (buf, " \t\n");
            while (tok)
            {
                if (alphaSize == max_alpha_size) ers << "Alphabet size overflow in " << fname << ThrowEx(BadMatrixFormat);
                if (tok [1] != 0) ers << "Multichar symbol found in alphabet in " << fname << ThrowEx(BadMatrixFormat);
                alphaBuf [alphaSize] = tok [0];
                alphaSize ++;
                tok = strtok (NULL , " \t\n");
            }
        }
        else
        {
            if (row >= alphaSize) ers << "Too many rows in " << fname << ThrowEx(BadMatrixFormat);
            strtok (buf, " \t\n");
            for (col = 0; col < alphaSize; col++)
            {
                char* tok = strtok (NULL, " \t\n");
                if (!tok) ers << "row is too short in " << fname << ThrowEx(BadMatrixFormat);
                if (sscanf (tok, "%d", &tmp) != 1) ers << "non-integer value in " << fname << ThrowEx(BadMatrixFormat);
                valueBuf [row * alphaSize + col] = tmp;
            }
            row ++;
        }
    }
    if (row < alphaSize) ers << "too few rows in " << fname << ThrowEx(BadMatrixFormat);

    return alphaSize;
}


