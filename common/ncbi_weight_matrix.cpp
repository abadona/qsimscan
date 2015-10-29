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


