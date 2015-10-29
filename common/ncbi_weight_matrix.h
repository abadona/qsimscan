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

#ifndef __ncbi_weight_matrix_h__
#define __ncbi_weight_matrix_h__

#include "rerror.h"

//#ifndef __ncbi_weight_matrix_cpp__
extern const char* ERR_BadMatrix;

#ifndef _MSC_VER
extern const int STD_PROT_MATRIX_SIZE;
#else
const int STD_PROT_MATRIX_SIZE = 24;
#endif

MAKE_ERROR_TYPE (BadMatrixFormat, ERR_BadMatrix);


// reads standard (square) NCBI weight matrix file;
// saves read alphabet and values into passed buffers.
// assumes alpha_buf is at least max_alpha_size bytes long and value_buf is at least max_alpha_size * max_alpha_size bytes long
// returns actual alphabet size
// on format error, throws BadMatrixFormat
unsigned readNcbiMatrix (const char* fname, unsigned max_alpha_size, char* alpha_buf, int* value_buf);

#endif // __ncbi_weight_matrix_h__
