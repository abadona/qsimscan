
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
