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

#ifndef __gaplibshim_h__
#define __gaplibshim_h__

#include "biosequence.h"
#include "weights.h"

typedef WEIGHTS <int, 24> WMatrixType;
#define WMatrixXY(m, x, y) m.mx[x][y]
typedef double ResFreqType[AANUM];


#endif /* __gaplibshim_h__ */
