
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

#ifndef __align_h__
#define __align_h__
#include <sciminmax.h>
#include <stdlib.h>
#include <iostream>
#include <rerror.h>
#include "sequtil.h"
#include "biosequence.h"
#include "weights.h"

#define ALIGN_DIAG 0
#define ALIGN_DOWN 1
#define ALIGN_LEFT 2
#define ALIGN_STOP 3

#define ALIGN_HSKIP 4
#define ALIGN_VSKIP 8

// #pragma pack (push, 4)

struct ALIGN_VECT
{
    int w;
    int h;
    int r;
} __attribute__ ((packed));

// #pragma pack (pop)

class ALIGN
{
    int xref, yref, xstep, bstep;
    int max_ylen, max_size;
    int max_w, max_x, max_y;
    char* max_bp;
    char* btrmx;
    ALIGN_VECT* ap_1;
    ALIGN_VECT* ap_3[3];
    WEIGHTS<int, 24>* w;

    /*
    member inner loops used for finding alignment
    fill backtrace matrix and store max score position
    */
    void align_y_loop (register ALIGN_VECT* ap, int *wp, int gip, int gep, char* bp, int x, int y, int len);
    void align_y_loop_n2a (register ALIGN_VECT* ap1, register ALIGN_VECT* ap3, int *wp, int gip, int gep, char* bp, int x, int y, int len);

public:

    ALIGN (WEIGHTS<int, 24>* w, int max_ylen, int max_size);
    ALIGN (WEIGHTS<int, 24>* w, int max_ylen);
    ~ALIGN ();

    /*
    calculates best local alignment between sequence pair of the same type
    returns maximum local alignment score
    */
    int  align (SEQ& xseq, SEQ& yseq);

    /*
    calculates best local alignment between sequence pair of different type using aa weighting
    returns maximum local alignment score
    */
    int  align_na (AA_SEQ& yseq, NA_SEQ& xseq);

    /*
    calculates best local alignment between sequence pair of the same type
    on a diagonal band (len, diag +- width)
    returns maximum local alignment score
    NOTE: batch xpos, ypos, len should be inside X and Y sequences, width > 0
    */
    int  align_band (SEQ& xseq, SEQ& yseq, int xpos, int ypos, int len, int width);

    /*
    follows backtrace matrix, fills BATCH array, returns number of batches
    */
    int  backtrace (BATCH *b_ptr, int max_cnt, unsigned width = 0);

    int get_max_x () const { return max_x; }
    int get_max_y () const { return max_y; }
};


/*
non-member inner loops used for search only
do not fill backtrace matrix, only finds score
*/
int align_y_loop_so (register ALIGN_VECT* ap, int *wp, int gip, int gep, int max_w, int len);
int align_y_loop_n2a_so (register ALIGN_VECT* ap1, register ALIGN_VECT* ap3, int *wp, int gip, int gep, int max_w, int len);


//in-place array element order reversal
template <class T>
void reverse (T* pb, int cnt)
{
    T tmp;
    T* pe = pb + cnt - 1;

    cnt >>= 1;
    while (cnt-- > 0)
        tmp = *pb, *pb++ = *pe, *pe-- = tmp;
}

#endif // __align_h__
