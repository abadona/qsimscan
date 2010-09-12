
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

#ifndef PRINT_BATCHES_H
#define PRINT_BATCHES_H

#include <ostream>
#include "biosequence.h"
#include "weights.h"


void print_batches   (SEQ& xseq, SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 78);
void print_batches_3 (AA_SEQ& xseq, AA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 78);
void print_batches   (AA_SEQ& xseq, NA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 76);
void print_batches_an (AA_SEQ& xseq, NA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 76);
int print_batches_cigar (const BATCH *b_ptr, int b_cnt, char* dest, unsigned destlen);
#endif
