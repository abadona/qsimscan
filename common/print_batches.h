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

#ifndef PRINT_BATCHES_H
#define PRINT_BATCHES_H

#include <ostream>
#include "biosequence.h"
#include "weights.h"


void print_batches   (SEQ& xseq, SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 78);
void print_batches_ascii_subj (SEQ& xseq, const char* asciisubj, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 78);
void print_batches_3 (AA_SEQ& xseq, AA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 78);
void print_batches   (AA_SEQ& xseq, NA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 76);
void print_batches_an (AA_SEQ& xseq, NA_SEQ& yseq, BATCH *b_ptr, int b_cnt, WEIGHTS<int, 24>* w, std::ostream& stream, int margin = 5, int width = 76);
int print_batches_cigar (const BATCH *b_ptr, int b_cnt, char* dest, unsigned destlen);
#endif
