
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

#ifndef __filters_h__
#define __filters_h__
#include "sequtil.h"
#include "weights.h"

/*
calculates nucleotide homology level in %
returned level is corrected by nucleotide composition probability
*/
float nu_score_b (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, int* matc);

/*
calculates aminoacid homology level in % in 6 frames
normalizes similarity to Y sequence self-similarity level (100%)
returns highest scoring Y - sequence frame number
scorexy contains all 6 frame scores

this function works with single batch only
*/
int aa_score (char* xseq, int xpos, char* yseq, int ypos, int len, int* na_wt[], int* scorexy);

/*
calculates protein alignment score on the protein similarity batches
also calculates auto-similarity scores on x and y fragment sets comprising similarity
*/
int prot_score (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, WEIGHTS<int, 24>* wm, int* x_auto = NULL, int* y_auto = NULL);

/*
calculates protein alignment score on the protein-to-nucleotide similarity batches
also calculates auto-similarity scores on x and y fragment sets comprising similarity
*/
int an_score (char* xaaseq, char* ynnseq, BATCH* b_ptr, int b_cnt, bool aafirst, WEIGHTS<int, 24>* wm, int* x_auto = NULL, int* y_auto = NULL);

/*
calculates chi^2 statistics for nucleotide similarity batch based on following assumptions:
	1. event space - 5 (a-match, g-match, c-match, t-match, mismatch), degrees of freedom - 4
	2. local (calculated on a batch) nucletide frequencies are used as reference
	(local composition bias is not considered meaningful by this hypothesis)
chi^2 values can be directly compared as long as both chi^2 values > chi^2 distribution mode
(pass1/pass2 selection heuristics almost guarantees that condition)
*/
float nu_matches_chi2_b (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt);

/*
calculates chi^2 statistics for nucleotide similarity batch based on following assumptions:
	1. event space - 16 (aa, ag, ac, at, ga, gg, ... tt), degrees of freedom - 15 (or 9)
	2. local (calculated on a batch) nucletide frequencies are used as reference
	(local composition bias is not considered meaningful by this hypothesis)
chi^2 values can be directly compared as long as both chi^2 values > chi^2 distribution mode
(pass1/pass2 selection heuristics almost guarantees that condition)
*/
float nu_all_chi2_b (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt);

/*
calculates an array of similarity scores (matches) with variable offset
between X and Y sequences (neighbor diagonals), 0 < offs < max_offs
returns offs != 0 with maximum score
*/
int nu_offs_score_b (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, int max_offs, int* offs_scores);

/*
calculates modulo-3 positional mismatch statistics for each y sequence
(does not initialise scores)
*/
int nu_mod3_score_accum_b  (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, int* y_mod3_scores);

/*
calculates full length positional mismatch statistics for each y sequence
(does not initialise scores)
*/
int nu_pos_score_accum_b  (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, int* mism, int* matc, bool revers = false, int ylen = 0);


//calculates nucleotide homology level in %
//returned level is corrected by nucleotide composition probability
int nu_score (char* xseq, int xpos, char* yseq, int ypos, int len, int* matc);

/*
calculates chi^2 statistics for nucleotide similarity batch based on following assumptions:
	1. event space - 5 (a-match, g-match, c-match, t-match, mismatch), degrees of freedom - 4
	2. local (calculated on a batch) nucletide frequencies are used as reference
	(local composition bias is not considered meaningful by this hypothesis)
chi^2 values can be directly compared as long as both chi^2 values > chi^2 distribution mode
(pass1/pass2 selection heuristics almost guarantees that condition)
*/
float nu_matches_chi2 (char* xseq, int xpos, char* yseq, int ypos, int len);

/*
calculates chi^2 statistics for nucleotide similarity batch based on following assumptions:
	1. event space - 16 (aa, ag, ac, at, ga, gg, ... tt), degrees of freedom - 15 (or 9)
	2. local (calculated on a batch) nucletide frequencies are used as reference
	(local composition bias is not considered meaningful by this hypothesis)
chi^2 values can be directly compared as long as both chi^2 values > chi^2 distribution mode
(pass1/pass2 selection heuristics almost guarantees that condition)
*/
float nu_all_chi2 (char* xseq, int xpos, char* yseq, int ypos, int len);

/*
calculates an array of similarity scores (matches) with variable offset
between X and Y sequences (neighbor diagonals), 0 < offs < max_offs
returns offs != 0 with maximum score
*/
int nu_offs_score (char* xseq, int xpos, char* yseq, int ypos, int len, int max_offs, int* offs_scores);


bool rep_filt (char* xseq, int xpos, char* yseq, int ypos, int len);


#endif
