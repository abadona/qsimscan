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

/* ARIADNE V.1.0

 GAPSTAT - functions for calculating the statistics of gapped
 alignment scores


 Copyright

 Richard Mott 2000

 Wellcome Trust Centre For Human Genetics
 Roosevelt Drive
 Oxford OX3 7AD

 Modified Sep 15, 2004 by Victor Joukov
 adapted to psimscan's matrix and frequency
 parameter convention


*/

// Modified Sep 15, 2004 by Victor Joukov
// adapted to psimscan's matrix and frequency
// parameter convention
#ifndef __gapstat_h__
#define __gapstat_h__

#include "gaplibshim.h"


void init_lambda( double lambda );
double e2lk( int k);


int KarlinAltschulStatistics(WMatrixType &matrix, ResFreqType &freq1, ResFreqType &freq2, double *lambda, double *Kminus, double *Kplus, double *H, double *r, double *s ) ;

double lambda_func( double lambda, double *h, int hmin, int hmax);

double solve_for_lambda( double *h, int hmin, int hmax  );

double entropy_H( double *h, int hmin, int hmax, double lambda ) ;

double HSP_length_correction( double H, double K, double lambda, int len1, int len2);

void iglehart( double *h, int hmin, int hmax, double lambda, double xmu, int delta, double *R1, double *R2, double *R3, double *Kplus, double *Kminus );

// double *get_h( int **matrix, double *freq1, double *freq2, int *hmin, int *hmax, double *mean );

double associated_dam_eqn( double *distribution, int min, int max, double lambda, double **transient, int *cn);

#ifdef USE_GEM_STATISTICS

    int gem_statistics( double gap_start, double gap_extend, double lambda, double s, double Kplus, double Kminus, double *alpha, double *theta1, double *theta2, double *K1, double *K2 );
    double upper_objective( double alpha, double beta );
    double lower_objective( double alpha, double beta );
    double gapped_theta( double alpha, double (*objective)( double, double) );

#endif

double GEM_length_correction( double L_HSP, double alpha, int len1, int len2 );

double regression_coeffs( double lambda, double Kmin, double H, double alpha, int len1, int len2, double *theta, double *kappa, double score );

double EmpiricalGEM( double lambda0, double K0, double H, double alpha, int len1, int len2, double *theta, double *logkappa, double score );

#endif /* __gapstat_h__ */
