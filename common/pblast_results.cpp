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

#include "pblast_results.h"
#include "gapstat.h"
#include <math.h>


SimMergerBase PblastResults :: null_merger_;

PblastResults::PblastResults (double min_score, int min_len, WMatrixType *w, bool eval_eval, unsigned res_per_qry, SimMergerBase& merger)
:
MergingResultStorage (merger, res_per_qry),
total_found_ (0),
min_score_ (min_score),
min_len_ (min_len),
w_(w),
uid_cached_(-1),
eval_eval_ (eval_eval)
{
}

const double ln2 = 0.69314718056;

// adapted from Ariadne's seq_util.c

static void RobinsonResidueFrequencies(ResFreqType &freq)
{
    /* these are are amino-acid frequencies used by blast, from Robinson & Robinson */
    int n;
    for (n = 0; n < AANUM; ++n) freq[n] = 0.0;

    freq[char2aa('A')] = 35155;
    freq[char2aa('R')] = 23105;
    freq[char2aa('N')] = 20212;
    freq[char2aa('D')] = 24161;
    freq[char2aa('C')] = 8669;
    freq[char2aa('Q')] = 19208;
    freq[char2aa('E')] = 28354;
    freq[char2aa('G')] = 33229;
    freq[char2aa('H')] = 9906;
    freq[char2aa('I')] = 23161;
    freq[char2aa('L')] = 40625;
    freq[char2aa('K')] = 25872;
    freq[char2aa('M')] = 10101;
    freq[char2aa('F')] = 17367;
    freq[char2aa('P')] = 23435;
    freq[char2aa('S')] = 32070;
    freq[char2aa('T')] = 26311;
    freq[char2aa('W')] = 5990;
    freq[char2aa('Y')] = 14488;
    freq[char2aa('V')] = 29012;

    // Normalize
    double t = 0.0;
    for (n = 0; n < AANUM; ++n) t += freq[n];
    for (n = 0; n < AANUM; ++n) freq[n] /= t;
}

static void PseudoResidueFrequencies(SEQ &seq, int pseudo_len, ResFreqType &freq)
{
    int n;

    if (seq.len > 0) {
        RobinsonResidueFrequencies(freq);

        // It comes normalized to 1, renormalize it to pseudo_len
        for (n = 0; n < AANUM; ++n) freq[n] *= pseudo_len;
        // Add real frequencies
        for (n = 0; n < seq.len; ++n) freq[seq.seq[n]]++;

        // Renormalize back to 1
        double t = 0.0;
        for (n = 0; n < AANUM; ++n) t += freq[n];
        for (n = 0; n < AANUM; ++n) freq[n] /= t;
    } else {
        for (n = 0; n < AANUM; ++n) freq[n] = 0.0;
    }
}


bool PblastResults::match_found  (SEQ& query_seq, SEQ& target_seq, BATCH* batches, int batch_no, double score, double q_auto_score, double t_auto_score)
{
    total_found_ ++;

    if (score < min_score_)
        return false;

    int len = 0;
    for (int b = 0; b < batch_no; b ++)
        len += batches [b].len;
    if (len < min_len_)
        return false;

    double pvalue, bitscore;

    if (eval_eval_) {
        // Following is copied over (with adaptation) from Ariadne's file ariadne.c

        // Mott's gap penalties, A+B*k
        double A = w_->gip-w_->gep;
        double B = w_->gep;
        if (target_seq.uid != uid_cached_) {
            PseudoResidueFrequencies(target_seq, 100, trg_freq_);
            ResFreqType freq_mean;
            RobinsonResidueFrequencies(freq_mean);
            double Kplus; // necessary output parameter, never used in this algorithm
            double r, s;
            K0mean_ = lambda0mean_ = Hmean_ = alpha_mean_ = 0.0;
            if (KarlinAltschulStatistics(*w_, trg_freq_, freq_mean, &lambda0mean_, &K0mean_, &Kplus, &Hmean_, &r, &s)) {
                alpha_mean_ = 2*s*exp(-lambda0mean_*(A+B))/(1-exp(-lambda0mean_*B));
            }
            uid_cached_ = target_seq.uid;
        }

        double lambda0, K0, Kplus, H, r, s, alpha;
        double theta, logkappa;
        double pval1=1.0, pval2=1.0;
        // TODO: may be better return e-value (before 1-exp(-x)) from EmpiricalGEM???
        pval1 = pval2 = EmpiricalGEM(lambda0mean_, K0mean_, Hmean_, alpha_mean_, query_seq.len, target_seq.len, &theta, &logkappa, score);
	    // TODO: real threshold
	    if (pval1 < 0.1 /* params->pthresh1 */) {
	      ResFreqType freq;
          PseudoResidueFrequencies(query_seq, 100, freq);
	      if (KarlinAltschulStatistics(*w_, trg_freq_, freq, &lambda0, &K0, &Kplus, &H, &r, &s)) {
	        alpha = 2*s*exp(-lambda0*(A+B))/(1-exp(-lambda0*B));
            // TODO: report invalid match here
    //	    if ( alpha > 0.25 )
    //	      IterativelyIncreaseGapPenalty( seq, pro, &A, &B, lambda0, K0, H, s, &alpha, &sw_score );
	        pval2 = EmpiricalGEM(lambda0, K0, H, alpha, query_seq.len, target_seq.len, &theta, &logkappa, score);
	      }
	      else {
    //	    printf("ERROR sequence %s has positive expected match score with profile %s - ignored\n", seq->name, pro->name );
	        pval2 = 1.0;
	      }
	    }
	    else {
	      pval2 = pval1;
	    }
    // TODO: handle this threshold instead of sw cutoff
    //	if (pval2 < params->pthresh2 ) return false;

        // for bitscore
        double logKappa = log(K0) + logkappa;
        double Lambda = lambda0*theta;
        bitscore = (Lambda * score - logKappa) / ln2;
        pvalue = pval2;
    } else {
        // Mock code with fixed K and lambda - use ONLY for benchmarking
        const double K = 0.041;
        const double lambda = 0.267;
        const double logK = -3.1941832123;

        bitscore = (lambda * score - logK) / ln2;
        double x = K*query_seq.len*target_seq.len*exp(-lambda*score);

        if (x < 1.0e-9) pvalue = x;
        else pvalue = 1.0-exp(-x);
    }

	// evaluate identity
	int al_length = 0;
	int mismatches = 0;
	int xp, yp;
    for (int i = 0; i < batch_no; i ++)
    {
        xp = batches [i].xpos;
        yp = batches [i].ypos;
        for (int p = 0; p < batches [i].len; p ++)
        {
            if (query_seq.seq [xp] != target_seq.seq [yp])
                mismatches ++;
            xp ++;
            yp ++;
        }
        al_length += batches [i].len;
    }
    float p_identity = (float (al_length) - float (mismatches)) * (float) 100.0 / float (al_length);


    return add_result (query_seq.uid, target_seq.uid, false, score, (int) score, p_identity, pvalue, bitscore, (int) q_auto_score, (int) t_auto_score, batch_no, batches, NULL, 0);
}

int PblastResults::totalFound ()
{
    return total_found_;
}

