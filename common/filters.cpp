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

#include <stdlib.h>
#include <sciminmax.h>
#include "filters.h"
#include "biosequence.h"

/*
calculates nucleotide homology level in %
returned level is corrected by nucleotide composition probability
*/
float nu_score_b (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, int* matc)
{
    int i, len, matches;
    int xn[4], yn[4];
    float p;

    //reset nucleotide composition
    for (i = 0; i < 4; i++)
        xn[i] = 0, yn[i] = 0;

    //calculate nucleotide composition
    for (len = 0, matches = 0; b_cnt > 0; b_cnt--, len += b_ptr->len, b_ptr++)
        for (i = 0; i < b_ptr->len; i++)
        {
            int xb = get_base(xseq, b_ptr->xpos + i);
            int yb = get_base(yseq, b_ptr->ypos + i);
            if (xb == yb) matches ++;
            xn[xb]++, yn[yb]++;
        }

    //calculate single match probability
    for (i = 0, p = 0.; i < 4; i++)
        p += (float) xn[i] * yn[i];

    //return adjusted homology level in %
    // ((matches * 100) / len) / ((4. * p) / (len * len))
    if (matc)
        *matc = matches;

    return ((float) matches * 100 * len / (p * 4));
}


/*
calculates aminoacid homology level in % in 6 frames
normalizes similarity to Y sequence self-similarity level (100%)
returns highest scoring Y - sequence frame number
scorexy contains all 6 frame scores

this function works with single batch only
*/
int aa_score (char* xseq, int xpos, char* yseq, int ypos, int len, int* na_wt[], int* scorexy)
{
    int i, j, tnum;
    int px, py, txy, tyy;
    int frame = ypos % 3;
    int max_score = 0;
    int scoreyy [6];
    int *wt_f = na_wt[0];
    int *wt_r = na_wt[1];

    //reset scores
    for (i = 0; i < 6; i++)
    {
        scorexy [i] = 0;
        scoreyy [i] = 0;
    }

    //calculate Y-X cross-similarity scores	and Y-Y self-similarity scores
    for (i = 0; i < len - 2; i += 3)
    {
        //get next 12 bases (good for up to 12 - 3 = 9 triplets)
        px = get_12_bases (xseq, xpos + i);
        py = get_12_bases (yseq, ypos + i);

        //calculate number of triplets that are inside the batch
        tnum = min_ (9, len - 2 - i);
        for (j = 0; j < tnum; j++)
        {
            //extract triplets
            txy = ((py << 6) & 0xfc0) | (px & 0x3f);
            tyy = ((py << 6) & 0xfc0) | (py & 0x3f);

            //calculate Y-X cross-similarity scores
            scorexy [frame] 	+= wt_f [txy];	//forward frames
            scorexy [frame + 3] += wt_r [txy];	//reverse frames

            //calculate Y-Y self-similarity scores
            scoreyy [frame] 	+= wt_f [tyy];	//forward frames
            scoreyy [frame + 3] += wt_r [tyy];	//reverse frames

            if (++frame > 2) frame = 0;
            px >>= 2, py >>= 2;
        }
    }


    //normalise X-Y scores with respect to Y-Y scores, convert to %
    for (i = 0; i < 6; i++)
        scorexy [i] = scorexy [i] * 100 / scoreyy [i];

    //select the highest scoring frame
    for (i = 0; i < 6; i++)
        if (scorexy[i] > max_score)
            max_score = scorexy[i], frame = i;

    return frame;
}

/*
calculates protein alignment score on the protein similarity batches
also calculates auto-similarity scores on x and y fragment sets comprising similarity
*/

int prot_score (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, WEIGHTS<int, 24>* wm, int* x_auto, int* y_auto)
{
    int pscore = 0;
    int xscore = 0;
    int yscore = 0;
    int prevxend = -1;
    int prevyend = -1;
    int gip = int (wm->gip);
    int gep = int (wm->gep);


    for (int idx = 0; idx < b_cnt; idx ++)
    {
        BATCH& batch = b_ptr [idx];

        // calc gaps
        if (prevxend != -1)
        {
            if (prevxend < batch.xpos)
                pscore -= (gip + (batch.xpos - prevxend) * gep);
            if (prevyend < batch.ypos)
                pscore -= (gip + (batch.ypos - prevyend) * gep);
        }
        // calc matching scores
        int xaa, yaa;
        for (int bpos = 0; bpos < batch.len; bpos ++)
        {
            xaa = xseq [bpos + batch.xpos];
            yaa = yseq [bpos + batch.ypos];

            xscore += wm->mx [xaa][xaa];
            yscore += wm->mx [yaa][yaa];
            pscore += wm->mx [yaa][xaa];
        }

        prevxend = batch.xpos + batch.len;
        prevyend = batch.ypos + batch.len;
    }
    if (x_auto) *x_auto = xscore;
    if (y_auto) *y_auto = yscore;
    return pscore;
}



/*
calculates protein alignment score on the protein-to-nucleotide similarity batches
also calculates auto-similarity scores on x and y fragment sets comprising similarity
*/
int an_score (char* xaaseq, char* ynnseq, BATCH* b_ptr, int b_cnt, bool aafirst, WEIGHTS<int, 24>* wm, int* x_auto, int* y_auto)
{
    // xaaseq is AA, ynnseq is NN

    // in the batch, x and y are switched

    int pscore = 0;
    int xscore = 0;
    int yscore = 0;
    int prevxend = -1;
    int prevyend = -1;
    int bxpos, bypos;
    int gip = int (wm->gip);
    int gep = int (wm->gep);

    for (int idx = 0; idx < b_cnt; idx ++)
    {
        BATCH& batch = b_ptr [idx];
        if (aafirst) bxpos = batch.xpos, bypos = batch.ypos;
        else bxpos = batch.ypos, bypos = batch.xpos;

        // calc gaps
        if (prevxend != -1)
        {
            if (prevxend < bxpos)
                pscore -= (gip + (bxpos - prevxend) * gep);

            if (prevyend < bypos)
                pscore -= (gip + (bypos - prevyend) * gep);

        }
        // calc matching scores
        int xaa, yaa;
        for (int bpos = 0; bpos < batch.len; bpos ++)
        {
            xaa = xaaseq [bpos + bxpos];

            int ycodone = (get_12_bases (ynnseq, bpos*3 + bypos)) & (0x3f);
            yaa = gcode_b [ycodone];

            xscore += wm->mx [xaa][xaa];
            yscore += wm->mx [yaa][yaa];
            pscore += wm->mx [yaa][xaa];
        }

        prevxend = bxpos + batch.len;
        prevyend = bypos + batch.len * 3;
    }
    if (x_auto) *x_auto = xscore;
    if (y_auto) *y_auto = yscore;
    return pscore;
}




/*
calculates chi^2 statistics for nucleotide similarity batch based on following assumptions:
    1. event space - 5 (a-match, g-match, c-match, t-match, mismatch), degrees of freedom - 4
    2. local (calculated on a batch) nucletide frequencies are used as reference
    (local composition bias is not considered meaningful by this hypothesis)
chi^2 values can be directly compared as long as both chi^2 values > chi^2 distribution mode
(pass1/pass2 selection heuristics almost guarantees that condition)
*/
float nu_matches_chi2_b (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt)
{
    int i, len;
    int xn[4], yn[4];
    int evt[5];
    float exp[5];
    float diff, chi2;

    //reset nucleotide composition and event count
    for (i = 0; i < 4; i++)
        xn[i] = 0, yn[i] = 0, evt[i] = 0;

    //calculate nucleotide composition and match count
    for (len = 0, evt[4] = 0; b_cnt > 0; b_cnt--, len += b_ptr->len, b_ptr++)
        for (i = 0; i < b_ptr->len; i++)
        {
            int xb = get_base(xseq, b_ptr->xpos + i);
            int yb = get_base(yseq, b_ptr->ypos + i);
            xn[xb]++, yn[yb]++;
            if (xb == yb) evt[xb]++, evt[4]--;
        }
    evt[4] += len;

    if (len == 0) return 0.;

    //calculate expected event counts based on nucleotide composition
    for (i = 0, exp[4] = (float) len; i < 4; i++)
    {
        exp[i] = (float) xn[i] * yn[i] / len;
        exp[4] -= exp[i];
    }

    //calculate chi^2
    //direct form: chi2 = sum{1..k} ((Ys - nPs) ^ 2 / nPs)
    for (i = 0, chi2 = 0.; i < 5; i++)
        if (exp[i] > 0.)
        {
            diff = exp[i] - evt[i];
            chi2 += diff * diff / exp[i];
        }

    return chi2;
}


/*
calculates chi^2 statistics for nucleotide similarity batch based on following assumptions:
    1. event space - 16 (aa, ag, ac, at, ga, gg, ... tt), degrees of freedom - 15 (or 9)
    2. local (calculated on a batch) nucletide frequencies are used as reference
    (local composition bias is not considered meaningful by this hypothesis)
chi^2 values can be directly compared as long as both chi^2 values > chi^2 distribution mode
(pass1/pass2 selection heuristics almost guarantees that condition)
*/
float nu_all_chi2_b (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt)
{
    int i, j, len;
    int xn[4], yn[4];
    int evt[4][4];
    float exp[4][4];
    float diff, chi2;

    //reset nucleotide composition and event count
    for (i = 0; i < 4; i++)
    {
        xn[i] = 0, yn[i] = 0;
        for (j = 0; j < 4; j++)
            evt[i][j] = 0;
    }

    //calculate nucleotide composition and match count
    for (len = 0; b_cnt > 0; b_cnt--, len += b_ptr->len, b_ptr++)
        for (i = 0; i < b_ptr->len; i++)
        {
            int xb = get_base(xseq, b_ptr->xpos + i);
            int yb = get_base(yseq, b_ptr->ypos + i);
            xn[xb]++, yn[yb]++, evt[xb][yb]++;
        }

    if (len == 0) return 0.;

    //calculate expected event counts based on nucleotide composition
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            exp[i][j] = ((float) xn[i] * yn[j]) / len;

    //calculate chi^2
    //direct form: chi2 = sum{1..k} ((Ys - nPs) ^ 2 / nPs)
    for (i = 0, chi2 = 0.; i < 4; i++)
        for (j = 0; j < 4; j++)
            if (exp[i][j] > 0.)
            {
                diff = exp[i][j] - evt[i][j];
                chi2 += diff * diff / exp[i][j];
            }

    return chi2;
}


/*
calculates an array of similarity scores (matches) with variable offset
between X and Y sequences (neighbor diagonals), 0 < offs < max_offs
returns offs != 0 with maximum score
*/
int nu_offs_score_b (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, int max_offs, int* offs_scores)
{
    int period;
    int max_score;
    int i, j, score;
    int rem;

    max_score = 0;
    for (int offs = 0; offs < max_offs; offs++)
    {
        //accumulate matches, step through 12 bases at a time
        score = 0;
        for (j = 0; j < b_cnt; j++)
            for (i = 0;;)
            {
                int r = ~get_12_bases (xseq, b_ptr[j].xpos + i) ^ get_12_bases (yseq, b_ptr[j].ypos + offs + i);
                r &= r >> 1;

                i += 12;
                rem = i - (b_ptr[j].len - offs);
                if (rem >= 0)
                {
                    //boundary case, shift out unused comparison bits
                    score += count_12 ((r << (rem << 1)) & 0x555555);
                    break;
                }
                score += count_12 (r & 0x555555);
            }
        offs_scores [offs] = score;
        if (score > max_score && offs > 0)
            max_score = score, period = offs;
    }
    return period;
}


/*
calculates modulo-3 positional mismatch statistics for each y sequence
(does not initialise scores)
*/
int nu_mod3_score_accum_b  (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, int* y_mod3_scores)
{
    int i;

    for (; b_cnt > 0; b_cnt--, b_ptr++)
        for (i = 0; i < b_ptr->len; i++)
            if (get_base (xseq, b_ptr->xpos + i) != get_base (yseq, b_ptr->ypos + i))
                y_mod3_scores [(b_ptr->ypos + i) % 3]++;

    return 0;
}

/*
calculates full length positional mismatch statistics for each y sequence
(does not initialise scores)
*/
int nu_pos_score_accum_b  (char* xseq, char* yseq, BATCH* b_ptr, int b_cnt, int* mism, int* matc, bool revers, int ylen)
{
    int i;

    for (; b_cnt > 0; b_cnt--, b_ptr++)
        for (i = 0; i < b_ptr->len; i++)
        {
            if (get_base (xseq, b_ptr->xpos + i) != get_base (yseq, b_ptr->ypos + i))
            {
                if (!revers)
                    mism [b_ptr->ypos + i] ++;
                else
                    mism [ylen - 1 - (b_ptr->ypos + i)] ++;
            }
            else
            {
                if (!revers)
                    matc [b_ptr->ypos + i] ++;
                else
                    matc [ylen - 1 - (b_ptr->ypos + i)] ++;
            }

        }

    return 0;
}

//calculates nucleotide homology level in %
//returned level is corrected by nucleotide composition probability
int nu_score (char* xseq, int xpos, char* yseq, int ypos, int len, int* matc)
{
    int i, matches;
    int xn[4], yn[4];
    float p;

    //reset nucleotide composition
    for (i = 0; i < 4; i++)
        xn[i] = 0, yn[i] = 0;

    //calculate nucleotide composition
    for (i = 0, matches = 0; i < len; i++)
    {
        int xb = get_base(xseq, xpos + i);
        int yb = get_base(yseq, ypos + i);
        if (xb == yb) matches ++;
        xn[xb]++, yn[yb]++;
    }

    //calculate single match probability
    for (i = 0, p = 0.; i < 4; i++)
        p += (float) xn[i] * yn[i];

    //return adjusted homology level in %
    // ((matches * 100) / len) / ((4. * p) / (len * len))
    if (matc)
        *matc = matches;

    return (int)((float) matches * 100. * len / (4. * p));
}


/*
calculates chi^2 statistics for nucleotide similarity batch based on following assumptions:
    1. event space - 5 (a-match, g-match, c-match, t-match, mismatch), degrees of freedom - 4
    2. local (calculated on a batch) nucletide frequencies are used as reference
    (local composition bias is not considered meaningful by this hypothesis)
chi^2 values can be directly compared as long as both chi^2 values > chi^2 distribution mode
(pass1/pass2 selection heuristics almost guarantees that condition)
*/
float nu_matches_chi2 (char* xseq, int xpos, char* yseq, int ypos, int len)
{
    int i;
    int xn[4], yn[4];
    int evt[5];
    float exp[5];
    float diff, chi2;

    if (len == 0) return 0.;

    //reset nucleotide composition and event count
    for (i = 0; i < 4; i++)
        xn[i] = 0, yn[i] = 0, evt[i] = 0;

    //calculate nucleotide composition and match count
    for (i = 0, evt[4] = len; i < len; i++)
    {
        int xb = get_base(xseq, xpos + i);
        int yb = get_base(yseq, ypos + i);
        xn[xb]++, yn[yb]++;
        if (xb == yb) evt[xb]++, evt[4]--;
    }

    //calculate expected event counts based on nucleotide composition
    for (i = 0, exp[4] = (float) len; i < 4; i++)
    {
        exp[i] = (float) xn[i] * yn[i] / len;
        exp[4] -= exp[i];
    }

    //calculate chi^2
    //direct form: chi2 = sum{1..k} ((Ys - nPs) ^ 2 / nPs)
    for (i = 0, chi2 = 0.; i < 5; i++)
        if (exp[i] > 0.)
        {
            diff = exp[i] - evt[i];
            chi2 += diff * diff / exp[i];
        }

    return chi2;
}


/*
calculates chi^2 statistics for nucleotide similarity batch based on following assumptions:
    1. event space - 16 (aa, ag, ac, at, ga, gg, ... tt), degrees of freedom - 15 (or 9)
    2. local (calculated on a batch) nucletide frequencies are used as reference
    (local composition bias is not considered meaningful by this hypothesis)
chi^2 values can be directly compared as long as both chi^2 values > chi^2 distribution mode
(pass1/pass2 selection heuristics almost guarantees that condition)
*/
float nu_all_chi2 (char* xseq, int xpos, char* yseq, int ypos, int len)
{
    int i, j;
    int xn[4], yn[4];
    int evt[4][4];
    float exp[4][4];
    float diff, chi2;

    if (len == 0) return 0.;

    //reset nucleotide composition and event count
    for (i = 0; i < 4; i++)
    {
        xn[i] = 0, yn[i] = 0;
        for (j = 0; j < 4; j++)
            evt[i][j] = 0;
    }

    //calculate nucleotide composition and match count
    for (i = 0; i < len; i++)
    {
        int xb = get_base(xseq, xpos + i);
        int yb = get_base(yseq, ypos + i);
        xn[xb]++, yn[yb]++, evt[xb][yb]++;
    }

    //calculate expected event counts based on nucleotide composition
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            exp[i][j] = (float) xn[i] * yn[j] / len;

    //calculate chi^2
    //direct form: chi2 = sum{1..k} ((Ys - nPs) ^ 2 / nPs)
    for (i = 0, chi2 = 0.; i < 4; i++)
        for (j = 0; j < 4; j++)
            if (exp[i][j] > 0.)
            {
                diff = exp[i][j] - evt[i][j];
                chi2 += diff * diff / exp[i][j];
            }

    return chi2;
}


/*
calculates an array of similarity scores (matches) with variable offset
between X and Y sequences (neighbor diagonals), 0 < offs < max_offs
returns offs != 0 with maximum score
*/
int nu_offs_score (char* xseq, int xpos, char* yseq, int ypos, int len, int max_offs, int* offs_scores)
{
    int period;
    int max_score = 0;
    int i, score;
    int rem;

    for (int offs = 0; offs < max_offs; offs++)
    {
        //accumulate matches, step through 12 bases at a time
        for (i = 0, score = 0;;)
        {
            int r = ~get_12_bases (xseq, xpos + i) ^ get_12_bases (yseq, ypos + offs + i);
            r &= r >> 1;

            i += 12;
            if ((rem = i - (len - offs)) >= 0)
            {
                //boundary case, shift out unused comparison bits
                score += count_12 ((r << (rem << 1)) & 0x555555);
                break;
            }
            score += count_12 (r & 0x555555);
        }

        offs_scores [offs] = score;
        if (score > max_score && offs > 0)
            max_score = score, period = offs;
    }

    return period;
}
