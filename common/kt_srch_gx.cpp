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
#include <string.h>
#include "kt_srch_gx.h"
#include <rerror.h>
#include <iostream>

#include <sciminmax.h>

#ifndef __min
#define __min min_
#endif
#ifndef __max
#define __max max_
#endif

static bool verbose = false;


KT_SEARCH::KT_SEARCH (int kt_sz, int *wts, int max_ynum, int max_total_ylen, int max_xlen, ResultReciever_blast* res_rec)
:
ki (NULL),
yi (NULL),
xseq (NULL),
di (NULL),
ke (NULL)
{
    batches_out = 0;
    batch = res_rec;
    max_x_len = max_xlen;
    max_y_num = max_ynum;
    max_y_len = max_total_ylen;
    kt_size = kt_sz;

    try
    {
        init_vars (max_ynum, wts);
    }
    catch (std::bad_alloc& exc)
    {
        ERR ("Insufficient memory for hash-based searcher");
    }
}

KT_SEARCH::KT_SEARCH (int kt_sz, int *wts, int max_ynum, int max_total_ylen, int max_xlen, ResultReciever_blast_batch* res_rec)
:
ki (NULL),
yi (NULL),
xseq (NULL),
di (NULL),
ke (NULL)
{
    batches_out = 1;
    batch_processor = res_rec;
    max_x_len = max_xlen;
    max_y_num = max_ynum;
    max_y_len = max_total_ylen;
    kt_size = kt_sz;

    init_vars (max_ynum, wts);
}


KT_SEARCH::~KT_SEARCH ()
{
    reset_y ();

    if (di)
        delete [] di;
    di = NULL;
    if (ki)
        delete [] ki;
    ki = NULL;
    if (yi)
        delete [] yi;
    yi = NULL;
    if (xseq)
        delete [] (xseq - 4);
    xseq = NULL;
}

void KT_SEARCH::init_vars (int max_ynum, int *wts)
{
    int i;
    // int ktuple info parameters
    kt_mask = 0xffffffff >> ((16 - kt_size) * 2);
    ki_len = 1 << (kt_size * 2);

    //init default search parameters
    max_offs = 0;
    rep_del = 0;
    // fast_mode = 0;
    xstep = 1;
    apprx_match = 0;
    k_thresh = 150;
    k_gip = 50;
    k_gep = 25;
    b_thresh_s = 90;
    b_thresh_l = 55;
    l_thresh = 30;
    g_period = 30;

    //init search variables
    ke_len = 0;
    ke = NULL;
    di_len = 0;
    di = NULL;
    kt_last = 0;
    y_tnum = 0;
    y_tlen = 0;
    k_matches = 0;
    b_matches = 0;

    //set status
    armed = 0;


    //allocate Y sequence info array
    yi = NULL;
    try
    {
        yi = new SEQ_INFO [max_ynum];
    }
    catch (std::bad_alloc&)
    {
    }
    if (!yi)
        ERR("unable to allocate Y sequence info array");

    // init seq_info array
    for (i = 0; i < max_ynum; i ++)
    {
        yi [i].seq = NULL;
        yi [i].uid = -1L;
    }


    //allocate k-tuple info array
    ki = NULL;
    try
    {
        ki = new KTUP_INFO [ki_len];
    }
    catch (std::bad_alloc&)
    {
    }
    if (!ki)
        ERR("unable to allocate k-tuple info array");

    //init k-tuple info array
    for (i = 0; i < ki_len; i++)
    {
        if (wts)
        ki[i].wt = wts [i];
        else
        ki[i].wt = DEF_K_WT;
        ki[i].cnt = 0;
        ki[i].ptr = 0;
    }

    //allocate x sequence (padd) buffer
    xseq = NULL;
    try
    {
        xseq = new char [((max_x_len+3)>>2)+8];
    }
    catch (std::bad_alloc&)
    {
    }
    if (!xseq)
        ERR ("unable to allocate x sequence buffer");
    xseq = xseq + 4;

    //allocate diag info array
    di = NULL;
    try
    {
        di = new DIAG_ENTRY [max_y_len + max_x_len + 1];
    }
    catch (std::bad_alloc)
    {
    }
    if (!di)
        ers << "Insufficient memory for diagonal info array, requested " << max_y_len + max_x_len + 1 << " (" << max_x_len << " max_subj_len + " << max_y_len << " max_qry_len + 1) DIAG_ENTRY structures (" << sizeof (DIAG_ENTRY) << " bytes each, " << (max_y_len + max_x_len + 1)*sizeof (DIAG_ENTRY) / (1024*1024) <<  " Mbytes total)" << Throw;


    // _ktups_found = 0;
    // _diags_found = 0;
    // _nuc_scanned = 0;
}

int KT_SEARCH::add_yseq (NN_SEQ& yseq)
{
    int bytes = (yseq.len + 3) >> 2;
    char *ptr;

    if (y_tnum >= max_y_num) ERR("number of Y sequences exceeded __max");

    yi[y_tnum].start = y_tlen;
    yi[y_tnum].len = yseq.len;
    yi[y_tnum].uid = yseq.uid;
    yi[y_tnum].rev = yseq.rev;

    //word padd sequence data for safe unaligned access
    ptr = NULL;
    try
    {
        ptr = new char [bytes + 8];
    }
    catch (std::bad_alloc&)
    {
    }
    if (!ptr)
        ERR ("unable to allocate Y sequence data");

    yi[y_tnum].seq = ptr + 4;
    memcpy (yi[y_tnum].seq, yseq.seq, bytes);

    y_tlen += yseq.len;
    if (y_tlen > max_y_len)
        ERR ("Maximum total length of Y sequences exceeded");
    y_tnum++;

    return 1;
}


int KT_SEARCH::init ()
{
    int i;

    //avoid multiple initialisations
    if (armed)
        return 0;
    armed = 1;

    //calculate batch delimiter parameters
    penalty = (l_thresh * (b_thresh_s - b_thresh_l)) / 100;
    offs_lev = __max (0, b_thresh_l - 5);
    // stop_lev = 3 * offs_lev;
    stop_lev = 3 * (100 - b_thresh_l);

    if (verbose) std::clog << "counting K-tuples ..." << std::endl << std::flush;

    for (i = 0; i < y_tnum; i++)
        ki_count (yi[i].seq, yi[i].len);

    if (verbose) std::clog << "Y seq K-tuple statistics: ki = " << ki_len << ", ke = " << ke_len << ", avg ke/ki = " << (double) ke_len / (double) ki_len << std::endl << std::flush;

    //create k-tuple entry array
    //allocate k-tuple entry array
    ke = NULL;
    try
    {
        ke = new KTUP_ENTRY [ke_len];
    }
    catch (std::bad_alloc&)
    {
    }
    if (!ke)
        ERR("unable to allocate k-tuple entry array");

    //init ki[].ptr, reset counters
    KTUP_ENTRY* ptr = ke;
    for (i = 0; i < ki_len; i++)
    {
        ki[i].ptr = ptr;
        ptr += ki[i].cnt;
        ki[i].cnt = 0;
    }

    //add k-tuples
    if (verbose) std::clog << "adding K-tuples ..." << std::endl;

    for (i = 0; i < y_tnum; i++)
        ke_add (yi[i].seq, yi[i].len, i, yi[i].start);

    //initialise diag info array
    di_len = y_tlen + max_x_len + 1;

    for (i = 0; i < di_len; i++)
    {
        di[i].xuid = -1L;
        di[i].ynum = -1;
        di[i].xpos = 0;
        di[i].score = 0;
        di[i].b_end = 0;
    }
    if (verbose) std::clog << std::endl << std::flush;
    return 1;
}


void KT_SEARCH::reset_y ()
{
    if (!armed) return;
    armed = 0;
    if (verbose) std::clog << "Resetting y" << std::endl << std::flush;
    //reset pointers to k-tuple entry array
    if (ki)
    {
        for (int i = 0; i < ki_len; i++)
        {
            ki[i].cnt = 0;
            ki[i].ptr = NULL;
        }
    }
    //deallocate k-tuple entry array
    if (ke)
        delete [] ke;
    ke = NULL;
    //if reset called after init deallocate yi array
    if (yi)
    {
        for (int i = 0; i < y_tnum; i++)
        {
            if (yi[i].seq)
            {
                delete [] (yi[i].seq - 4);
                yi[i].seq = NULL;
            }
        }
    }

    //init search variables
    ke_len = 0;
    ke = NULL;
    di_len = 0;
    kt_last = 0;
    y_tnum = 0;
    y_tlen = 0;
    k_matches = 0;
    b_matches = 0;

    //_ktups_found = 0;
    //_diags_found = 0;
    //_nuc_scanned = 0;
}


void KT_SEARCH::ki_count (char* seq, int len)
{
    int i, j, pos;
    unsigned ktup, ktup_m, mask, skip;
    //count k-tuples
    for (pos = 0; pos <= len - kt_size; pos++)
    {
        ktup = get_ktup (seq, pos, kt_mask);
        //filter repeating K-tuples
        prev_ktup[kt_last] = ktup, kt_last = (kt_last + 1) & 0x7;
        for (i = 0, skip = 0; i < rep_del; i++)
        if (prev_ktup[(kt_last - i - 2) & 0x7] == ktup && i < pos)
        {
            skip = 1;
            break;
        }
        if (skip) continue;
        ki[ktup].cnt++;
        ke_len++;
        if (apprx_match)
        {
            //if approximate matching is enabled run through kt_size*4 possible mutations
            for (i = 0; i < (kt_size << 1); i += 2)
            {
                mask = ~(3 << i);
                for (j = 0; j < 4; j++)
                {
                    ktup_m = (ktup & mask) | (j << i);
                    if (ktup_m == ktup) continue;
                    ki[ktup_m].cnt++;
                    ke_len++;
                }
            }
        }
    }
}


void KT_SEARCH::ke_add (char* seq, int len, int ynum, int ystart)
{
    int i, j, pos, ypos, k_cnt;
    unsigned ktup, ktup_m, mask, skip;

    //initialise k-tuple entries
    for (pos = 0; pos <= len - kt_size; pos++)
    {
        ktup = get_ktup (seq, pos, kt_mask);
        ypos = ystart + pos + max_x_len;

        //filter repeating K-tuples
        prev_ktup[kt_last] = ktup, kt_last = (kt_last + 1) & 0x7;
        for (i = 0, skip = 0; i < rep_del; i++)
        if (prev_ktup[(kt_last - i - 2) & 0x7] == ktup && i < pos)
        {
            skip = 1;
            break;
        }
        if (skip) continue;

        //init entry
        k_cnt = ki[ktup].cnt;
        ki[ktup].ptr[k_cnt].ynum = ynum;
        ki[ktup].ptr[k_cnt].ypos = ypos;
        ki[ktup].cnt++;
        if (apprx_match)
        {
            //if approximate matching is enabled run through kt_size*4 possible mutations
            for (i = 0; i < (kt_size << 1); i += 2)
            {
                mask = ~(3 << i);
                for (j = 0; j < 4; j++)
                {
                    ktup_m = (ktup & mask) | (j << i);

                    if (ktup_m == ktup) continue;

                    //init entry
                    k_cnt = ki[ktup_m].cnt;
                    ki[ktup_m].ptr[k_cnt].ynum = ynum;
                    ki[ktup_m].ptr[k_cnt].ypos = ypos;
                    ki[ktup_m].cnt++;
                }
            }
        }
    }
}


int KT_SEARCH::search (NN_SEQ& xseq)
{
    if (xseq.len > max_x_len)
    {
        std::clog << std::endl <<  "X sequence length " << xseq.len << " exceeds declared length " << max_x_len << std::flush;
        ERR("X sequence length exceeds declared length");
    }

    if (!armed) init ();

    int bytes = (xseq.len + 3) >> 2;
    memcpy (KT_SEARCH::xseq, xseq.seq, bytes);

    xlen = xseq.len;
    xuid = xseq.uid;
    xrev = xseq.rev;

    if (xstep == 4)
        scan_l1_fast ();
    else if (xstep == 1)
        scan_l1 ();
    else
        scan_l1_stepped ();
    return 1;
}


//pass 1
//K-tuple match processing
void KT_SEARCH::scan_l1 ()
{
    int maxx, dist, offs, pnl, score, k_wt, k_cnt;
    int xend = xlen - kt_size;
    char* s = xseq;
    unsigned seqw, ktup;
    KTUP_ENTRY* k_ptr;
    register DIAG_ENTRY* d_ptr;

    for (xpos = 0; xpos < xend; s++)
    {
        seqw = GET32_U (s);
        maxx = __min(xpos + 4, xend);
          //run up to 4 times (byte level loop)
        for (; xpos < maxx; xpos++)
        {
            //get next K-tuple
            ktup = seqw & kt_mask, seqw >>= 2;

            // _nuc_scanned ++;

            //filter repeating K-tuples
            kt_last = (kt_last + 1) & 0x7;
            prev_ktup[kt_last] = ktup;
            for (int i = 1; i <= __min (rep_del, xpos); i++)
                if (prev_ktup[(kt_last - i) & 0x7] == ktup)
                    goto SKIP_KTUP;

            //get K-tuple info
            k_cnt = ki[ktup].cnt;
            k_wt  = ki[ktup].wt;
            k_ptr = ki[ktup].ptr;

            // if (k_cnt)
            //    _ktups_found ++;

            //for all matched K-tuples
            while (k_cnt--)
            {
                ynum = k_ptr->ynum;
                ypos = k_ptr->ypos;
                k_ptr++;

                diag = ypos - xpos;
                d_ptr = di + diag;

                //look at previous match on the same diagonal
                if ((d_ptr->ynum == ynum) && (d_ptr->xuid == xuid) && (d_ptr->xrev == xrev))
                {
                    //previous match in same sequence - adjust score
                    dist = xpos - d_ptr->xpos - kt_size;
                    if (dist <= -kt_size)
                        //batch already identified - do not update diag info
                        continue;
                    if (dist < 0)
                        //overlapping k-tuple matches - score increases with distance
                        score = __max (0, d_ptr->score) + (k_wt >> (-dist >> 2));
                    else
                        //non-overlapping k-tuple matches - score decreases with distance
                        score = __max (k_wt, k_wt + d_ptr->score - dist);
                }
                else
                {
                    //previous match in different sequence - reinitialise score

                    d_ptr->ynum = ynum;
                    d_ptr->xuid = xuid;
                    d_ptr->xrev = xrev;
                    d_ptr->b_end = 0;
                    score = k_wt;
                }

                //look at previous matches on adjacent diagonals, pick the best score
                //check non-overlapping k-tuple matches only, penalise offset
                for (offs = 1, pnl = k_gip; offs <= max_offs; offs++, pnl += k_gep)
                {
                    //go left
                    register DIAG_ENTRY* ad_ptr = d_ptr + offs;
                    if ((ad_ptr->ynum == ynum) && (ad_ptr->xuid == xuid) && (ad_ptr->xrev == xrev) && ((dist = xpos - ad_ptr->xpos - kt_size - offs) >= 0))
                        score = __max (score, k_wt + ad_ptr->score - dist - pnl);

                    //go down
                    ad_ptr = d_ptr - offs;
                    if ((ad_ptr->ynum == ynum) && (ad_ptr->xuid == xuid) && (ad_ptr->xrev == xrev)&& ((dist = xpos - ad_ptr->xpos - kt_size) >= 0))
                        score = __max (score, k_wt + ad_ptr->score - dist - pnl);
                }

                d_ptr->xpos = xpos;
                d_ptr->score = score;

                if (score > k_thresh)
                {
                    //convert ypos global to ypos local
                    ypos -= yi[ynum].start + max_x_len;
                    yseq = yi[ynum].seq;
                    ylen = yi[ynum].len;

                    //delimit / filter batch
                    // _diags_found ++;
                    if (!batches_out)
                        scan_l2_g ();
                    else
                        scan_l2_gb_1 ();
                    k_matches++;
                }
            }
SKIP_KTUP:
            continue;
        }
    }
}


//alternative pass 1
//K-tuple match processing
//trades off sensitivity for increased speed
//(stepping through 4 bases (8 bits) at a time)
//use for high similarity level (>80%) search only !!!
void KT_SEARCH::scan_l1_fast ()
{
    int dist, offs, pnl, score, k_wt, k_cnt;
    int xend = xlen - kt_size;
    char* s = xseq;
    unsigned ktup;
    KTUP_ENTRY* k_ptr;
    register DIAG_ENTRY* d_ptr;

    for (xpos = 0; xpos < xend; s++, xpos += 4)
    {
        //get next K-tuple
        ktup = GET32_U (s) & kt_mask;

        //get K-tuple info
        k_cnt = ki[ktup].cnt;
        k_wt  = ki[ktup].wt;
        k_ptr = ki[ktup].ptr;

        //for all matched K-tuples
        while (k_cnt--)
        {
            ynum = k_ptr->ynum;
            ypos = k_ptr->ypos;
            k_ptr++;

            diag = ypos - xpos;
            d_ptr = di + diag;

            //look at previous match on the same diagonal
            if ((d_ptr->ynum == ynum) && (d_ptr->xuid == xuid) && (d_ptr->xrev == xrev))
            {
                //previous match in same sequence - adjust score
                dist = xpos - d_ptr->xpos - kt_size;
                if (dist <= -kt_size)
                //batch already identified - do not update diag info
                continue;
                if (dist < 0)
                //overlapping k-tuple matches - score increases with distance
                score = __max (0, d_ptr->score) + (k_wt >> (-dist >> 2));
                else
                //non-overlapping k-tuple matches - score decreases with distance
                score = __max (k_wt, k_wt + d_ptr->score - dist);
            }
            else
            {
                //previous match in different sequence - reinitialise score
                d_ptr->ynum = ynum;
                d_ptr->xuid = xuid;
                d_ptr->xrev = xrev;
                d_ptr->b_end = 0;
                score = k_wt;
            }

            //look at previous matches on adjacent diagonals, pick the best score
            //check non-overlapping k-tuple matches only, penalise offset
            for (offs = 1, pnl = k_gip; offs <= max_offs; offs++, pnl += k_gep)
            {
                //go left
                register DIAG_ENTRY* ad_ptr = d_ptr + offs;
                if ((ad_ptr->ynum == ynum) && (ad_ptr->xuid == xuid) && (ad_ptr->xrev == xrev) && ((dist = xpos - ad_ptr->xpos - kt_size - offs) >= 0))
                score = __max (score, k_wt + ad_ptr->score - dist - pnl);

                //go down
                ad_ptr = d_ptr - offs;
                if ((ad_ptr->ynum == ynum) && (ad_ptr->xuid == xuid) && (ad_ptr->xrev == xrev) && ((dist = xpos - ad_ptr->xpos - kt_size) >= 0))
                score = __max (score, k_wt + ad_ptr->score - dist - pnl);
            }

            d_ptr->xpos = xpos;
            d_ptr->score = score;

            if (score > k_thresh)
            {
                //convert ypos global to ypos local
                ypos -= yi[ynum].start + max_x_len;
                yseq = yi[ynum].seq;
                ylen = yi[ynum].len;

                //delimit / filter batch
                if (!batches_out)
                    scan_l2_g ();
                else
                    scan_l2_gb_1 ();
                k_matches++;
            }
        }
    }
}

//alternative pass 1
//K-tuple match processing
//trades off sensitivity for increased speed
//(stepping through subject sequence STEP bits at a time
void KT_SEARCH::scan_l1_stepped ()
{
    int dist, offs, pnl, score, k_wt, k_cnt;
    int xend = xlen - kt_size;
    char* s = xseq;
    unsigned ktup;
    KTUP_ENTRY* k_ptr;
    register DIAG_ENTRY* d_ptr;
    int bytebound = 0;
    int dword_reminder = 16 - kt_size;
    int shift = 0;
    int baseoff = 0;
    DWORD seqw = GET32_U (s);

    for (xpos = 0; xpos < xend; xpos += xstep)
    {
        // find if last boundary is still good
        while (xpos > ((bytebound + shift) << 2) + dword_reminder)
            ++shift;
        if (shift)
        {
            bytebound += shift;
            s += shift;
            shift = 0;
            baseoff = xpos - (bytebound << 2);
            seqw = GET32_U (s);
            seqw >>= (baseoff << 1);
        }

        //get next K-tuple
        ktup = seqw;
        ktup &= kt_mask;
        seqw >>= 2;


        //get K-tuple info
        k_cnt = ki[ktup].cnt;
        k_wt  = ki[ktup].wt;
        k_ptr = ki[ktup].ptr;

        //for all matched K-tuples
        while (k_cnt--)
        {
            ynum = k_ptr->ynum;
            ypos = k_ptr->ypos;
            k_ptr++;

            diag = ypos - xpos;
            d_ptr = di + diag;

            //look at previous match on the same diagonal
            if ((d_ptr->ynum == ynum) && (d_ptr->xuid == xuid) && (d_ptr->xrev == xrev))
            {
                //previous match in same sequence - adjust score
                dist = xpos - d_ptr->xpos - kt_size;
                if (dist <= -kt_size)
                //batch already identified - do not update diag info
                continue;
                if (dist < 0)
                //overlapping k-tuple matches - score increases with distance
                score = __max (0, d_ptr->score) + (k_wt >> (-dist >> 2));
                else
                //non-overlapping k-tuple matches - score decreases with distance
                score = __max (k_wt, k_wt + d_ptr->score - dist);
            }
            else
            {
                //previous match in different sequence - reinitialise score
                d_ptr->ynum = ynum;
                d_ptr->xuid = xuid;
                d_ptr->xrev = xrev;
                d_ptr->b_end = 0;
                score = k_wt;
            }

            //look at previous matches on adjacent diagonals, pick the best score
            //check non-overlapping k-tuple matches only, penalise offset
            for (offs = 1, pnl = k_gip; offs <= max_offs; offs++, pnl += k_gep)
            {
                //go left
                register DIAG_ENTRY* ad_ptr = d_ptr + offs;
                if ((ad_ptr->ynum == ynum) && (ad_ptr->xuid == xuid) && (ad_ptr->xrev == xrev) && ((dist = xpos - ad_ptr->xpos - kt_size - offs) >= 0))
                score = __max (score, k_wt + ad_ptr->score - dist - pnl);

                //go down
                ad_ptr = d_ptr - offs;
                if ((ad_ptr->ynum == ynum) && (ad_ptr->xuid == xuid) && (ad_ptr->xrev == xrev) && ((dist = xpos - ad_ptr->xpos - kt_size) >= 0))
                score = __max (score, k_wt + ad_ptr->score - dist - pnl);
            }

            d_ptr->xpos = xpos;
            d_ptr->score = score;

            if (score > k_thresh)
            {
                //convert ypos global to ypos local
                ypos -= yi[ynum].start + max_x_len;
                yseq = yi[ynum].seq;
                ylen = yi[ynum].len;

                //delimit / filter batch
                if (!batches_out)
                    scan_l2_g ();
                else
                    scan_l2_gb_1 ();
                k_matches++;
            }
        }
    }
}

//pass 2
//fast and crude batch delimiter / filter
void KT_SEARCH::scan_l2 ()
{
    char *sx, *sy;
    unsigned r, rr;
    int shx = (xpos & 3) << 1;
    int shy = (ypos & 3) << 1;
    int r_len, l_len, b_len, max_len, len;
    int score, max_score, matches;
    int m, msum = 0;

        //scan in forward direction (in 12-base steps)
    sx = xseq + (xpos >> 2);
    sy = yseq + (ypos >> 2);

    max_len = __min (xlen - xpos, ylen - ypos);
    len = 0, r_len = 0, rr = 0;
    score = 0, max_score = 0, matches = 0;
    while (score > max_score - stop_lev)
    {
        //compare 12 bases
        r = ~(GET32_U (sx) >> shx) ^ (GET32_U (sy) >> shy);
        r &= (r >> 1) & 0x555555;
        sx += 3, sy += 3, len += 12;

        //sequence boundary case
        if (len > max_len)
        {
            int rem = len - max_len;
            r <<= rem << 1;
            r &= 0xffffff;
            m = count_12 (r);
            msum += m, score += (m << 3) - ((12 - rem) << 2);
            if (score > max_score)
                max_score = score, matches = msum, rr = ~r, r_len = max_len;
            break;
        }

        m = count_12 (r);
        msum += m, score += (m << 3) - offs_lev;
        if (score > max_score)
        max_score = score, matches = msum, rr = ~r, r_len = len;
    }

    //trim edge mismatches
    while (r_len > 0 && (rr & 0x400000)) r_len --, rr <<= 2;

    //scan in reverse direction (in 12-base steps)
    sx = xseq + ((xpos - 12) >> 2);
    sy = yseq + ((ypos - 12) >> 2);

    //do not expand beyond previously found batch on the same diagonal
    max_len = __min (xpos - di[diag].b_end, ypos);
    //  max_len = __min (xpos, ypos);
    len = 0, l_len = 0, rr = 0;
    score = 0, max_score = 0, msum = matches;
    while (score > max_score - stop_lev)
    {
        //compare 12 bases
        r = ~(GET32_U (sx) >> shx) ^ (GET32_U (sy) >> shy);
        r &= (r >> 1) & 0x555555;
        sx -= 3, sy -= 3, len += 12;

        //sequence boundary case
        if (len > max_len)
        {
            int rem = len - max_len;
            r >>= rem << 1;
            m = count_12 (r);
            msum += m, score += (m << 3) - ((12 - rem) << 2);
            if (score > max_score)
                max_score = score, matches = msum, rr = ~r, l_len = max_len;
            break;
        }

        m = count_12 (r);
        msum += m, score += (m << 3) - offs_lev;
        if (score > max_score)
        max_score = score, matches = msum, rr = ~r & 0xffffff, l_len = len;
    }

    //trim edge mismatches
    while (l_len > 0 && (rr & 1)) l_len --, rr >>= 2;

    b_len = l_len + r_len;
    di[diag].score = 0;
    di[diag].xpos = xpos + r_len;

    //characterise homology level in % and compare to the threshold
    //penalise short batches
    if (((matches - penalty) * 100 > b_thresh_l * b_len) && (b_len >= l_thresh))
    {
        b_matches++;
        NN_SEQ sx_seq, sy_seq;

        sx_seq.uid = xuid;
        sx_seq.rev = xrev;
        sx_seq.len = xlen;
        sx_seq.seq = xseq;

        sy_seq.uid = yi[ynum].uid;
        sy_seq.rev = yi[ynum].rev;
        sy_seq.len = yi[ynum].len;
        sy_seq.seq = yseq;

        batch -> match_found (sx_seq, xpos - l_len, sy_seq, ypos - l_len, b_len, matches);
        di[diag].b_end = xpos + r_len;
    }
}

//pass 2 scoring routine
//allows short gaps
void KT_SEARCH::scan_l2_g ()
{
    unsigned r, rr, rx, ry;
    int x, y, xoffs, yoffs, offs, rem;
    int r_len, l_len, b_len, len;
    int score, max_score, matches, gaps;
    int s, m, max_am, msum = 0, gsum = 0;
    int cur_diag = diag;

    //scan in forward direction (in 12-base steps)
    x = xpos, y = ypos;
    len = 0, r_len = 0, rr = 0;
    rem = __max (x - xlen, y - ylen);
    score = 0, max_score = 0, matches = 0, gaps = 0;
    while (score > max_score - stop_lev)
    {
        //check current diag
        ry = get_12_bases (yseq, y);
        rx = get_12_bases (xseq, x);
        r = ~rx ^ ry, r &= (r >> 1) & 0x555555;
        len += 12, rem += 12;

        //sequence boundary case
        if (rem > 0)
        {
            r <<= rem << 1;
            r &= 0xffffff;
            m = count_12 (r);
            msum += m, score += (m << 3) - ((12 - rem) << 2);
            if (score > max_score)
                max_score = score, matches = msum, gaps = gsum, rr = ~r, r_len = len - rem;
            break;
        }

        //if local 12-tuple score below offs_lev check adjacent diags
        m = count_12 (r);
        s = (m << 3) - offs_lev;
        if (s < 0 && rem + max_offs < 0)
        {
            //merge 12-tuple matches on adjacent and current diags
            //self correlate both match patterns (r and ar)
            r &= (r << 2) | 1;
            for (offs = 1, max_am = 0; offs <= max_offs; offs++)
            {
                int am, ar;

                //go right from current diag
                ar = ~get_12_bases (xseq, x + offs) ^ ry;
                ar &= (ar >> 1) & 0x555555;
                ar &= (ar >> 2) | 0x400000;
                if ((am = count_12 (ar | r)) > max_am)
                max_am = am, xoffs = offs, yoffs = 0;

                //go up from current diag
                ar = ~get_12_bases (yseq, y + offs) ^ rx;
                ar &= (ar >> 1) & 0x555555;
                ar &= (ar >> 2) | 0x400000;
                if ((am = count_12 (ar | r)) > max_am)
                max_am = am, xoffs = 0, yoffs = offs;
            }

            //check if diagonal change makes sense
            if ((max_am << 3) - offs_lev >= 0)
            {
                //score += (max_am << 3) - offs_lev;
                msum += max_am;
                x += xoffs + 12, y += yoffs + 12, rem = __max (x - xlen, y - ylen);
                di[cur_diag].xpos = x;
                di[cur_diag].score = 0;
                cur_diag += yoffs - xoffs;
                gsum++;
                continue;
            }
        }

        score += s, msum += m;
        x += 12, y += 12;
        if (score > max_score)
        max_score = score, matches = msum, gaps = gsum, rr = ~r, r_len = len;
    }
    di[cur_diag].xpos = x;
    di[cur_diag].score = 0;

    //trim edge mismatches
    while (r_len > 0 && (rr & 0x400000)) r_len --, rr <<= 2;

    //scan in reverse direction (in 12-base steps)
    //do not expand beyond previously found batch on the same diagonal
    x = xpos - 12, y = ypos - 12;
    len = 0, l_len = 0, rr = 0;
    rem = __max (- x, - y);
    score = 0, max_score = 0, msum = matches, gsum = gaps;
    while (score > max_score - stop_lev)
    {
        //check current diag
        ry = get_12_bases (yseq, y);
        rx = get_12_bases (xseq, x);
        r = ~rx ^ ry, r &= (r >> 1) & 0x555555;
        len += 12;

        //sequence boundary case
        if (rem > 0)
        {
            r >>= rem << 1;
            m = count_12 (r);
            msum += m, score += (m << 3) - ((12 - rem) << 2);
            if (score > max_score)
                max_score = score, matches = msum, gaps = gsum, rr = ~r, l_len = len - rem;
            break;
        }

        //if local 12-tuple score below offs_lev check adjacent diags
        m = count_12 (r);
        s = (m << 3) - offs_lev;
        if (s < 0 && rem + max_offs < 0)
        {
            //merge 12-tuple matches on adjacent and current diags
            //self correlate both match patterns (r and ar)
            r &= (r >> 2) | 0x400000;
            for (offs = 1, max_am = 0; offs <= max_offs; offs++)
            {
                int am, ar;

                //go left from current diag
                ar = ~get_12_bases (xseq, x - offs) ^ ry;
                ar &= (ar >> 1) & 0x555555;
                ar &= (ar << 2) | 1;
                if ((am = count_12 (ar | r)) > max_am)
                max_am = am, xoffs = offs, yoffs = 0;

                //go down from current diag
                ar = ~get_12_bases (yseq, y - offs) ^ rx;
                ar &= (ar >> 1) & 0x555555;
                ar &= (ar << 2) | 1;
                if ((am = count_12 (ar | r)) > max_am)
                max_am = am, xoffs = 0, yoffs = offs;
            }

            //check if diagonal change makes sense
            if ((max_am << 3) - offs_lev >= 0)
            {
                //score += (max_am << 3) - offs_lev;
                msum += max_am;
                x -= xoffs + 12, y -= yoffs + 12, rem = __max (- x, - y);
                gsum++;
                continue;
            }
        }

        score += s, msum += m;
        x -= 12, y -= 12, rem += 12;
        if (score > max_score)
        max_score = score, matches = msum, gaps = gsum, rr = ~r & 0xffffff, l_len = len;
    }

    //trim edge mismatches
    while (l_len > 0 && (rr & 1)) l_len --, rr >>= 2;

    b_len = l_len + r_len;

    //characterise homology level in % and compare to the threshold
    //penalise short batches
    if (((matches - penalty) * 100 > b_thresh_l * b_len) && (b_len >= l_thresh) && b_len > gaps * g_period) //gaps can be replaced with (gaps + 1)
    {
        b_matches++;
        NN_SEQ sx_seq, sy_seq;

        sx_seq.uid = xuid;
        sx_seq.rev = xrev;
        sx_seq.len = xlen;
        sx_seq.seq = xseq;

        sy_seq.uid = yi[ynum].uid;
        sy_seq.rev = yi[ynum].rev;
        sy_seq.len = yi[ynum].len;
        sy_seq.seq = yseq;

        batch -> match_found (sx_seq, xpos - l_len, sy_seq, ypos - l_len, b_len, matches);
        di[diag].b_end = xpos + r_len;
    }
}


//pass 2 scoring routine
//allows short gaps
//fills bathces[] array
//merges 12-tuple matches around gap
void KT_SEARCH::scan_l2_gb_1 ()
{
  int i, j;
  unsigned r, rr, rx, ry;
  int x, y, lx, ly, xoffs, yoffs, offs, rem;
  int b_len, len;
  int score, max_score, matches;
  int s, m, max_am, msum;
  int num, l_num, r_num;
  int cur_diag = diag;
  BATCH batches[101];

    //scan in forward direction (in 12-base steps)
  r_num = num = 50;
  lx = x = xpos, ly = y = ypos;
  len = 0, b_len = 0, rr = 0;
  rem = __max (x - xlen, y - ylen);
  score = 0, max_score = 0, matches = msum = 0;
  while (score > max_score - stop_lev)
  {
    //check current diag
    ry = get_12_bases (yseq, y);
    rx = get_12_bases (xseq, x);
    r = ~rx ^ ry, r &= (r >> 1) & 0x555555;
    len += 12, rem += 12;

    //sequence boundary case
    if (rem > 0)
    {
      r <<= rem << 1;
      r &= 0xffffff;
      m = count_12 (r);
      msum += m, score += (m << 3) - ((12 - rem) << 2);
      if (score > max_score)
        matches = msum, r_num = num, rr = ~r, b_len = len - rem;
      break;
    }

    //if local 12-tuple score below offs_lev check adjacent diags
    m = count_12 (r);
    s = (m << 3) - offs_lev;
    if (s < 0 && rem + max_offs < 0 && num < 100)
    {
      //merge 12-tuple matches on adjacent and current diags
      //self correlate both match patterns (r and ar)
      r &= (r << 2) | 1;
      for (offs = 1, max_am = 0; offs <= max_offs; offs++)
      {
        int am, ar;

        //go right from current diag
        ar = ~get_12_bases (xseq, x + offs) ^ ry;
        ar &= (ar >> 1) & 0x555555;
        ar &= (ar >> 2) | 0x400000;
        if ((am = count_12 (ar | r)) > max_am)
          max_am = am, xoffs = offs, yoffs = 0;

        //go up from current diag
        ar = ~get_12_bases (yseq, y + offs) ^ rx;
        ar &= (ar >> 1) & 0x555555;
        ar &= (ar >> 2) | 0x400000;
        if ((am = count_12 (ar | r)) > max_am)
          max_am = am, xoffs = 0, yoffs = offs;
      }

      //check if diagonal change makes sense
      if ((max_am << 3) - offs_lev >= 0)
      {
//        score += (max_am << 3) - offs_lev;
        msum += max_am;
        batches[++num].xpos = x + xoffs;
        batches[num].ypos = y + yoffs;
        x += xoffs + 12, y += yoffs + 12, rem = __max (x - xlen, y - ylen);
        di[cur_diag].xpos = x;
        di[cur_diag].score = 0;
        cur_diag += yoffs - xoffs;
        continue;
      }
    }

    score += s, msum += m;
    if (score > max_score)
      max_score = score, matches = msum, r_num = num, rr = ~r, b_len = len;
    x += 12, y += 12;
  }
  di[cur_diag].xpos = x;
  di[cur_diag].score = 0;

  //trim edge mismatches
  while (b_len > 0 && (rr & 0x400000)) b_len --, rr <<= 2;

  //scan in reverse direction (in 12-base steps)
  l_num = num = 50;
  x = xpos - 12, y = ypos - 12;
  len = b_len, rr = 0;
  rem = __max (- x, - y);
  score = 0, max_score = 0, msum = matches;
  while (score > max_score - stop_lev)
  {
    //check current diag
    ry = get_12_bases (yseq, y);
    rx = get_12_bases (xseq, x);
    r = ~rx ^ ry, r &= (r >> 1) & 0x555555;
    len += 12;

    //sequence boundary case
    if (rem > 0)
    {
      r >>= rem << 1;
      m = count_12 (r);
      msum += m, score += (m << 3) - ((12 - rem) << 2);
      if (score > max_score)
        matches = msum, l_num = num, rr = ~r, b_len = len - rem, lx = x + rem, ly = y + rem;
      break;
    }

    //if local 12-tuple score below offs_lev check adjacent diags
    m = count_12 (r);
    s = (m << 3) - offs_lev;
    if (s < 0 && rem + max_offs < 0 && num > 0)
    {
      //merge 12-tuple matches on adjacent and current diags
      //self correlate both match patterns (r and ar)
      r &= (r >> 2) | 0x400000;
      for (offs = 1, max_am = 0; offs <= max_offs; offs++)
      {
        int am, ar;

        //go left from current diag
        ar = ~get_12_bases (xseq, x - offs) ^ ry;
        ar &= (ar >> 1) & 0x555555;
        ar &= (ar << 2) | 1;
        if ((am = count_12 (ar | r)) > max_am)
          max_am = am, xoffs = offs, yoffs = 0;

        //go down from current diag
        ar = ~get_12_bases (yseq, y - offs) ^ rx;
        ar &= (ar >> 1) & 0x555555;
        ar &= (ar << 2) | 1;
        if ((am = count_12 (ar | r)) > max_am)
          max_am = am, xoffs = 0, yoffs = offs;
      }

      //check if diagonal change makes sense
      if ((max_am << 3) - offs_lev >= 0)
      {
//        score += (max_am << 3) - offs_lev;
        msum += max_am;
        batches[num].xpos = x;
        batches[num--].ypos = y;
        x -= xoffs + 12, y -= yoffs + 12, rem = __max (- x, - y);
        continue;
      }
    }

    score += s, msum += m;
    if (score > max_score)
      max_score = score, matches = msum, l_num = num, rr = ~r & 0xffffff, b_len = len, lx = x, ly = y;
    x -= 12, y -= 12, rem += 12;
  }

  //trim edge mismatches
  while (b_len > 0 && (rr & 1)) b_len--, lx++, ly++, rr >>= 2;

  batches[l_num].xpos = lx;
  batches[l_num].ypos = ly;

  //characterise homology level in % and compare to the threshold
  //penalise short batches

  int approx_score = b_len;

  if (((matches - penalty) * 100 > b_thresh_l * b_len) && (b_len >= l_thresh) && b_len > (r_num - l_num) * g_period)
  {
    b_matches++;

    //now delimit batches with sub 12-tuple accuracy
    for (i = l_num; i < r_num; i++)
    {
      unsigned int r1, r2;
      int cnt, max_cnt, max_j;

      len = __min (batches[i + 1].xpos - batches[i].xpos, batches[i + 1].ypos - batches[i].ypos);

      r1 = ~get_12_bases (xseq, batches[i].xpos + len) ^ get_12_bases (yseq, batches[i].ypos + len);
      r1 &= (r1 >> 1) & 0x555555;

      r2 = ~get_12_bases (xseq, batches[i + 1].xpos) ^ get_12_bases (yseq, batches[i + 1].ypos);
      r2 &= (r2 >> 1) & 0x555555;

      //locate highest score gap position
      cnt = count_12 (r2);
      max_cnt = cnt, max_j = 0;
      for (j = 1; j <= 12; j++, r1 >>= 2, r2 >>= 2)
      {
        cnt -= r2 & 1;
        cnt += r1 & 1;
        if (cnt > max_cnt)
          max_cnt = cnt, max_j = j;
      }

      //adjust batch coordinates and lengths
      batches[i + 1].xpos += max_j;
      batches[i + 1].ypos += max_j;
      batches[i].len = len + max_j;
      b_len -= len + max_j;
    }
    batches[r_num].len = b_len;

    //process batches here ...
    //total number of batches = r_num - l_num + 1
    //first (leftmost) batch = batches[l_num];

    NN_SEQ sx_seq, sy_seq;

    sx_seq.uid = xuid;
    sx_seq.rev = xrev;
    sx_seq.len = xlen;
    sx_seq.seq = xseq;

    sy_seq.uid = yi[ynum].uid;
    sy_seq.rev = yi[ynum].rev;
    sy_seq.len = yi[ynum].len;
    sy_seq.seq = yseq;

    // batch_processor is allowed to modify batches, so adjust diag. position before calling it
    di[diag].b_end = batches [r_num].xpos + batches [r_num].len;
    batch_processor -> match_found (sx_seq, sy_seq, batches + l_num, r_num - l_num + 1, matches);
  }
}


//pass 2 scoring routine
//allows short gaps
//fills bathces[] array
//does not merges 12-tuple matches
void KT_SEARCH::scan_l2_gb_2 ()
{
  int i, j;
  unsigned r, rr, rx, ry;
  int x, y, lx, ly, xoffs, yoffs, offs, rem;
  int b_len, len;
  int score, max_score, matches;
  int s, m, max_am, msum;
  int num, l_num, r_num;
  int cur_diag = diag;
  BATCH batches [101];

    //scan in forward direction (in 12-base steps)
  r_num = num = 50;
  lx = x = xpos, ly = y = ypos;
  len = 0, b_len = 0, rr = 0;
  rem = __max (x - xlen, y - ylen);
  score = 0, max_score = 0, matches = msum = 0;
  while (score > max_score - stop_lev)
  {
    //check current diag
    ry = get_12_bases (yseq, y);
    rx = get_12_bases (xseq, x);
    r = ~rx ^ ry, r &= (r >> 1) & 0x555555;
    len += 12, rem += 12;

    //sequence boundary case
    if (rem > 0)
    {
      r <<= rem << 1;
      r &= 0xffffff;
      m = count_12 (r);
      msum += m, score += (m << 3) - ((12 - rem) << 2);
      if (score > max_score)
        matches = msum, r_num = num, rr = ~r, b_len = len - rem;
      break;
    }

    //if local 12-tuple score below offs_lev check adjacent diags
    m = count_12 (r);
    s = (m << 3) - offs_lev;
    if (s < 0 && rem + max_offs < 0 && num < 100)
    {
      for (offs = 1, max_am = 0; offs <= max_offs; offs++)
      {
        int am, ar;

        //go right from current diag
        ar = ~get_12_bases (xseq, x + offs) ^ ry;
        ar &= (ar >> 1) & 0x555555;
        if ((am = count_12 (ar)) > max_am)
          max_am = am, xoffs = offs, yoffs = 0;

        //go up from current diag
        ar = ~get_12_bases (yseq, y + offs) ^ rx;
        ar &= (ar >> 1) & 0x555555;
        if ((am = count_12 (ar)) > max_am)
          max_am = am, xoffs = 0, yoffs = offs;
      }

      //check if diagonal change makes sense
      if ((max_am << 3) - offs_lev >= 0)
      {
//        score += (max_am << 3) - offs_lev;
        msum += max_am;
        batches[++num].xpos = x + xoffs;
        batches[num].ypos = y + yoffs;
        x += xoffs + 12, y += yoffs + 12, rem = __max (x - xlen, y - ylen);
        di[cur_diag].xpos = x;
        di[cur_diag].score = 0;
        cur_diag += yoffs - xoffs;
        continue;
      }
    }

    score += s, msum += m;
    if (score > max_score)
      max_score = score, matches = msum, r_num = num, rr = ~r, b_len = len;
    x += 12, y += 12;
  }
  di[cur_diag].xpos = x;
  di[cur_diag].score = 0;

  //trim edge mismatches
  while (b_len > 0 && (rr & 0x400000)) b_len --, rr <<= 2;

  //scan in reverse direction (in 12-base steps)
  l_num = num = 50;
  x = xpos - 12, y = ypos - 12;
  len = b_len, rr = 0;
  rem = __max (- x, - y);
  score = 0, max_score = 0, msum = matches;
  while (score > max_score - stop_lev)
  {
    //check current diag
    ry = get_12_bases (yseq, y);
    rx = get_12_bases (xseq, x);
    r = ~rx ^ ry, r &= (r >> 1) & 0x555555;
    len += 12;

    //sequence boundary case
    if (rem > 0)
    {
      r >>= rem << 1;
      m = count_12 (r);
      msum += m, score += (m << 3) - ((12 - rem) << 2);
      if (score > max_score)
        matches = msum, l_num = num, rr = ~r, b_len = len - rem, lx = x + rem, ly = y + rem;
      break;
    }

    //if local 12-tuple score below offs_lev check adjacent diags
    m = count_12 (r);
    s = (m << 3) - offs_lev;
    if (s < 0 && rem + max_offs < 0 && num > 0)
    {
      for (offs = 1, max_am = 0; offs <= max_offs; offs++)
      {
        int am, ar;

        //go left from current diag
        ar = ~get_12_bases (xseq, x - offs) ^ ry;
        ar &= (ar >> 1) & 0x555555;
        if ((am = count_12 (ar)) > max_am)
          max_am = am, xoffs = offs, yoffs = 0;

        //go down from current diag
        ar = ~get_12_bases (yseq, y - offs) ^ rx;
        ar &= (ar >> 1) & 0x555555;
        if ((am = count_12 (ar)) > max_am)
          max_am = am, xoffs = 0, yoffs = offs;
      }

      //check if diagonal change makes sense
      if ((max_am << 3) - offs_lev >= 0)
      {
//        score += (max_am << 3) - offs_lev;
        msum += max_am;
        batches[num].xpos = x;
        batches[num--].ypos = y;
        x -= xoffs + 12, y -= yoffs + 12, rem = __max (- x, - y);
        continue;
      }
    }

    score += s, msum += m;
    if (score > max_score)
      max_score = score, matches = msum, l_num = num, rr = ~r & 0xffffff, b_len = len, lx = x, ly = y;
    x -= 12, y -= 12, rem += 12;
  }

  //trim edge mismatches
  while (b_len > 0 && (rr & 1)) b_len--, lx++, ly++, rr >>= 2;

  batches[l_num].xpos = lx;
  batches[l_num].ypos = ly;

  //characterise homology level in % and compare to the threshold
  //penalise short batches
  if (((matches - penalty) * 100 > b_thresh_l * b_len) && (b_len >= l_thresh) && b_len > (r_num - l_num) * g_period)
  {
    b_matches++;

    //now delimit batches with sub 12-tuple accuracy
    for (i = l_num; i < r_num; i++)
    {
      unsigned int r1, r2;
      int cnt, max_cnt, max_j;

      len = __min (batches[i + 1].xpos - batches[i].xpos, batches[i + 1].ypos - batches[i].ypos);

      r1 = ~get_12_bases (xseq, batches[i].xpos + len) ^ get_12_bases (yseq, batches[i].ypos + len);
      r1 &= (r1 >> 1) & 0x555555;

      r2 = ~get_12_bases (xseq, batches[i + 1].xpos) ^ get_12_bases (yseq, batches[i + 1].ypos);
      r2 &= (r2 >> 1) & 0x555555;

      //locate highest score gap position
      cnt = count_12 (r2);
      max_cnt = cnt, max_j = 0;
      for (j = 1; j <= 12; j++, r1 >>= 2, r2 >>= 2)
      {
        cnt -= r2 & 1;
        cnt += r1 & 1;
        if (cnt > max_cnt)
          max_cnt = cnt, max_j = j;
      }

      //adjust batch coordinates and lengths
      batches[i + 1].xpos += max_j;
      batches[i + 1].ypos += max_j;
      batches[i].len = len + max_j;
      b_len -= len + max_j;
    }
    batches[r_num].len = b_len;

    //process batches here ...
    //total number of batches = r_num - l_num + 1
    //first (leftmost) batch = batches[l_num];

    NN_SEQ sx_seq, sy_seq;

    sx_seq.uid = xuid;
    sx_seq.rev = xrev;
    sx_seq.len = xlen;
    sx_seq.seq = xseq;

    sy_seq.uid = yi[ynum].uid;
    sy_seq.rev = yi[ynum].rev;
    sy_seq.len = yi[ynum].len;
    sy_seq.seq = yseq;

    // batch_processor is allowed to modify batches, so adjust diag. position before calling it
    di[diag].b_end = batches [r_num].xpos + batches [r_num].len;
    batch_processor -> match_found (sx_seq, sy_seq, batches + l_num, r_num - l_num + 1, matches);
  }
}
