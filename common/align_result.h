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

#ifndef __align_result_h__
#define __align_result_h__
#include <platform.h>
#include <resource.h>
#include <sciminmax.h>
#include "align_batch.h"
#include <vector>

struct AlignResult
{
    longlong    uid_;
    bool        reverse_;
    float       al_score_;
    int         score_;
    float       chi2_;
    double      evalue_;
    double      bitscore_;
    int         q_auto_score_;
    int         t_auto_score_;
    int         batch_no_;
    MemWrapper <BATCH> batches_;
    MemWrapper <char> subject_;
    QWORD       subjid_;
    DWORD       subjlen_;

    AlignResult ();
    AlignResult (longlong uid, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, const BATCH* batches, const char* binsubj = NULL, QWORD subjid = 0, DWORD subjlen = 0);
    bool operator < (const AlignResult& other) const;
    bool operator == (const AlignResult& other) const;
    bool overlaps (const AlignResult& other, unsigned max_orphan = 0) const;
    bool intersects (const AlignResult& other) const;
    bool can_continue (const AlignResult& other, unsigned max_ovl, bool& runoff) const;
    unsigned gap (const AlignResult& next) const;
};

typedef std::vector <AlignResult> ARVect;

inline AlignResult::AlignResult ()
:
uid_ (0),
reverse_ (false),
al_score_ (0.0),
score_ (0),
chi2_ (0.0),
evalue_ (0.0),
bitscore_ (0.0),
q_auto_score_ (0),
t_auto_score_ (0),
batch_no_ (0),
subjid_ (0),
subjlen_ (0)
{
}

inline bool AlignResult::operator < (const AlignResult& other) const
{
    if (al_score_ < other.al_score_) return true;
    if (al_score_ > other.al_score_) return false;
    if (evalue_ < other.evalue_) return true;
    if (evalue_ > other.evalue_) return false;
    if (chi2_ < other.chi2_) return true;
    return false;
}

inline bool AlignResult::operator == (const AlignResult& other) const
{
    return (al_score_ == other.al_score_ &&
            evalue_ == other.evalue_ &&
            chi2_ == other.chi2_);

}

inline bool AlignResult::overlaps (const AlignResult& other, unsigned max_orphan) const
{
    const BATCH& fb1 = batches_ [0];
    const BATCH& lb1 = batches_ [batch_no_ - 1];
    const BATCH& fb2 = other.batches_ [0];
    const BATCH& lb2 = other.batches_ [other.batch_no_ - 1];
    int bx1 = fb1.xpos;
    int ex1 = lb1.xpos + lb1.len;
    int by1 = fb1.ypos;
    int ey1 = lb1.ypos + lb1.len;
    int bx2 = fb2.xpos;
    int ex2 = lb2.xpos + lb2.len;
    int by2 = fb2.ypos;
    int ey2 = lb2.ypos + lb2.len;
    if ((bx1 > ex2 || bx2 > ex1) &&
        (by1 > ey2 || by2 > ey1))
        return false; // no overlap on either axes
    // shorter similarity should not have unique zone longer then max_orphan
    // process x
    // check uniq zone on shorter sim
    int uzlen;
    if (bx1 < ex2 && bx2 < ex1) // x overlaps
    {
        if (ex1 - bx1 < ex2 - bx2) // 1 shorter
            uzlen = max_ (0, max_ (ex1 - ex2, bx2 - bx1));
        else
            uzlen = max_ (0, max_ (ex2 - ex1, bx1 - bx2));
        if (uzlen > max_orphan)
            return false;
    }
    // process y
    // check uniq zone on shorter sim
    if (by1 < ey2 && by2 < ey1) // x overlaps
    {
        if (ey1 - by1 < ey2 - by2) // 1 shorter
            uzlen = max_ (0, max_ (ey1 - ey2, by2 - by1));
        else
            uzlen = max_ (0, max_ (ey2 - ey1, by1 - by2));
        if (uzlen > max_orphan)
            return false;
    }
    return true;
}

inline bool AlignResult::can_continue (const AlignResult& other, unsigned max_ovl, bool& runoff) const
{
    const BATCH& fb1 = other.batches_ [0];
    const BATCH& lb1 = other.batches_ [other.batch_no_ - 1];
    const BATCH& fb2 = batches_ [0];

    // check if 1 can preceed 2
    if (lb1.xpos > fb2.xpos + max_ovl)
    {
        runoff = true;
        return false;
    }
    else runoff = false;
    if (fb1.xpos >= fb2.xpos) // 1 starts after 2 - cannot preceed
        return false;
    if (lb1.xpos + lb1.len > fb2.xpos + max_ovl)
        return false; // prev ends too far by X (x overlap too long)
    if (fb1.ypos >= fb2.ypos) // 1 starts after 2 - cannot preceed
        return false;
    if (lb1.ypos + lb1.len > fb2.ypos + max_ovl)
        return false; // prev ends too far by Y (y overlap too long)
    return true;
}

inline unsigned AlignResult::gap (const AlignResult& next) const
{
    const BATCH& b1 = batches_ [batch_no_ - 1];
    const BATCH& b2 = next.batches_ [0];
    unsigned xgap = max_ (0, int (b2.xpos) - int (b1.xpos + b1.len));
    unsigned ygap = max_ (0, int (b2.ypos) - int (b1.ypos + b1.len));
    return xgap + ygap;
}

#endif
