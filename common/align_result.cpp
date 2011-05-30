
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

#include "align_result.h"
#include "sequtil.h"
#include <rerror.h>
#include <algorithm>
#include <cstring>

AlignResult::AlignResult (longlong uid, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, const BATCH* batches, const char* binsubj, QWORD subjid, DWORD subjlen)
:
uid_ (uid),
reverse_ (reverse),
al_score_ (al_score),
score_ (score),
chi2_ (chi2),
evalue_ (evalue),
bitscore_ (bitscore),
q_auto_score_ (q_auto_score),
t_auto_score_ (t_auto_score),
batch_no_ (batch_no),
subjid_ (subjid),
subjlen_ (subjlen)
{
    if (batch_no_)
    {
        batches_ = new BATCH [batch_no_];
        std::copy (batches, batches + batch_no, (BATCH*) batches_);
        if (binsubj)
        {
            unsigned subj_len = (batches + batch_no - 1)->ypos + (batches + batch_no - 1)->len - batches->ypos;
            subject_ = new char [subj_len + 1];
            n_binary2ascii (subject_, subj_len + 1, binsubj, batches->ypos, subj_len);
        }
    }
}


bool AlignResult::intersects (const AlignResult& other) const
{
    for (int b1 = 0; b1 != batch_no_; b1 ++)
        for (int b2 = 0; b2 != other.batch_no_; b2++)
            if (batches_ [b1].intersects (other.batches_ [b2]))
                return true;
    return false;
}


