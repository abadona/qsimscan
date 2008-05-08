
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

AlignResult::AlignResult (longlong uid, bool reverse, float al_score, int score, float chi2, double evalue, double bitscore, int q_auto_score, int t_auto_score, int batch_no, BATCH* batches)
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
batch_no_ (batch_no)
{
    if (batch_no_)
        if (!(batches_ = new BATCH [batch_no_]))
            Error (MemoryRerror);
    std::copy (batches, batches + batch_no, (BATCH*) batches_);
}

AlignResult::AlignResult ()
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
batch_no_ (0)
{
}

bool AlignResult::operator < (const AlignResult& other) const
{
    if (al_score_ < other.al_score_) return true;
    if (al_score_ > other.al_score_) return false;
    if (evalue_ < other.evalue_) return true;
    if (evalue_ > other.evalue_) return false;
    if (chi2_ < other.chi2_) return true;
    return false;
}
bool AlignResult::operator == (const AlignResult& other) const
{
    return (al_score_ == other.al_score_ &&
            evalue_ == other.evalue_ &&
            chi2_ == other.chi2_);

}



