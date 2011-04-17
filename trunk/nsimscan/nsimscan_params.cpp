
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

#include "nsimscan_params.h"

static const char* HEADER = "Tool for searching for nucleotide sequence similarities (based on QSIMSCAN algorithm by SciDM)";
extern const char* VERSION;

Nsimscan_params::Nsimscan_params ()
:
KTSearch_params (HEADER, NULL, VERSION)
{
}

Process_params* process_params_factory ()
{
    return new Nsimscan_params ();
}
