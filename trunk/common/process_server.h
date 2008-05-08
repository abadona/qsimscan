
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

#ifndef __process_server_h__
#define __process_server_h__

#pragma warning (disable: 4786)

// in order to implement an AppServer, the following should be defined:
// - Process* process_factory ();
// - Process_params* params_factory ();
// - define classes inherited from Process (process_thread.h) and Process_params (process_params.h)
// Also should link to to console_util_main.obj

#include "process_thread.h"
#include "process_params.h"

Process* process_factory ();
Process_params* process_params_factory ();
const char* process_name ();


#endif
