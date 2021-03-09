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

#ifndef __process_server_h__
#define __process_server_h__

#ifdef _MSC_VER
#pragma warning (disable: 4786)
#endif

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
