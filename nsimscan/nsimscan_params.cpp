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

#include "nsimscan_params.h"

static const char* HEADER = "Tool for searching for nucleotide sequence similarities (based on YABLAST algorithm by SciDM)";
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
