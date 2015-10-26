//////////////////////////////////////////////////////////////////////////////
// This software module is developed by SCIDM team in 1999-2015.
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
// For any questions please contact SciDM team by email at scidmteam@yahoo.com
//
// This module is modified for use in PDOP MexSat project by Hughes Network Systems, LLC
//
//////////////////////////////////////////////////////////////////////////////
#include "tracer.h"
#include <ctime>

#include <iostream>

static const int TMBUF_SZ = 64;
std::ostream& operator << (std::ostream& e, Trace__ t)
{
    struct tm* timeinfo;
    char tmbuffer [TMBUF_SZ];

    timeinfo = localtime (&(t.time_));
    strftime (tmbuffer, TMBUF_SZ, "[%x %X %Z]", timeinfo);

    if (t.func_)  e << t.func_ << ": ";
    if (t.fname_) e << t.fname_ << ":" << t.lno_ << ": ";
    if (t.func_) e << "thr " << t.thread_ << ":";
    e << tmbuffer << " ";
    return e;
}

Logger::Logger (bool enabled)
:
logger_on_ (enabled),
o_ (std::cerr)
{
}

Logger::Logger (std::ostream& o, bool enabled)
:
logger_on_ (enabled),
o_ (o)
{
}


Logger trclog (false);
Logger dbglog (false);
Logger info (false);
Logger warnlog (true);
Logger errlog (true);

void set_logging_level (Logger::LEVEL level)
{
    trclog.enable   (level >= Logger::TRACE);
    dbglog.enable   (level >= Logger::DEBUG);
    info.enable  (level >= Logger::INFO);
    warnlog.enable  (level >= Logger::WARNING);
    errlog.enable   (level >= Logger::ERROR);
}

