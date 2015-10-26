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
// For any questions please contact SciDM team by email at team@scidm.org
//
// This module is modified for use in PDOP MexSat project by Hughes Network Systems, LLC
//
//////////////////////////////////////////////////////////////////////////////
#ifndef __tracer_h__
#define __tracer_h__

#include <ostream>
#include <iomanip>
#include <ctime>
#include <pthread.h>
#include "common_str.h"


class Trace__
{
public:
    const char* fname_;
    int lno_;
    const char* func_;
    time_t time_;
    pthread_t thread_;

    Trace__ (const char* fname = NULL, int lno = 0, const char* funcname = NULL)
    :
    fname_ (fname),
    lno_ (lno),
    func_ (funcname),
    time_ (time (NULL)),
    thread_ (pthread_self ())
    {
    }
};


std::ostream& operator << (std::ostream& e, Trace__ t);

#define Trace Trace__(__FILE__, __LINE__, __FUNCTION__)

class Logger
{
    bool logger_on_;
public:

    enum LEVEL
    {
        CRITICAL = 0,
        ERROR   = 10,
        WARNING = 20,
        INFO    = 30,
        DEBUG   = 40,
        TRACE   = 50,
    };

    std::ostream& o_;

    Logger (bool enabled = false);
    Logger (std::ostream& o, bool enabled = true);

    void enable  (bool op = true) { logger_on_ = op;  }
    void disable () { logger_on_ = false; }
    bool enabled () const { return logger_on_; }
};

const Logger::LEVEL loglevels [] = {Logger::CRITICAL, Logger::ERROR, Logger::WARNING, Logger::INFO, Logger::DEBUG, Logger::TRACE};


// output operator support
template <class TT>
Logger& operator << (Logger& logger, const TT& operand)
{
    if (logger.enabled ())
        logger.o_ << operand;
    return logger;
}
// ostream manipulators support
inline Logger& operator << (Logger& logger, std::ostream& (*op) (std::ostream&))
{
    if (logger.enabled ())
        logger.o_ << op;
    return logger;
}

#define clreol "\x1B[K"

#ifndef __tracer_cpp__
extern Logger trclog;   // common tracing, off by default
extern Logger dbglog;
extern Logger info;
extern Logger warnlog;
extern Logger errlog;
#endif

#define trc (trclog << Trace)
#define dbg (dbglog << Trace << "Debug: ")
#define errl (errlog << Trace << "Error: ")
#define warn (warnlog << "Warning: ")

// convinience function setting logging level
void set_logging_level (Logger::LEVEL level);

#define IND(x) std::setw (x) << EMPTY_STR << std::setw (0)


#endif
