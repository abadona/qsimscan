
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
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// For any questions please contact SciDM team by email at scidmteam@yahoo.com
//////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <exception>
#include "process_thread.h"
#include "i64out.h"
#include "rerror.h"

// #define NOTHREADS

//#ifdef NOTHREADS
#include <iostream>
// #endif


// define NOTHREADS to switch off multithreading
// #define NOTHREADS
// define DONT_CATCH_UNKNOWN to pass unknown exceptions to OS handlers

#ifdef _DEBUG
#define DONT_CATCH_UNKNOWN 1
#endif

static const int default_report_interval = 2; // 2 seconds!

#ifndef NOTHREADS
#define CRBEG_P process->sync.lock ();
#define CREND_P process->sync.unlock ();
#define CRBEG sync.lock ();
#define CREND sync.unlock ();;
#else
#define CRBEG_P ;
#define CREND_P ;
#define CRBEG ;
#define CREND ;
#endif


void process_thread_function (void* arglist)
{
    Process* process = (Process*) arglist;
    bool finish = false;
    time_t prev_check = time (NULL);
    time_t report_interval = process->report_interval_;

    // record start time
    CRBEG_P
    process->start_time_ = time (NULL);
    CREND_P

    do
    {
        time_t cur_check = time (NULL);
        if (cur_check - prev_check > report_interval)
        {
            prev_check = cur_check;
            CRBEG_P
            finish = process->interrupt_;
            CREND_P
        }
        if (finish) break;

        // process next
        if (!process->next ())
        {
            CRBEG_P
            finish = true;
            CREND_P
        }
        if (process->error_) // no need for synchronisation here
            process->save_results_ = false;
    }
    while (!finish);

    // record end time
    CRBEG_P
    process->end_time_ = time (NULL);
    CREND_P

    // report results and close
    if (process->save_results_)
        process->report_results ();

    process->close ();
}

Process::Process ()
:
last_current_ (0),
report_interval_ (default_report_interval),
interrupt_ (false),
save_results_ (true),
error_ (false),
state_ (APS_NOT_INIT),
total_ (0),
current_ (0),
subphase_start_current_ (0),
skipped_ (0),
res_no_ (0),
phase_ (0),
subphase_ (0)
{
}

Process::~Process ()
{
}

bool Process::run (Process_params* params)
{
    if (!init (params)) return false;
    if (!start ()) return false;
    return true;
}

bool Process::wait (unsigned msec)
{
    if (state () != APS_RUNNING) return true;
    // WaitForSingleObject ((HANDLE) processing_thread_, msec);
#ifndef NOTHREADS
    sync.wait (msec);
    if (state () != APS_RUNNING) return true;
#endif
    return false;
}

bool Process::start ()
{
    if (state () != APS_READY) return false;
    CRBEG
    state_ = APS_RUNNING;
    CREND

    start_time_sec_ = time (NULL);
    subphase_start_time_ = start_time_ = last_time_ = time (NULL);
    subphase_start_current_ = current_;

#ifndef NOTHREADS
    sync.begin (process_thread_function, this);
#else
    process_thread_function (this);
#endif
    return true;
}

bool Process::stop (bool save_results)
{
    CRBEG
    interrupt_ = true;
    save_results_ = save_results;
    CREND
    return true;
}

bool Process::init (Process_params* params)
{
    STATE curstate = state ();
    if (curstate == APS_RUNNING || curstate == APS_READY) return false;
    CRBEG
    total_ = 0;
    current_ = 0;
    subphase_start_current_ = 0;
    skipped_ = 0;
    res_no_ = 0;
    phase_ = 0;
    subphase_ = 0;
    interrupt_ = false;
    error_ = false;
    error_string_.erase ();
    CREND
    bool toRet = false;
    try
    {
        toRet = init_handler (params);
    }
    catch (Rerror& e)
    {
        std::ostringstream str;
        str << "Error during INITIALIZATION stage: " << (const char*) e << std::endl << std::flush;
        error (str);
    }
#if !defined (DONT_CATCH_UNKNOWN)
    catch (std::exception& e)
    {
        std::ostringstream str;
        str << "Standard exception: " << e.what () << " caught during INITIALIZATION stage " << std::endl << std::flush;
        error (str);
    }
    catch (...)
    {
        std::ostringstream str;
        str << "Unhandled exception caught during INITIALIZATION stage " << std::endl << std::flush;
        error (str);
    }
#endif
    if (toRet)
    {
        CRBEG
        save_results_ = true;
        state_ = APS_READY;
        CREND
    }
    else
        state_ = APS_ERROR;
    return toRet;
}

bool Process::next ()
{
    // process next
    bool toRet = false;
    try
    {
        toRet = next_handler ();
    }
    catch (Rerror& e)
    {
        std::ostringstream str;
        str << "Error during PROCESSING stage: " << (const char*) e << std::endl
            << "phase " << phase () << " (" << phase_name (phase ())
            << ") , subphase " << subphase () << " (" << subphase_name (phase (), subphase ())
            << ") at " << processing_item_name () << " " << current () << " of " << total () << "." << std::endl << std::flush;
        error (str);
    }
#if !defined (DONT_CATCH_UNKNOWN)
    catch (std::exception& e)
    {
        std::ostringstream str;
        str << "Standard exception: " << e.what () << " caught during PROCESSING stage " << std::endl
            << "phase " << phase () << " (" << phase_name (phase ())
            << ") , subphase " << subphase () << " (" << subphase_name (phase (), subphase ())
            << ") at " << processing_item_name () << " " << current () << " of " << total () << "." << std::endl << std::flush;
        error (str);
    }
    catch (...)
    {
        std::ostringstream str;
        str << "Unhandled exception caught during PROCESSING stage " << std::endl
            << "phase " << phase () << " (" << phase_name (phase ())
            << ") , subphase " << subphase () << " (" << subphase_name (phase (), subphase ())
            << ") at " << processing_item_name () << " " << current () << " of " << total () << "." << std::endl << std::flush;
        error (str);
    }
#endif
    // toRet could be false in the case of error or in the case of finished processing.
    // They are distinguished by the state of _error member of ProcessThread
    return toRet;
}

bool Process::report_results ()
{
    end_time_sec_ = time (NULL);
    bool toRet = false;

    try
    {
        toRet = report_results_handler ();
    }
    catch (Rerror& e)
    {
        std::ostringstream str;
        str << "Error during REPORT stage: " << (const char*) e << std::endl << std::flush;
        error (str);
    }
#if !defined (DONT_CATCH_UNKNOWN)
    catch (std::exception& e)
    {
        std::ostringstream str;
        str << "Exception: " << e.what () << " caught during results reproting" << std::endl << std::flush;
        error (str);
    }
    catch (...)
    {
        std::ostringstream str;
        str << "Unhandled exception caught during reporting results" << std::endl << std::flush;
        error (str);
    }
#endif
    return toRet;
}

bool Process::close ()
{
    STATE curstate = state ();
    if (curstate != APS_RUNNING)
    {
        error ("State error: closing running process");
        return false;
    }
    bool toRet = false;
    try
    {
        toRet = close_handler ();
    }
    catch (Rerror& e)
    {
        std::ostringstream str;
        str << "Error during CLOSE stage: " << (const char*) e << std::endl << std::flush;
        error (str);
    }
#if !defined (DONT_CATCH_UNKNOWN)
    catch (std::exception& e)
    {
        error (e.what ());
    }
    catch (...)
    {
        std::ostringstream str;
        str << "Unhandled exception caught during CLOSE stage " << std::endl << std::flush;
        error (str);
    }
#endif

    CRBEG
    if (!toRet)
    {
        state_ = APS_ERROR;
    }
    else
    {
        if (error_) state_ = APS_ERROR;
        else if (interrupt_) state_ = APS_INTERRUPTED;
        else state_ = APS_COMPLETE;
    }
    CREND
    return toRet;
}

longlong Process::total   ()
{
    CRBEG
    longlong toRet = total_;
    CREND
    return toRet;
}

longlong Process::current ()
{
    CRBEG
    longlong toRet = current_;
    CREND
    return toRet;
}

longlong Process::skipped   ()
{
    CRBEG
    longlong toRet = skipped_;
    CREND
    return toRet;
}

longlong Process::res_no  ()
{
    CRBEG
    longlong toRet = res_no_;
    CREND
    return toRet;
}

long Process::phase ()
{
    CRBEG
    long toRet = phase_;
    CREND
    return toRet;
}

long Process::subphase ()
{
    CRBEG
    long toRet = subphase_;
    CREND
    return toRet;
}

void Process::total   (longlong val)
{
    CRBEG
    total_ = val;
    CREND
}

void Process::current (longlong val)
{
    CRBEG
    current_ = val;
    CREND
}

void Process::skipped   (longlong val)
{
    CRBEG
    skipped_ = val;
    CREND
}

void Process::res_no  (longlong val)
{
    CRBEG
    res_no_ = val;
    CREND
}

void Process::phase (long phase_no)
{
    CRBEG
    phase_ = phase_no;
    subphase_ = 0;
    subphase_start_time_ = time (NULL);
    subphase_start_current_ = current_;
    CREND
}

void Process::subphase (long subphase_no)
{
    CRBEG
    subphase_ = subphase_no;
    subphase_start_time_ = time (NULL);
    subphase_start_current_ = current_;
    CREND
}

longlong Process::incrcur (longlong incr)
{
    CRBEG
    longlong toR = current_ += incr;
    CREND
    return toR;
}

longlong  Process::incrskipped (longlong incr)
{
    CRBEG
    longlong toR = skipped_ += incr;
    CREND
    return toR;
}

longlong Process::incrresno (longlong incr)
{
    CRBEG
    longlong toR = res_no_ += incr;
    CREND
    return toR;
}

long Process::incrphase ()
{
    bool err = false;
    CRBEG
    long toR = phase_ ++;
    if (phase_ >= phases_total ()) err = true;
    subphase_= 0;
    CREND
    if (err) ERR("Phase number overflow");
    return toR;
}

long Process::incrsubphase ()
{
    bool err = false;
    CRBEG
    long toR = subphase_ ++;
    if (subphase_ >= subphases_total (phase_)) err = true;
    CREND
    if (err) ERR("Subphase number overflow");
    return toR;
}

void Process::report_state (std::ostream& o)
{
    o << "\r" << phase_name (phase ()) << " : " << subphase_name (phase (), subphase ()) << ": "
      << processing_item_name () << " " << current () << " of " << total ()
      << " (" << std::setprecision (2) << std::setiosflags (std::ios::fixed) << percent_done () << "%), speed is "
      << average_speed () << " " << processing_item_name () << "/sec; "
      << res_no () << " " << result_item_name () << " produced.    "
      << std::flush;
}

void Process::report_self (std::ostream& o)
{
    //report_header (o);
}

void Process::report_phase_switch (std::ostream& o)
{
    o << std::endl;
}

void Process::report_header (std::ostream& o, const char* procname)
{
    o << (procname ? procname : process_name ()) << " version " << version () << ", by Scientific Data Management (SciDM.org)" << std::endl << "  " << description () << std::endl << std::flush;
}

void Process::report_final (std::ostream& o)
{
    o << std::endl << process_name () << " completed succesfully." << std::endl;
}

double Process::percent_done ()
{
    longlong tot = total ();
    longlong cur = current ();
    if (tot < 1)
        return 100.0;
    if (cur < 0)
        cur = 0;
    return double (cur) * 100.0 / tot;
}

double Process::average_speed ()
{
    time_t cur_time = time (NULL);
    int t_diff = cur_time - subphase_start_time_;
    if (!t_diff) t_diff = 1;
    longlong cur_current = current () - subphase_start_current_;
    double speed = double (cur_current) / double (t_diff);
    return speed;
}

double Process::last_speed ()
{
    time_t cur_time = time (NULL);
    int t_diff = cur_time - last_time_;
    if (!t_diff) t_diff = 1;
    longlong cur_current = current ();
    double speed = double (cur_current - last_current_) / double (t_diff);
    last_current_ = cur_current;
    last_time_ = cur_time;
    return speed;
}

int Process::report_interval ()
{
    CRBEG
    int toRet = report_interval_;
    CREND
    return toRet;
}

void Process::error (const char* errstr)
{
    if (error_string_.length ())
      error_string_.append ("\n");
    error_string_.append (errstr);
    error_ = true;
}

void Process::error (std::string& errstr)
{
    if (error_string_.length ())
      error_string_.append ("\n");
    error_string_.append (errstr);
    error_ = true;
}
void Process::error (std::ostringstream& errstr)
{
    if (error_string_.length ())
      error_string_.append ("\n");
    error_string_.append (errstr.str ());
    error_ = true;
}

Process::STATE Process::state ()
{
    CRBEG
    STATE toRet = state_;
    CREND
    return toRet;
}

const char* Process::getError ()
{
    return error_string_.c_str ();
}

time_t Process::cur_time_sec () const
{
    return time (NULL);
}
