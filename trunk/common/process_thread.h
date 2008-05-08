
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

#ifndef __process_thread_h__
#define __process_thread_h__

#include <time.h>
#include <ostream>
#include <sstream>
#include <string>
#include "platform.h"
#include "sync_control.h"


class Process_params;

class Process
{
public:
    typedef enum
    {
        APS_NOT_INIT,
        APS_READY,
        APS_RUNNING,
        APS_INTERRUPTED,
        APS_ERROR,
        APS_COMPLETE
    }
    STATE;

//private:

//    unsigned long processing_thread_;

protected:

    // CRITICAL_SECTION critical_section_;
#ifndef NOTHREADS
    SyncControl sync;
#endif

    STATE state_;
    bool interrupt_;     // shared
    bool save_results_;

    bool error_;
    std::string error_string_;

    time_t start_time_sec_;
    time_t end_time_sec_;
    time_t start_time_;
    time_t subphase_start_time_;
    time_t end_time_;
    time_t last_time_;

    longlong subphase_start_current_;
    longlong last_current_;

    longlong total_;     // shared
    longlong current_;   // shared
    longlong skipped_;   // shared
    longlong res_no_;    // shared
    long phase_;         // shared
    long subphase_;      // shared

    time_t report_interval_;


    // handles
    virtual bool init_handler (Process_params* params) = 0;
    virtual bool next_handler () = 0;
    virtual bool report_results_handler () = 0;
    virtual bool close_handler () = 0;


    // process stages
    virtual bool start ();
    virtual bool init (Process_params* params);
    virtual bool next ();
    virtual bool report_results ();
    virtual bool close ();

    // status control
    void total   (longlong val);
    void current (longlong val);
    void skipped (longlong val);
    void res_no  (longlong val);
    void phase   (long phase_no);
    void subphase (long subphase_no);
    longlong incrcur ();
    longlong incrskipped ();
    longlong incrresno ();
    long incrphase ();
    long incrsubphase ();
    void error   (const char* errstr = NULL);
    void error   (std::string& errstr);
    void error   (std::ostringstream& errstr);

public:

    Process ();
    virtual ~Process ();

    virtual bool run (Process_params* params);
    virtual bool stop (bool save_results);
    virtual bool wait (unsigned msec);

    virtual STATE state ();
    virtual const char* getError ();

    // naming convention
    virtual const char* process_name () = 0;
    virtual const char* description () = 0;
    virtual const char* version () = 0;
    virtual const char* processing_item_name () = 0;
    virtual const char* result_item_name () = 0;

    // process description
    virtual const char* phase_name (int phase_no) = 0;
    virtual int phases_total () = 0;
    virtual const char* subphase_name (int phase_no, int subphase_no) = 0;
    virtual int subphases_total (int phase_no) = 0;


    // state and speed
    int         stage       ();
    longlong    total       ();
    longlong    current     ();
    longlong    skipped     ();
    longlong    res_no      ();
    long        phase       ();
    long        subphase    ();

    virtual void report_state (std::ostream& o);
    virtual void report_self (std::ostream& o);
    virtual void report_header (std::ostream& o);
    virtual void report_phase_switch (std::ostream& o);
    virtual void report_final (std::ostream& o);
    double percent_done () ;
    double average_speed () ;
    double last_speed () ;
    int    report_interval (); // in ticks; 1 tick = 1/CLOCK_PER_SEC seconds

    time_t start_time_sec () const {return start_time_sec_;}
    time_t subphase_start_time_sec () const {return subphase_start_time_;}
    time_t end_time_sec () const {return end_time_sec_;}
    time_t cur_time_sec () const;

    friend void process_thread_function (void*);
};

#endif // __process_thread_h__
