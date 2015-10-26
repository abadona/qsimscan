
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

#include "process_thread.h"
#include "process_params.h"
#include "resource.h"
#include "rerror.h"
#include "tracer.h"
#include <signal.h>
#include <iostream>
#include <time.h>
#include "process_server.h"

static Process* proc = NULL;

static void terminate (int sig)
{
    if (proc)
        proc->stop (true);
}

int console_main (int argc, char* argv [])
{
    ObjWrapper <Process_params> params = process_params_factory ();

    if (!params->parseCmdline (argc, argv))
    {
        if (!params->help_mode ())
        {
            params->cmdline ()->reportErrors (std::cerr);
            return -1;
        }
        else
            return 0;
    }
    if (!params->process ())
        ers << "Error pocessing parameters" << Throw;

    // set the logging level according to verbose / debug parameters
    Logger::LEVEL level = params->verbose () ? Logger::INFO : Logger::WARNING;
    if (params->debug () < 0) 
        params->debug (0);
    if (params->debug () >= (sizeof (loglevels) / sizeof (*loglevels))) 
        params->debug ((sizeof (loglevels) / sizeof (*loglevels)) - 1);
    if (level < loglevels [params->debug ()])
        level = loglevels [params->debug ()];
    set_logging_level (level);

    ObjWrapper <Process> process = process_factory ();
    if (!process)
        return -1;
    proc = process;

    // set cbrk handler
    signal (SIGINT, terminate);

    if (params->verbose ())
        process->report_header (std::cerr, params->procname ());
    if (process->run (params))
    {
        int timeslice = process->report_interval ();
        int prev_phase = process->phase ();
        int prev_subphase = process->subphase ();
        if (params->verbose ())
            process->report_self (std::cerr);
        time_t prevtime = 0;
        while (!process->wait (timeslice * 1000)) // timeslice is expressed in seconds
        {
            if (params->verbose ())
            {
                time_t curtime = time (NULL);
                if (curtime - prevtime >= timeslice)
                {
                    int cur_phase = process->phase ();
                    int cur_subphase = process->subphase ();
                    if (cur_phase != prev_phase || cur_subphase != prev_subphase)
                    {
                        process->report_phase_switch (std::cerr);
                        prev_phase = cur_phase, prev_subphase = cur_subphase;
                    }
                    process->report_state (std::cerr);
                    prevtime = curtime;
                }
            }
        }
    }
    else
    {
        std::cerr << std::endl << "Error while initiating process:" << std::endl << (const char*) process->getError () << std::endl;
        return -1;
    }

    int rv = 0;
    Process::STATE fin_state = process->state ();
    switch (fin_state)
    {
    case Process::APS_COMPLETE:
        if (params->verbose ())
            process->report_final (std::cerr);
        break;
    case Process::APS_INTERRUPTED:
        if (params->verbose ())
        {
            std::cerr << std::endl << "Interrupted by user." << std::endl << std::flush;
        }
        rv = 2;
        break;
    case Process::APS_ERROR:
        std::cerr << std::endl << "Terminated due to error(s):" << std::endl;
        std::cerr << process->getError () << std::endl << std::flush;
        rv = 1;
        break;
    default:
        std::cerr << std::endl << "Abnormal process termination, state = " << fin_state << std::endl << std::flush;
        rv = -2;
    }

	return rv;
}

int main (int argc, char *argv[])
{
    int rv = -1;
    try
    {
        rv = console_main (argc, argv);
    }
    catch (std::bad_alloc& e)
    {
        std::cerr << "Insufficient memory: " << e.what () << std::endl;
    }
    catch (std::exception& e)
    {
        std::cerr << "Error : Standatd Library Exception:" << e.what () << std::endl;
    }
    catch (Rerror& e)
    {
        std::cerr << "Error : " << (const char*) e << std::endl;
    }
#ifndef DONT_CATCH_UNKNOWN
    catch (...)
    {
        std::cerr << "Unhandled exception caught" << std::endl << std::flush;
    }
#endif
    return rv;
}
