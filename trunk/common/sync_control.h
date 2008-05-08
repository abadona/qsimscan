
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

#ifndef __sync_control_h__
#define __sync_control_h__


#if defined(_WIN32) && defined(_MSC_VER)
#include <windows.h>
#else
#include <pthread.h>
#endif


// this class serves one thread / one lock scenario only!

#if defined(_WIN32) && defined(_MSC_VER)

#else
struct TFWrapper_params
{
	void		   (*thread_func) (void *);
	void*			 params;
	pthread_cond_t*  cond;
	pthread_mutex_t* term_mutex;
	bool*			 term_flag;
};
#endif


class SyncControl
{
	#if defined(_WIN32) && defined(_MSC_VER)

		unsigned long		processing_thread_;
		CRITICAL_SECTION	critical_section_;

	#else

		TFWrapper_params	wrapper_params_;
		pthread_t			processing_thread_;
		pthread_mutex_t		term_mutex_;
		pthread_mutex_t		critical_section_;
		pthread_cond_t		exit_condition_;

	#endif

		bool				terminated_;

public:
	SyncControl		();
	~SyncControl	();
	void begin		(void (*thread_func) (void*), void* arg);
	bool wait		(unsigned int msec);
	bool isRunning	();

	void lock		();
	void unlock		();
};


#endif
