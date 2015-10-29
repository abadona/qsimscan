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

#include "sync_control.h"

#if defined(_WIN32) && defined(_MSC_VER)
#include <process.h>
#endif


static const int default_stack_size = 0x200000; // 2 Mb


SyncControl::SyncControl ()
{
#if defined(_WIN32) && defined(_MSC_VER)

    InitializeCriticalSection (&critical_section_);

#else

    pthread_mutex_init (&critical_section_, NULL);
    pthread_mutex_init (&term_mutex_, NULL);
    pthread_cond_init (&exit_condition_, NULL);


#endif

    terminated_ = false;
}

SyncControl::~SyncControl ()
{
#if defined(_WIN32) && defined(_MSC_VER)

    DeleteCriticalSection (&critical_section_);

#else

    pthread_mutex_destroy (&critical_section_);
    pthread_mutex_destroy (&term_mutex_);
    pthread_cond_destroy (&exit_condition_);

#endif
}

#if defined(_WIN32) && defined(_MSC_VER)

#else

void* tfwrapper (void* params)
{
    TFWrapper_params* p = (TFWrapper_params*) params;
    // this assumes tha thread_func will not terminate the thread by itself, but will always either or throw exception
    try
    {
        p->thread_func (p->params);
    }
    catch (...)
    {
        pthread_mutex_lock (p->term_mutex);
        *p->term_flag = true;
        pthread_mutex_unlock (p->term_mutex);
        pthread_cond_signal (p->cond);

        throw;
    }


    pthread_mutex_lock (p->term_mutex);
    *p->term_flag = true;
    pthread_mutex_unlock (p->term_mutex);
    pthread_cond_signal (p->cond);

    return NULL;
}
#endif

void SyncControl::begin		(void (*thread_func) (void*), void* arg)
{
#if defined(_WIN32) && defined(_MSC_VER)

    processing_thread_ = _beginthread (thread_func, default_stack_size, arg);

#else

    wrapper_params_.thread_func = thread_func;
    wrapper_params_.params = arg;
    wrapper_params_.cond = &exit_condition_;
    wrapper_params_.term_mutex = &term_mutex_;
    wrapper_params_.term_flag = &terminated_;

    pthread_create (&processing_thread_, NULL, tfwrapper, &wrapper_params_);

#endif
}

bool SyncControl::wait (unsigned int msec)
{
#if defined(_WIN32) && defined(_MSC_VER)

    return WAIT_OBJECT_0 == WaitForSingleObject ((HANDLE) processing_thread_, msec);

#else

    // wait for conditional variable for desired time;
    // if signalled, enter thread waiting state. otherwise fail

    timespec tspec;
    tspec.tv_sec = time (NULL) + msec / 1000;
    tspec.tv_nsec = (msec % 1000) * 1000;

    pthread_mutex_lock (&term_mutex_);
    while (!terminated_)
    {
        if (0 != pthread_cond_timedwait (&exit_condition_, &term_mutex_, &tspec))
            break;
    }
    bool toR = terminated_;
    pthread_mutex_unlock (&term_mutex_);

    return toR;

#endif
}

bool SyncControl::isRunning	()
{
#if defined(_WIN32) && defined(_MSC_VER)

    return !terminated_;

#else

    pthread_mutex_lock (&term_mutex_);
    bool toR = terminated_;
    pthread_mutex_unlock (&term_mutex_);

    return toR;

#endif
}

void SyncControl::lock		()
{
#if defined(_WIN32) && defined(_MSC_VER)

    EnterCriticalSection (&critical_section_);

#else

    pthread_mutex_lock (&critical_section_);


#endif
}

void SyncControl::unlock	()
{
#if defined(_WIN32) && defined(_MSC_VER)

    LeaveCriticalSection (&critical_section_);

#else

    pthread_mutex_unlock (&critical_section_);

#endif
}
