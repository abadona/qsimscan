
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
// For any questions please contact SciDM team by email at team@scidm.org
//////////////////////////////////////////////////////////////////////////////

#ifndef __sync_control_h__
#define __sync_control_h__

#if defined(_WIN32) && defined(_MSC_VER)
#include <windows.h>
#else
#include <pthread.h>
#include <time.h>
#include <cstring>
#include <cerrno>
#endif

#include "rerror.h"
#include "tracer.h"

#if defined(_WIN32) && defined(_MSC_VER)
typedef CRITICAL_SECTION MUTEX;
#else
typedef pthread_mutex_t MUTEX;
#endif

// The SyncControl class serves one thread / one lock scenario only!

#if defined(_WIN32) && defined(_MSC_VER)

#else
struct TFWrapper_params
{
    void             (*thread_func) (void *);
    void*            params;
    pthread_cond_t*  cond;
    pthread_mutex_t* term_mutex;
    bool*            term_flag;
};
#endif

class SyncControl
{
    #if defined(_WIN32) && defined(_MSC_VER)
        unsigned long        processing_thread_;
        CRITICAL_SECTION     critical_section_;
    #else
        TFWrapper_params     wrapper_params_;
        pthread_t            processing_thread_;
        pthread_mutex_t      term_mutex_;
        pthread_mutex_t      critical_section_;
        pthread_cond_t       exit_condition_;
    #endif
        bool                 terminated_;

public:
    SyncControl     ();
    ~SyncControl    ();
    void begin      (void (*thread_func) (void*), void* arg);
    bool wait       (unsigned int msec);
    bool isRunning  ();
    void lock       ();
    void unlock     ();
};

/////////////////////////////////////////////////////////////////////
// Wrapper for a in-process mutex ('critical section')
class Lock
{
        MUTEX mutex_;
        MUTEX* mutex_ptr_;
        bool own_;
    public:
        Lock ();
        Lock (MUTEX& mutex);
        ~Lock ();
        void lock ();
        bool try_lock ();
        void unlock ();
        void unlock_if_locked ();
        MUTEX* get_mutex ()
        {
            return &mutex_;
        };
        bool locked_;
    friend class Engage;
};

/////////////////////////////////////////////////////////////////////
// RIA mechanism for a Lock instance to ensure unlocking on scope exit
// !!! Engage is not thread-safe itself - it should never be shared!
class Engage
{
        Lock& lock_;
        bool armed_;
    public:
        Engage (Lock& lock)
        :
        lock_ (lock),
        armed_ (false)
        {
            lock_.lock ();
            armed_ = true;
        }
        void disarm ()
        {
            if (armed_)
            {
                lock_.unlock ();
                armed_ = false;
            }
        }
        ~Engage ()
        {
            disarm ();
        }
};

class Notification
{
#if defined(_WIN32) && defined(_MSC_VER)
    #error Windows Signalling not implemented
    // TODO: implement using windows thread and CreateEvent / WaitFor{Single/Multiple}Objects / SetEvent mechanism
#else
        unsigned wait_counter_;
        pthread_cond_t alarm_;
#endif
        Lock alarm_lock_;

    public:
        Notification ();
        Notification (Lock& alarm);
        ~Notification ();
        void lock (); // acquires exclusivity for conditions checking. HAS TO be called before 'wait';
        void unlock (); // gives up exclusivity; HAS TO be called after 'wait' (as it acquires exclusivity before returning)
        void wait (); // waits for notification to arrive from other thread. 'lock' should be called befor 'wait'.
        bool timed_wait (unsigned nsec); // waits for notification to arrive from other thread. 'lock' should be called befor 'wait'.. Returns false on timeout
        void send (); // sends notification to one waiting thread
        void broadcast (); // broadcasts notification to all waiting threads
        unsigned number_waiting (); // returns number of threads waiting for notification. Calling in 'locked' state is safe (result may not change bwefore unlock)
        Lock& get_lock () { return alarm_lock_; }
};

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// Implementations

/////////////////////////////////////////////////////////////////////
// Lock
/////////////////////////////////////////////////////////////////////

inline Lock :: Lock ()
:
mutex_ptr_ (&mutex_),
own_ (true),
locked_ (false)
{
#if defined(_WIN32) && defined(_MSC_VER)
    InitializeCriticalSection (&mutex_);
#else
    int r = pthread_mutex_init (&mutex_, NULL);
    if (r)
        ers << "Mutex init error: " << strerror (r) << Throw;
#endif
}

inline Lock :: Lock (MUTEX& ext_mutex)
:
mutex_ptr_ (&ext_mutex),
own_ (false),
locked_ (false)
{
}

inline Lock :: ~Lock ()
{
    if (own_)
    {
#if defined(_WIN32) && defined(_MSC_VER)
        DeleteCriticalSection (&mutex_);
#else
        int r = pthread_mutex_destroy (&mutex_);
        if (r)
        {
            errlog << "Mutex destroy error: " << strerror (r) << std::endl;
        }
#endif
    }
}

inline bool Lock :: try_lock ()
{
#if defined(_WIN32) && defined(_MSC_VER)
    return TryEnterCriticalSection (mutex_ptr_);
#else
    int r = pthread_mutex_trylock (mutex_ptr_);
    if (r && r != EBUSY)
        ers << "Mutex try_lock error: " << strerror (r) << Throw;
     locked_ = (!r || r != EBUSY);
     return !r;
#endif
}

inline void Lock :: lock ()
{
#ifdef LOCKING_TRACEOUT
    trc << "Lock("<< mutex_ptr_ <<")::lock" << std::endl;
#endif
#if defined(_WIN32) && defined(_MSC_VER)
    EnterCriticalSection (mutex_ptr_);
#else
    int r = pthread_mutex_lock (mutex_ptr_);
    if (r)
        ers << "Mutex lock error: " << strerror (r) << Throw;
#endif
    locked_ = true;
}

inline void Lock :: unlock ()
{
#ifdef LOCKING_TRACEOUT
    trc << "Lock("<< mutex_ptr_ <<")::unlock" << std::endl;
#endif
#if defined(_WIN32) && defined(_MSC_VER)
    LeaveCriticalSection (mutex_ptr_);
#else
    int r = pthread_mutex_unlock (mutex_ptr_);
    if (r)
        ers << "Mutex unlock error: " << strerror (r) << Throw;
#endif
    locked_ = false;
}

inline void Lock :: unlock_if_locked ()
{
    if (locked_)
        unlock ();
}

/////////////////////////////////////////////////////////////////////
// Notification
/////////////////////////////////////////////////////////////////////

inline Notification :: Notification ()
:
wait_counter_ (0)
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    int r = pthread_cond_init (&alarm_, NULL);
    if (r)
        ers << "Conditional init error: " << strerror (r) << Throw;
#endif
}

inline Notification :: Notification (Lock& alarm)
:
wait_counter_ (0),
alarm_lock_ (*alarm.get_mutex ())
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    int r = pthread_cond_init (&alarm_, NULL);
    if (r)
        ers << "Conditional init error: " << strerror (r) << Throw;
#endif
}

inline Notification :: ~Notification ()
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    int r = pthread_cond_destroy (&alarm_);
    if (r)
        ers << "Conditional destroy error: " << strerror (r) << Throw;
#endif
}

inline void Notification :: lock ()
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    alarm_lock_.lock ();
#endif
}

inline void Notification :: unlock ()
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    alarm_lock_.unlock ();
#endif
}


inline void Notification :: wait ()
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    // assumes locked state, otherwise disaster pending. (Testing for that not always trivial, so it is not done to avoid complication)
    wait_counter_ ++;
    alarm_lock_.locked_ = false; // this is not synchronized. Debug puroposes only
#ifdef LOCKING_TRACEOUT
    trc << "Unlocking " << alarm_lock_.get_mutex () << " due to cond_wait call" << std::endl;
#endif
    int r = pthread_cond_wait (&alarm_, alarm_lock_.get_mutex ());
    wait_counter_ --;
    if (r)
        ers << "Notification wait error: " << strerror (r) << Throw;
    else
        alarm_lock_.locked_ = true;
#endif
}

inline bool Notification :: timed_wait (unsigned nsec)
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    // assumes locked state, otherwise disaster pending. (Testing for that not always trivial, so it is not done to avoid complication)
    wait_counter_ ++;
    alarm_lock_.locked_ = false; // this is not synchronized. Debug puroposes only
#ifdef LOCKING_TRACEOUT
    trc << "Unlocking " << alarm_lock_.get_mutex () << " due to cond_wait call" << std::endl;
#endif
    struct timespec abstime;
    clock_gettime (CLOCK_REALTIME, &abstime);
    abstime.tv_nsec += nsec;
    if (abstime.tv_nsec >= 1000000000)
    {
        abstime.tv_sec += abstime.tv_nsec / 1000000000;
        abstime.tv_nsec %= 1000000000;
    }
    int r = pthread_cond_timedwait (&alarm_, alarm_lock_.get_mutex (), &abstime);
    wait_counter_ --;
    if (r &&  (r != ETIMEDOUT))
        ers << "Notification wait error: " << strerror (r) << Throw;
    else
        alarm_lock_.locked_ = true;
    return r != ETIMEDOUT;
#endif
}

inline void Notification :: send ()
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    int r = pthread_cond_signal (&alarm_);
    if (r)
        ers << "Notification send error: " << strerror (r) << Throw;
#endif
}

inline void Notification :: broadcast ()
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    int r = pthread_cond_broadcast (&alarm_);
    if (r)
        ers << "Notification broadcast error: " << strerror (r) << Throw;
#endif
}

inline unsigned Notification :: number_waiting ()
{
#if defined(_WIN32) && defined(_MSC_VER)
#else
    return wait_counter_;
#endif
}

#endif
