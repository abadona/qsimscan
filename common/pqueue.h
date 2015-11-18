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

#ifndef __pqueue_h__
#define __pqueue_h__

// This is a customised version of heap-based priority queue
// with limited capacity. It uses minimum based PQ rather than
// more usual maximum to be able to remove the least element when
// queue is full and larger element to be inserted.

// Elements should implement operator< and operator== and have operator=
// Default operator= often does right

#include <algorithm>
#include <vector>

class PQueueError {};
class PQueueUnderflow:   public PQueueError {};
class PQueueNotAHeap:    public PQueueError {};
class PQueueRestructure: public PQueueError {};

template <class T>
bool greater (T const& v1, T const& v2)
{
    return v2 < v1;
}

template <class T>
class PQueue
{
public:
    typedef T Elem;
    typedef std::vector <T> ElemVec;
private:

	class GreaterType
	{
	public:
		bool operator () (const T& a, const T& b) const
		{
			if (a < b)
				return false;
			else
				return true;
		}
	};
	GreaterType Greater;

    ElemVec     values_;
    unsigned    n_ : 31;
    bool        heap_ok_ : 1;

    void heapify ()
    {
        std::make_heap (values_.begin (), values_.begin () + n_, Greater);
        heap_ok_ = true;
    }

public:
    PQueue() // default constructor - enabling use of containers on PQueues
    :
    n_ (0),
    heap_ok_ (true)
    {
    }
    explicit PQueue (unsigned maxN)
    :
    n_ (0),
    heap_ok_ (true)
    {
        values_.resize (maxN);
    }
    T const& top () const throw (PQueueUnderflow, PQueueNotAHeap)
    {
        if (!heap_ok_) throw PQueueNotAHeap ();
        if (!n_) throw PQueueUnderflow ();
        return values_.front ();
    }
    bool push (T const& val) throw (PQueueNotAHeap)
    {
        bool toR = false;
        if (!heap_ok_) throw PQueueNotAHeap ();
        if (n_ == values_.size ())
        {
            // Queue full, so we will insert element only
            // if the smallest element (at head) is smaller then it.
            if (values_.front () < val)
            {
                std::pop_heap (values_.begin (), values_.end (), Greater);
                values_.back () = val;
                std::push_heap (values_.begin (), values_.end (), Greater);
                toR = true;
            }
        }
        else
        {
            values_ [n_++] = val;
            std::push_heap (values_.begin (), values_.begin () + n_, Greater);
            toR = true;
        }
        return toR;
    }
    void pop () throw (PQueueUnderflow, PQueueNotAHeap)
    {
        if (!heap_ok_) throw PQueueNotAHeap ();
        if (!n_) throw PQueueUnderflow ();
        std::pop_heap (values_.begin (), values_.begin () + n_--, Greater);
    }
    bool empty ()
    {
        return !n_;
    }
    unsigned size () const
    {
        return n_;
    }
    unsigned capacity () const
    {
        return values_.size ();
    }
    void setCapacity (int maxN) throw (PQueueRestructure)
    {
        if (n_ > maxN) throw PQueueRestructure ();
        values_.resize (maxN);
        heapify ();
    }
    ElemVec& sort ()
    {
        // this is not sort_heap, just sort - because different predicate then heapifying one is used
        values_.resize (n_);
        // std::sort (values_.begin (), values_.end ());
        std::sort_heap (values_.begin (), values_.end (), Greater);
        heap_ok_ = false;
        return values_;
    }
    const ElemVec& data () const
    {
        return values_;
    }
};

#endif // __pqueue_h__
