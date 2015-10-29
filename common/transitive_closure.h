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

#ifndef __transitive_closure_h__
#define __transitive_closure_h__

#include <vector>
#include <algorithm>
#include "common_typedefs.h"
#include "resource.h"

typedef std::vector <UIntVect> Clusters;

class Closure
{
    enum STAGE { INIT, FILL, COLLECT };
    MemWrapper <unsigned> closure_;
    unsigned sz_;
    unsigned clno_;
    unsigned singlet_no_;
    STAGE stage_;
public:
    Closure ();
    Closure (unsigned sz);
    void init (unsigned sz);
    void add (unsigned idx1, unsigned idx2);
    void finalize ();
    void fillClusters (Clusters& dest, bool include_singletons = false);
    Clusters* collect (bool include_singletons = false);
    unsigned clusters_number () const
    {
        return clno_;
    }
    unsigned size () const
    {
        return sz_;
    }
    unsigned singletons_number () const
    {
        return singlet_no_;
    }
    unsigned non_trivial_number () const
    {
        return sz_ - singlet_no_;
    }
    const unsigned* result ();
};

#endif
