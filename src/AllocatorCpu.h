/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
/** \ingroup PsimagLite */
/*@{*/

/*!
 *
 *
 */
#ifndef ALLOCATOR_CPU_H
#define ALLOCATOR_CPU_H

#include <vector>
#include <stdexcept>
#include <cstdlib>
#include "MemoryCpu.h"

namespace PsimagLite {

template<typename T,int templateParamFlags> class AllocatorCpu : public std::allocator<T> {
	typedef typename std::allocator<T> BaseType;

	typedef MemoryCpu MemoryCpuType;
	typedef unsigned char ByteType;

public:

	template <class U> struct rebind { typedef AllocatorCpu<U,templateParamFlags>
									   other; };

	typename BaseType::pointer allocate(typename BaseType::size_type n,void* = 0)
	{
		if (n>this->max_size()) throw
			std::runtime_error("Bad allocation\n");

		typename BaseType::pointer x = (typename BaseType::pointer) globalMemoryCpu.allocate(n*sizeof(T));

//		x = (typename BaseType::pointer) malloc(n*sizeof(T));
		return static_cast<T*>(x);
	}

	void deallocate(typename BaseType::pointer p, typename BaseType::size_type n)
	{
		globalMemoryCpu.deallocate(p);
	}

}; // class AllocatorCpu

template<typename T>
class Allocator {
public:
#ifdef USE_CUSTOM_ALLOCATOR
	typedef AllocatorCpu<T,1> Type;
#else
	typedef std::allocator<T> Type;
#endif
}; // class Allocator

} // namespace PsimagLite

/*@}*/
#endif // ALLOCATOR_CPU_H
