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

#include <string>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <sstream>
#ifdef USE_CUSTOM_ALLOCATOR
#include "MemoryCpu.h"
#endif

#ifdef USE_LONG
typedef size_t SizeType;
#else
typedef unsigned int SizeType;
#endif

namespace PsimagLite {

#ifdef USE_CUSTOM_ALLOCATOR
template<typename T,int templateParamFlags> class AllocatorCpu : public std::allocator<T> {
	typedef typename std::allocator<T> BaseType;

	typedef MemoryCpu MemoryCpuType;
	typedef unsigned char ByteType;

public:

	template<typename U>
	struct rebind {
		typedef AllocatorCpu<U,templateParamFlags> other;
	}; // struct rebind

	AllocatorCpu() {}

	// FIXME: needs template constraints here
	template <typename OtherType>
	AllocatorCpu(const OtherType& x)
	    : std::allocator<T>(x)
	{}

	typename BaseType::pointer allocate(typename BaseType::size_type n,void* = 0)
	{
		if (n>this->max_size()) throw
			std::runtime_error("Bad allocation\n");

		typename BaseType::pointer x = (typename BaseType::pointer) globalMemoryCpu.allocate(n*sizeof(T));

		return static_cast<T*>(x);
	}

	void deallocate(typename BaseType::pointer p, typename BaseType::size_type n)
	{
		globalMemoryCpu.deallocate(p);
	}

}; // class AllocatorCpu
#endif

template<typename T>
class Allocator {
public:
#ifdef USE_CUSTOM_ALLOCATOR
	typedef AllocatorCpu<T,1> Type;
#else
	typedef std::allocator<T> Type;
#endif
}; // class Allocator

template<bool b,typename T>
class EnableIf {
};

template<typename T>
class EnableIf<true,T> {
public:
	typedef T Type;
};

template<typename T>
struct RemoveConst {
	typedef T Type;
};

template<typename T>
struct RemoveConst<const T> {
	typedef T Type;
};

typedef std::basic_string<char,std::char_traits<char>,Allocator<char>::Type> String;
typedef std::basic_istringstream<char,std::char_traits<char>,Allocator<char>::Type> IstringStream;
typedef std::basic_ostringstream<char,std::char_traits<char>,Allocator<char>::Type> OstringStream;

class RuntimeError : public std::runtime_error {
public:
  explicit RuntimeError (const String& what_arg)
	    : std::runtime_error(what_arg.c_str())
	{}
};

class RangeError : public std::range_error {
public:
  explicit RangeError(const String& what_arg)
	    : std::range_error(what_arg.c_str())
	{}
};

class LogicError : public std::logic_error {
public:
  explicit LogicError (const String& what_arg)
	    : std::logic_error(what_arg.c_str())
	{}
};

} // namespace PsimagLite

/*@}*/
#endif // ALLOCATOR_CPU_H

