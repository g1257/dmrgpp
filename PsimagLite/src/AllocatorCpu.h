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

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifdef USE_CUSTOM_ALLOCATOR
#include "MemoryCpu.h"
#endif

#ifndef USE_SHORT
using SizeType = long unsigned int;
#else
using SizeType = uint32_t;
#endif

namespace PsimagLite {

template <typename T, bool B = std::is_enum<T>::value> struct IsEnumClass : std::false_type { };

template <typename T>
struct IsEnumClass<T, true>
    : std::integral_constant<
          bool,
          !std::is_convertible<T, typename std::underlying_type<T>::type>::value> { };

#ifdef USE_CUSTOM_ALLOCATOR
template <typename T, int templateParamFlags> class AllocatorCpu : public std::allocator<T> {
	using BaseType = typename std::allocator<T>;

	using MemoryCpuType = MemoryCpu;
	using ByteType      = unsigned char;

public:

	template <typename U> struct rebind {
		using other = AllocatorCpu<U, templateParamFlags>;
	}; // struct rebind

	AllocatorCpu() { }

	// FIXME: needs template constraints here
	template <typename OtherType>
	AllocatorCpu(const OtherType& x)
	    : std::allocator<T>(x)
	{ }

	typename BaseType::pointer allocate(typename BaseType::size_type n, void* = 0)
	{
		if (n > this->max_size())
			throw std::runtime_error("Bad allocation\n");

		typename BaseType::pointer x
		    = (typename BaseType::pointer)globalMemoryCpu.allocate(n * sizeof(T));

		return static_cast<T*>(x);
	}

	void deallocate(typename BaseType::pointer p, typename BaseType::size_type n)
	{
		globalMemoryCpu.deallocate(p);
	}

}; // class AllocatorCpu
#endif

template <typename T> class Allocator {
public:

#ifdef USE_CUSTOM_ALLOCATOR
	using Type = AllocatorCpu<T, 1>;
#else
	using Type = std::allocator<T>;
#endif
}; // class Allocator

template <bool b, typename T> class EnableIf { };

template <typename T> class EnableIf<true, T> {
public:

	using Type = T;
};

template <typename T> struct RemoveConst {
	using Type = T;
};

template <typename T> struct RemoveConst<const T> {
	using Type = T;
};

template <typename T> struct IsStringLike {
	enum
	{
		True = false
	};
};

template <typename A> struct IsStringLike<std::basic_string<char, std::char_traits<char>, A>> {
	enum
	{
		True = true
	};
};

using String        = std::basic_string<char, std::char_traits<char>, Allocator<char>::Type>;
using IstringStream = std::basic_istringstream<char, std::char_traits<char>, Allocator<char>::Type>;

class OstringStream {
public:

	typedef std::basic_ostringstream<char, std::char_traits<char>, Allocator<char>::Type>
	    OstringStreamType;

	OstringStream(int prec) { data_.precision(prec); }

	int precision(int prec) { return data_.precision(prec); }

	OstringStreamType& operator()() { return data_; }

private:

	OstringStreamType data_;
};

class RuntimeError : public std::runtime_error {
public:

	explicit RuntimeError(const String& what_arg)
	    : std::runtime_error(what_arg.c_str())
	{ }
};

class RangeError : public std::range_error {
public:

	explicit RangeError(const String& what_arg)
	    : std::range_error(what_arg.c_str())
	{ }
};

class LogicError : public std::logic_error {
public:

	explicit LogicError(const String& what_arg)
	    : std::logic_error(what_arg.c_str())
	{ }
};

} // namespace PsimagLite

/*@}*/
#endif // ALLOCATOR_CPU_H
