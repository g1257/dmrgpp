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
#ifndef STRING_HEADER_H
#define STRING_HEADER_H
#include <string>
#include "AllocatorCpu.h"

namespace PsimagLite {

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
#endif // STRING_HEADER_H
