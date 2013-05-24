// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
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
// END LICENSE BLOCK

/** \ingroup PsimagLite */
/*@{*/

/*! \file RandomForTests.h
 * 
 *  THIS RNG IS FOR TESTS, DO NOT USE FOR PRODUCTION, IT'S NOT RANDOM ENOUGH!!
 *
 */

#ifndef RANDOM_FOR_TESTS_H
#define RANDOM_FOR_TESTS_H
#include <cstdlib>
#include <iostream>
#include <string>
#include <stdexcept>

namespace PsimagLite {
	template<typename T>
	class  RandomForTests {
	public:
		typedef long int LongType;
		typedef T value_type; // legacy name
		RandomForTests(int seed) //LongType seed = 127773,SizeType rank = 0,SizeType nprocs = 1)
		: next_(seed)
		{
		}

		T operator()()
		{
			next_ = 16807 * (next_ % 127773) - 2836 * (next_ / 127773);
			if (next_ <= 0) next_ += 2147483647;
			if (next_ >= 2147483647) next_=1;
			return static_cast<T>(next_) / 2147483647.0;
		}
	private:
		int next_;
	}; // RandomForTests
} // namespace PsimagLite

/*@}*/
#endif // RANDOM_FOR_TESTS_H
