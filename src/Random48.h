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
#ifndef RANDOM48_H
#define RANDOM48_H
#include <cstdlib>

namespace PsimagLite {
	template<typename T>
	class  Random48 {
	public:
		typedef long int LongType;
		typedef T value_type; // legacy name
		static void seed(LongType seed)
		{
			srand48(seed);
		}

		static T random()
		{
			return drand48();
		}
	}; // Random48
} // namespace PsimagLite

#endif // RANDOM48_H

