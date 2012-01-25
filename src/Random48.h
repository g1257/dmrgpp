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
#include <iostream>
namespace PsimagLite {
	template<typename T>
	class  Random48 {
	public:
		typedef long int LongType;
		typedef T value_type; // legacy name
		Random48(LongType seed,size_t rank = 0,size_t nprocs = 1)
		{
			srand48(seed);
			std::vector<LongType> vOfSeeds(nprocs);
			for (size_t i=0;i<vOfSeeds.size();i++) vOfSeeds[i] = LongType(10000*random());
			seed_=vOfSeeds[rank];
			srand48(seed_);
			
		}

		T random() // deprecated!!! use operator() instead
		{
			return static_cast<T>(drand48());
		}

		T operator()() const { return static_cast<T>(drand48()); }

		LongType seed() const { return seed_; }

		void seed(const LongType& seed) { seed_ = seed; }

	private:
	
		LongType seed_;
	}; // Random48
} // namespace PsimagLite

#endif // RANDOM48_H

