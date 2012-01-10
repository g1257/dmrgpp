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

// A class to profile a scope of code
#ifndef PROFILING_H_
#define PROFILING_H_

#include <string>
#include <iostream>
#include "Rusage.h"

namespace PsimagLite {

		std::ostream& operator<<(std::ostream& os,const std::pair<size_t,size_t>& p)
		{

			os<<p.first<<" "<<p.second<<" ";
			return os;
		}

	class  Profiling {

		typedef std::pair<size_t,size_t> PairType;

		double diff(const PairType& p1,const PairType& p2) const
		{
			size_t sec = p1.first -p2.first;
			size_t usec = p1.second - p2.second;
			return double(sec) + double(usec)/1000000.0; 
		}

		public:
			Profiling(const std::string& s,std::ostream& os = std::cout) 
			: message_(s),start_(rusage_.time()),os_(os)
			{
				os_<<"Profiling: Starting clock for "<<s<<"\n";
			}

			~Profiling()
			{
				PairType end = rusage_.time();
				double elapsed = diff(end,start_);
				os_<<"Profiling: Stoping clock for "<<message_<<" elapsed="<<elapsed<<"\n";
				std::cout<<" start="<<start_<<" end="<<end<<"\n";
			}
		private:
			Rusage rusage_;
			std::string message_;
			PairType start_;	
			std::ostream& os_;
	}; // Profiling
} // PsimagLite

#endif


