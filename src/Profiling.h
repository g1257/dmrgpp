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

namespace PsimagLite {

	class  Profiling {
		public:
			Profiling(const std::string& s) : message_(s),start_(clock())
			{
				std::cout<<"Profiling: Starting clock for: "<<s<<"\n";
			}

			~Profiling()
			{
				double elapsed = clock()-start_;
				elapsed /= CLOCKS_PER_SEC; 
				std::cout<<"Profiling: Stoping clock for: "<<message_<<" elapsed="<<elapsed<<"\n";
			}
		private:
			std::string message_;
			clock_t start_;	
	}; // Profiling
} // PsimagLite

#endif


