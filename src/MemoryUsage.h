// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."
 
*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
// END LICENSE BLOCK
/** \ingroup PsimagLite */
/*@{*/

/*! \file MemoryUsage.h
 *
 * getrusage system call
 */
  
#ifndef MEMORY_USAGE_H
#define MEMORY_USAGE_H
#include <iostream>
#include <fstream>
#include <string>

namespace PsimagLite {
	class MemoryUsage {
		static const char *MY_SELF_FILE;
		static const size_t MY_MAX_LINE = 40240;
	public:
		MemoryUsage() : data_("")
		{
			update();
		}
		
		void update()
		{
			std::ifstream ifp(MY_SELF_FILE);
			char tmp[MY_MAX_LINE];
			data_ = "";
			while (!ifp.eof()) {
				ifp.getline(tmp,MY_MAX_LINE);
				data_ += std::string(tmp);
			}
			ifp.close();
		}
		
		long vmSize(bool needsUpdate = true)
		{
			if (needsUpdate) update();
			return findEntry("VmSize:");
		}

		long vmPeak(bool needsUpdate = true)
		{
			if (needsUpdate) update();
			return findEntry("VmPeak:");
		}
		
	private:
		long findEntry(const std::string& label)
		{
			size_t x = data_.find(label);
			if (x==std::string::npos) {
				std::string s = "MemoryUsage::findInString(...) failed for label=" + label + "\n";
				throw std::runtime_error(s.c_str());
			}
			size_t y = data_.find(" ",x);
			if (y==std::string::npos) {
				std::string s = "MemoryUsage::findInString(...) no value for label=" + label + "\n";
				throw std::runtime_error(s.c_str());
			}
			std::string buffer = "";
			for (size_t i=y;i<data_.length();i++) {
				char c = data_[i];
				if (c==' ') continue;
				if (c=='.') {
					buffer += data_[i];
					continue;
				}
				if (c>=48 && c<58) {
					buffer += data_[i];
					continue;
				}
				break;
			}
			return atol(buffer.c_str());
		}

		std::string data_;
	}; // class MemoryUsage

	const char *MemoryUsage::MY_SELF_FILE = "/proc/self/status";
} // namespace PsimagLite 

/*@}*/	
#endif // MEMORY_USAGE_H
