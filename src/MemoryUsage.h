/* Copyright (c) 2009-2013, UT-Battelle, LLC
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
/** \ingroup PsimagLite */
/*@{*/

/*! \file MemoryUsage.h
 *
 */

#ifndef MEMORY_USAGE_H
#define MEMORY_USAGE_H
#include <iostream>
#include <fstream>
#include <time.h>
#include "String.h"

namespace PsimagLite {
class MemoryUsage {

	static const SizeType MY_MAX_LINE = 40240;

public:

	MemoryUsage(const String& myself="")
	    : data_(""),myself_(myself),startTime_(::time(0))
	{
		if (myself_=="") myself_="/proc/self/status";
		update();
	}

	void update()
	{
		std::ifstream ifp(myself_.c_str());
		if (!ifp || !ifp.good() || ifp.bad()) return;
		char tmp[MY_MAX_LINE];
		data_ = "";
		while (!ifp.eof()) {
			ifp.getline(tmp,MY_MAX_LINE);
			data_ += String(tmp);
			data_ += String("\n");
		}

		ifp.close();
	}

	String findEntry(const String& label)
	{
		long unsigned int x = data_.find(label);
		if (x==String::npos) {
			return "NOT_FOUND";
		}
		x += label.length();
		long unsigned int y = data_.find("\n",x);
		SizeType len = y-x;
		if (y==String::npos) len = data_.length()-x;
		String s2 = data_.substr(x,len);
		x = 0;
		for (SizeType i=0;i<s2.length();i++) {
			x++;
			if (s2.at(i)==' ' || s2.at(i)=='\t') continue;
			else break;
		}
		if (x>0) x--;
		len = s2.length()-x;
		return s2.substr(x,len);
	}

	double time() const
	{
		return ::time(0)-startTime_;
	}

private:

	String data_,myself_;
	time_t startTime_;
}; // class MemoryUsage

} // namespace PsimagLite 

/*@}*/	
#endif // MEMORY_USAGE_H

