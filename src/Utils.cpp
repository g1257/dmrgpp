/*
Copyright (c) 2009-2016, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

#include "Utils.h"

// Utility functions that are still needed
namespace utils {

PsimagLite::String basename(PsimagLite::String pathname)
{
	return PsimagLite::String(std::find_if(pathname.rbegin(),
	                                       pathname.rend(),
	                                       UnixPathSeparator()).base(),pathname.end());
}

PsimagLite::String pathPrepend(PsimagLite::String pre,PsimagLite::String pathname)
{
	bool addDotDot = false;
	PsimagLite::String path1("");
	if (pathname.length() > 2 && pathname[0] == '.' && pathname[1] == '.') {
		addDotDot = true;
		path1 = pathname.substr(2,pathname.length());
	}

	if (pathname.length() > 1 && pathname[0] == '/') path1 = pathname;

	if (path1 == "") return pre + pathname;

	size_t index = path1.find_last_of("/");

	index++;
	PsimagLite::String ret = path1.substr(0,index) + pre +
	        path1.substr(index,path1.length());
	return (addDotDot) ? ".." + ret : ret;
}

SizeType exactDivision(SizeType a, SizeType b)
{
	SizeType c = static_cast<SizeType>(a/b);
	if (c * b != a)
		throw PsimagLite::RuntimeError("exactDivision expected\n");

	return c;
}

SizeType bitSizeOfInteger(SizeType x)
{
	SizeType counter = 0;
	while (x) {
		counter++;
		x >>= 1;
	}

	return counter;
}

SizeType powUint(SizeType x, SizeType n)
{
	if (n == 0) return 1;
	if (n == 1) return x;
	SizeType ret = x;
	for (SizeType i = 1; i < n; ++i) ret *= x;
	return ret;
}

} //namespace utils

