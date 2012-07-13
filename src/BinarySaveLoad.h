/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file BinarySaveLoad.h
 *
 *  This class contains saves and loads for native types and stl classes
 *  Don't add here saves and adds for custom classes, those need to be
 *  added to the class itself as companions
 *
 */
  
#ifndef BINARY_LOAD_SAVE_H
#define BINARY_LOAD_SAVE_H

#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include "TypeToString.h"

namespace PsimagLite {

class BinarySaveLoad {

public:

	static void save(int fd,int x)
	{
		mywrite(fd,(const void *)&x,sizeof(x));
	}

	static void save(int fd,size_t x)
	{
		mywrite(fd,(const void *)&x,sizeof(x));
	}

	static void save(int fd,const std::vector<int>& vec)
	{
		vsave_(fd,vec);
	}

	static void save(int fd,const std::vector<size_t>& vec)
	{
		vsave_(fd,vec);
	}

	static void save(int fd,const std::vector<double>& vec)
	{
		vsave_(fd,vec);
	}

	static void save(int fd,const std::vector<float>& vec)
	{
		vsave_(fd,vec);
	}

	static void save(int fd,const std::vector<std::complex<double> >& vec)
	{
		vsave_(fd,vec);
	}

	template<typename NonNativeType>
	static void save(int fd,const NonNativeType& composed)
	{
		composed.save(fd);
	}

	static void mywrite(int fd,const void *buf,size_t count)
	{
		ssize_t ret = write(fd,buf,count);
		failIfNegative(ret,__FILE__,__LINE__);
	}

private:

	template<typename NativeType>
	static void vsave_(int fd,const std::vector<NativeType>& vec)
	{
		size_t length = vec.size();
		mywrite(fd,(const void *)&length,sizeof(length));

		const NativeType dummy = 0;

		for (size_t i=0;i<vec.size();i++) {
			const NativeType* const ptr = &(vec[i]);
			mywrite(fd, (const void *)ptr,sizeof(dummy));
		}
	}


	static void failIfNegative(const ssize_t& x,const std::string& thisFile,int lineno)
	{
		if (x>=0) return;
		std::string str(thisFile);
		str += " " + ttos(lineno) + "\n";
		str += "Read or Write has failed\n";
		throw std::runtime_error(str.c_str());
	}
}; // BinarySaveLoad


} // namespace PsimagLite

/*@}*/	
#endif
