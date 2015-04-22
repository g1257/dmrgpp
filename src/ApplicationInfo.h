/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file HostInfo.h
 *
 * Information about the host computer
 */

#ifndef APPLICATION_INFO_H
#define APPLICATION_INFO_H

#include <sys/time.h>
#include <time.h>

namespace PsimagLite {

class ApplicationInfo {

public:

	ApplicationInfo(const PsimagLite::String& name)
	    : firstCall_(true),name_(name)
	{
		std::cout<<(*this);
	}

	~ApplicationInfo()
	{
		std::cout<<getTimeDate();
	}

	String getTimeDate() const
	{
		struct timeval tv;
		time_t tt;

		gettimeofday(&tv,0);
		tt=tv.tv_sec; /* seconds since 1970 */
		return asctime(localtime(&tt));
	}

	friend std::ostream& operator<<(std::ostream& os,const ApplicationInfo& ai);

private:

	mutable bool firstCall_;
	PsimagLite::String name_;
}; // class ApplicationInfo

std::ostream& operator<<(std::ostream& os,const ApplicationInfo& ai)
{
	os<<ai.getTimeDate();
	if (!ai.firstCall_) return os;
	os<<ai.name_<<" sizeof(SizeType)="<<sizeof(SizeType)<<"\n";
#ifdef USE_FLOAT
	os<<ai.name_<<" using float\n";
#else
	os<<ai.name_<<" using double\n";
#endif
	return os;
}

} // namespace PsimagLite

/*@}*/
#endif // APPLICATION_INFO_H

