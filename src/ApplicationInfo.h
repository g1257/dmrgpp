/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.]
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
#include "Vector.h"
#include "MersenneTwister.h"
#include <sys/types.h>
#include <unistd.h>
#include "BitManip.h"
#include <cstdlib>
#include <cassert>
#include "Io/IoSerializerStub.h"

namespace PsimagLite {

class ApplicationInfo {

public:

	typedef String RunIdType;

	ApplicationInfo(const PsimagLite::String& name)
	    : name_(name),
	      pid_(getpid()),
	      runId_(runIdInternal()),
	      isFinalized_(false)
	{}
	
	void finalize() { isFinalized_ = true;}

	time_t unixTime(bool arg  = false) const
	{
		struct timeval tv;
		gettimeofday(&tv,0);
		return (arg) ? tv.tv_usec : tv.tv_sec;
	}
	
	String getTimeDate() const
	{
		time_t tt = unixTime();
		return asctime(localtime(&tt));
	}

	String hostname() const
	{
		int len = 1024;
		char* name = new char[len];
		int ret = gethostname(name,len);
		String retString;
		if (ret != 0) {
			retString = "UNKNOWN";
		} else {
			retString = name;
		}

		delete[] name;

		return retString;
	}

	const RunIdType runId() const
	{
		return runId_;
	}

	unsigned int pid() const { return pid_; }

	void write(String label, IoSerializer& serializer) const
	{
		String root = label;
		if (!isFinalized_) {
			serializer.createGroup(root);

			serializer.write(root + "/Name", name_);
			serializer.write(root + "/RunId", runId_);
			serializer.write(root + "/UnixTimeStart", unixTime(false));
		} else {
			serializer.write(root + "/UnixTimeEnd", unixTime(false));
		}
	}

	static void setEnv(String name, String value)
	{
		int ret = setenv(name.c_str(), value.c_str(), true);
		if (ret != 0)
			throw RuntimeError("Could not setenv " + name + "=" + value + "\n");
		std::cout<<"Set "<<name<<"="<<value<<"\n";
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const ApplicationInfo& ai)
	{
		if (ai.isFinalized_)
			printFinalLegacy(os, ai);
		else
			printInit(os, ai);

		return os;
	}

private:

	static void printInit(std::ostream& os, const ApplicationInfo& ai)
	{
		os<<ai.getTimeDate();
		os<<"Hostname: "<<ai.hostname()<<"\n";
		os<<"RunID="<<ai.runId_<<"\n";
		os<<"UnixTimeStart="<<ai.unixTime(false)<<"\n";
	}

	static void printFinalLegacy(std::ostream& os, const ApplicationInfo& ai)
	{
		OstringStream msg(std::cout.precision());
		msg()<<ai.name_<<"\nsizeof(SizeType)="<<sizeof(SizeType)<<"\n";
#ifdef USE_FLOAT
		msg()<<ai.name_<<" using float\n";
#else
		msg()<<ai.name_<<" using double\n";
#endif
		msg()<<"UnixTimeEnd="<<ai.unixTime(false)<<"\n";
		msg()<<ai.getTimeDate();
		os<<msg().str();
	}

	RunIdType runIdInternal() const
	{
		unsigned int p = getpid();
		time_t tt = unixTime(true);
		MersenneTwister mt(tt + p);
		unsigned int x = tt ^ mt.random();
		OstringStream msgg(std::cout.precision());
		OstringStream::OstringStreamType& msg = msgg();
		msg<<x;
		x = p ^ mt.random();
		msg<<x;
		unsigned long int y = atol(msg.str().c_str());
		y ^= mt.random();
		x = BitManip::countKernighan(y);
		OstringStream msgg2(std::cout.precision());
		OstringStream::OstringStreamType& msg2 = msgg2();
		msg2<<y;
		if (x < 10) msg2<<"0";
		msg2<<x;
		return msg2.str();
	}

	long unsigned int convertToLuint(PsimagLite::String str) const
	{
		long unsigned int sum = 0;
		long unsigned int prod = 1;
		int l = str.length();
		assert(l < 20);

		for (int i = 0; i < l; ++i) {
			unsigned int c = str[l - i - 1] - 48;
			sum += prod*c;
			prod *= 10;
		}

		return sum;
	}

	PsimagLite::String name_;
	unsigned int pid_;
	const RunIdType runId_;
	bool isFinalized_;
}; // class ApplicationInfo

std::ostream& operator<<(std::ostream& os,const ApplicationInfo& ai);

} // namespace PsimagLite

/*@}*/
#endif // APPLICATION_INFO_H

