/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file ProgressIndicator.h
 *
 *  This class handles output to a progress indicator (usually the terminal)
 */

#ifndef PROGRESS_INDICATOR_H
#define PROGRESS_INDICATOR_H

#include <iostream>
#include <sstream>
#include <vector>
#include "Concurrency.h"
#include "MemoryUsage.h"
#include <sys/types.h>
#include <unistd.h>
#include "TypeToString.h"

namespace PsimagLite {

class ProgressIndicator {

	static MemoryUsage musage_;
	static OstringStream buffer_;
	static bool bufferActive_;

public:

	ProgressIndicator(String caller,SizeType threadId = 0)
	    : threadId_(threadId), rank_(0)
	{
		if (threadId_ != 0) return;

		caller_ = caller;
		rank_ = Concurrency::rank();
	}

	static void updateBuffer(int signal)
	{
		if (bufferActive_) {
			pid_t p = getpid();
			String outName("buffer");
			outName += ttos(p);
			outName += ".txt";
			std::ofstream fout(outName.c_str());
			fout<<buffer_.str()<<"\n";
			fout.close();
			buffer_.str("");
		}

		bufferActive_ = !bufferActive_;

		String bufferActive = (bufferActive_) ? "active" : "inactive";
		std::cerr<<"ProgressIndicator: signal "<<signal<<" received.";
		std::cerr<<" buffer is now "<<bufferActive<<"\n";

	}

	template<typename SomeOutputType>
	void printline(const String &s,SomeOutputType& os) const
	{
		if (threadId_ != 0) return;
		if (rank_!=0) return;
		prefix(os);
		os<<s<<"\n";

		if (!bufferActive_) return;

		prefix(buffer_);
		buffer_<<s<<"\n";
	}

	void printline(OstringStream &s,std::ostream& os) const
	{
		if (threadId_ != 0) return;
		if (rank_!=0) return;
		prefix(os);
		os<<s.str()<<"\n";
		s.seekp(std::ios_base::beg);

		if (!bufferActive_) return;

		prefix(buffer_);
		buffer_<<s.str()<<"\n";
		s.seekp(std::ios_base::beg);
	}

	void print(const String& something,std::ostream& os) const
	{
		if (threadId_ != 0) return;
		if (rank_!=0) return;
		prefix(os);
		os<<something;

		if (!bufferActive_) return;

		prefix(buffer_);
		buffer_<<something;
	}

	void printMemoryUsage()
	{
		musage_.update();
		String vmPeak = musage_.findEntry("VmPeak:");
		String vmSize = musage_.findEntry("VmSize:");
		OstringStream msg;
		msg<<"Current virtual memory is "<<vmSize<<" maximum was "<<vmPeak;
		printline(msg,std::cout);

		if (!bufferActive_) return;

		buffer_<<"Current virtual memory is "<<vmSize<<" maximum was "<<vmPeak;
		printline(buffer_,std::cout);
	}

	static MemoryUsage::TimeHandle time() { return musage_.time(); }

private:

	template<typename SomeOutputStreamType>
	void prefix(SomeOutputStreamType& os) const
	{
		const MemoryUsage::TimeHandle t = musage_.time();
		const double seconds = t.millis();
		const SizeType prec = os.precision(3);
		os<<caller_<<" "<<"["<<std::fixed<<seconds<<"]: ";
		os.precision(prec);
	}

	String caller_;
	SizeType threadId_;
	SizeType rank_;
}; // ProgressIndicator

} // namespace PsimagLite

/*@}*/
#endif

