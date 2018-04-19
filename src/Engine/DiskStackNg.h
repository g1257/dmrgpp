/*
Copyright (c) 2009-2015-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

#ifndef DISKSTACK_NG_H
#define DISKSTACK_NG_H

// All these includes are in PsimagLite
#include "Stack.h"
#include "Io/IoNg.h"
#include "ProgressIndicator.h"
#include <exception>

// A disk stack, similar to std::stack but stores in disk not in memory
namespace Dmrg {
template<typename DataType>
class DiskStack {

	typedef typename PsimagLite::IoNg::In IoInType;
	typedef typename PsimagLite::IoNg::Out IoOutType;

public:

	DiskStack(const PsimagLite::String name1,
	          const PsimagLite::String name2,
	          bool hasLoad,
	          bool isObserveCode)
	    : ioIn_(name1),
	      ioOut_(name2),
	      isObserveCode_(isObserveCode),
	      total_(0),
	      progress_("DiskStack"),
	      dt_(0)
	{
		if (!hasLoad) return;
		err("DiskStackNg does not support reading from file yet\n");
	}

	bool inDisk() const { return true; }

	void push(const DataType& d)
	{
		try {
			d.write(ioOut_, "/" + ttos(total_), DataType::SAVE_ALL);
		} catch (std::exception&) {
			d.overwrite(ioOut_, "/" + ttos(total_), DataType::SAVE_ALL);
		}

		stack_.push(total_++);
	}

	void pop()
	{
		stack_.pop();
	}

	const DataType& top() const
	{
		assert(stack_.size() > 0);
		const SizeType dummy = 0;
		delete dt_;
		dt_ = new DataType(ioIn_, "/" + ttos(stack_.top()), dummy, isObserveCode_);
		return *dt_;
	}

	SizeType size() const { return stack_.size(); }

	void finalize() {}

private:

	mutable IoInType ioIn_;
	IoOutType ioOut_;
	bool isObserveCode_;
	int total_;
	PsimagLite::ProgressIndicator progress_;
	PsimagLite::Stack<int>::Type stack_;
	mutable DataType* dt_;
}; // class DiskStack

} // namespace Dmrg

#endif

