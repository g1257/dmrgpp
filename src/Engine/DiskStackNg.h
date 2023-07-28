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
#include "BasisTraits.hh"
#include "Io/IoNg.h"
#include "ProgressIndicator.h"
#include "Stack.h"
#include <exception>

// A disk stack, similar to std::stack but stores in disk not in memory
namespace Dmrg
{
template <typename DataType>
class DiskStack
{

	typedef typename PsimagLite::IoNg::In IoInType;
	typedef typename PsimagLite::IoNg::Out IoOutType;

public:

	DiskStack(const PsimagLite::String filename,
	    bool needsToRead,
	    PsimagLite::String label,
	    const BasisTraits& basisTraits)
	    : ioOut_((needsToRead) ? 0 : new IoOutType(filename, PsimagLite::IoNg::ACC_RDW))
	    , ioIn_((needsToRead) ? new IoInType(filename) : 0)
	    , label_("DiskStack" + label)
	    , basisTraits_(basisTraits)
	    , total_(0)
	    , progress_("DiskStack")
	    , dt_(0)
	{
		if (!needsToRead) {
			ioOut_->createGroup(label_);
			ioOut_->write(total_, label_ + "/Size");
			return;
		}

		ioIn_->read(total_, label_ + "/Size");
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Read from file " + filename + " succeeded";
		progress_.printline(msgg, std::cout);
	}

	~DiskStack()
	{
		delete dt_;
		dt_ = 0;
		delete ioIn_;
		ioIn_ = 0;
		delete ioOut_;
		ioOut_ = 0;
	}

	void flush()
	{
		assert(ioOut_);
		ioOut_->flush();
	}

	bool inDisk() const { return true; }

	void push(const DataType& d)
	{
		assert(ioOut_);

		try {
			d.write(*ioOut_,
			    label_ + "/" + ttos(total_),
			    IoOutType::Serializer::NO_OVERWRITE,
			    DataType::SaveEnum::ALL);
		} catch (...) {
			d.write(*ioOut_,
			    label_ + "/" + ttos(total_),
			    IoOutType::Serializer::ALLOW_OVERWRITE,
			    DataType::SaveEnum::ALL);
		}

		++total_;

		ioOut_->write(total_,
		    label_ + "/Size",
		    IoOutType::Serializer::ALLOW_OVERWRITE);
	}

	void pop()
	{
		if (total_ == 0)
			err("Can't pop; the stack is empty!\n");

		--total_;

		if (!ioOut_)
			return;

		ioOut_->write(total_,
		    label_ + "/Size",
		    IoOutType::Serializer::ALLOW_OVERWRITE);
	}

	void restore(SizeType total)
	{
		total_ = total;
		if (!ioOut_)
			return;

		ioOut_->write(total_,
		    label_ + "/Size",
		    IoOutType::Serializer::ALLOW_OVERWRITE);
	}

	const DataType& top() const
	{
		if (!ioIn_)
			err("DiskStack::top() called with ioIn_ as nullptr\n");

		assert(total_ > 0);
		delete dt_;
		dt_ = 0;
		dt_ = new DataType(*ioIn_,
		    label_ + "/" + ttos(total_ - 1),
		    basisTraits_);
		return *dt_;
	}

	SizeType size() const { return total_; }

private:

	DiskStack(const DiskStack&);

	DiskStack& operator=(const DiskStack&);

	IoOutType* ioOut_;
	IoInType* ioIn_;
	PsimagLite::String label_;
	const BasisTraits& basisTraits_;
	int total_;
	PsimagLite::ProgressIndicator progress_;
	mutable DataType* dt_;
}; // class DiskStack

} // namespace Dmrg

#endif
