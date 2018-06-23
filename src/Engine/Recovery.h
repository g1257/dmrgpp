/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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
/** \ingroup DMRG */
/*@{*/

/*! \file Recovery.h
 *
 *
 */

#ifndef DMRG_RECOVER_H
#define DMRG_RECOVER_H

#include "Checkpoint.h"
#include "Vector.h"
#include "ProgramGlobals.h"
#include "ProgressIndicator.h"
#include <fstream>

namespace Dmrg {

template<typename ParametersType,typename TargetingType>
class Recovery  {

public:

	enum {SYSTEM = ProgramGlobals::SYSTEM, ENVIRON = ProgramGlobals::ENVIRON};

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef Checkpoint<ParametersType,TargetingType> CheckpointType;
	typedef typename TargetingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename TargetingType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename CheckpointType::IoType IoType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename CheckpointType::MemoryStackType MemoryStackType;
	typedef typename CheckpointType::DiskStackType DiskStackType;

	Recovery(const CheckpointType& checkpoint,
	         const WaveFunctionTransfType& wft,
	         const BasisWithOperatorsType& pS,
	         const BasisWithOperatorsType& pE)
	    : progress_("Recovery"),
	      checkpoint_(checkpoint),
	      wft_(wft),
	      pS_(pS),
	      pE_(pE),
	      counter_(0)
	{}

	~Recovery()
	{
		if (checkpoint_.parameters().options.find("recoveryNoDelete") !=
		        PsimagLite::String::npos) return;

		for (SizeType i = 0; i < files_.size(); ++i)
			unlink(files_[i].c_str());
	}

	void write(const TargetingType& psi,
	           const VectorSizeType& v,
	           int lastSign,
	           bool isObserveCode,
	           typename IoType::Out& ioOutCurrent) const
	{
		if (checkpoint_.parameters().recoverySave == "0")
			return;

		PsimagLite::String prefix("Recovery");
		prefix += ttos(counter_++);
		PsimagLite::String savedName(prefix + checkpoint_.parameters().filename);
		files_.push_back(savedName);
		ioOutCurrent.flush();

		std::ifstream source(ioOutCurrent.filename().c_str(), std::ios::binary);
		std::ofstream dest(savedName.c_str(), std::ios::binary);
		dest << source.rdbuf();
	    source.close();
	    dest.close();

		typename IoType::Out ioOut(savedName, IoType::ACC_RDW);

		// taken from end of finiteDmrgLoops
		checkpoint_.write(pS_, pE_, ioOut);
		ioOut.createGroup("FinalPsi");
		psi.write(v, ioOut, "FinalPsi");
		ioOut.write(lastSign, "LastLoopSign");

		// wft dtor
		wft_.write(ioOut);

		ioOut.close();
		// checkpoint stacks
		checkpoint_.checkpointStacks(savedName);
	}

private:

	PsimagLite::ProgressIndicator progress_;
	const CheckpointType& checkpoint_;
	const WaveFunctionTransfType& wft_;
	const BasisWithOperatorsType& pS_;
	const BasisWithOperatorsType& pE_;
	mutable SizeType counter_;
	mutable VectorStringType files_;
};     //class Recovery

} // namespace Dmrg
/*@}*/
#endif

