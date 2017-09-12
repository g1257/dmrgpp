/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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
	      flag_m_(false)
	{}

	~Recovery()
	{
		if (checkpoint_.parameters().options.find("recoveryNoDelete") !=
		        PsimagLite::String::npos) return;

		for (SizeType i = 0; i < files_.size(); ++i)
			unlink(files_[i].c_str());
	}

	void save(const TargetingType& psi,
	          VectorSizeType vsites,
	          int lastSign,
	          bool isObserveCode) const
	{
		if (checkpoint_.parameters().recoverySave == "0")
			return;

		PsimagLite::String prefix("Recovery");
		prefix += (flag_m_) ? "1" : "0";
		PsimagLite::String rootName(prefix + checkpoint_.parameters().filename);
		typename IoType::Out ioOut(rootName);
		ioOut<<checkpoint_.parameters();
		checkpoint_.save(pS_,pE_,ioOut);
		psi.save(vsites,ioOut);
		PsimagLite::OstringStream msg;
		msg<<"#LastLoopSign="<<lastSign<<"\n";
		ioOut<<msg.str();
		files_.push_back(rootName);

		saveStacksForRecovery(rootName,isObserveCode);

		wft_.save(rootName);
		wft_.appendFileList(files_,rootName);
		flag_m_ = !flag_m_;
	}

private:

	void saveStacksForRecovery(PsimagLite::String rootWriteFile,
	                           bool isObserveCode) const
	{
		PsimagLite::OstringStream msg;
		msg<<"Writing sys. and env. stacks to disk (for recovery)...";
		progress_.printline(msg,std::cout);
		PsimagLite::String sysWriteFile = utils::pathPrepend(checkpoint_.SYSTEM_STACK_STRING,
		                                                     rootWriteFile);
		PsimagLite::String envWriteFile = utils::pathPrepend(checkpoint_.ENVIRON_STACK_STRING,
		                                                     rootWriteFile);
		PsimagLite::String sysReadFile = "/dev/null";
		PsimagLite::String envReadFile = "/dev/null";

		{
			MemoryStackType systemStackCopy(checkpoint_.memoryStack(SYSTEM));
			DiskStackType systemDiskTemp(sysReadFile,sysWriteFile,false,isObserveCode);
			files_.push_back(sysWriteFile);
			CheckpointType::loadStack(systemDiskTemp,systemStackCopy);
		}

		{
			MemoryStackType envStackCopy(checkpoint_.memoryStack(ENVIRON));
			DiskStackType envDiskTemp(envReadFile,envWriteFile,false,isObserveCode);
			files_.push_back(envWriteFile);
			CheckpointType::loadStack(envDiskTemp,envStackCopy);
		}
	}

	PsimagLite::ProgressIndicator progress_;
	const CheckpointType& checkpoint_;
	const WaveFunctionTransfType& wft_;
	const BasisWithOperatorsType& pS_;
	const BasisWithOperatorsType& pE_;
	mutable bool flag_m_;
	mutable VectorStringType files_;
};     //class Recovery

} // namespace Dmrg
/*@}*/
#endif

