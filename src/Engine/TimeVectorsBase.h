/*
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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

/*! \file TimeVectorsBase.h
 *
 *
 */

#ifndef TIME_VECTORS_BASE
#define TIME_VECTORS_BASE
#include <iostream>
#include "Vector.h"
#include "ProgramGlobals.h"
#include "Wft/WftHelper.h"

namespace Dmrg {

template<typename TargetParamsType,
         typename ModelType,
         typename WaveFunctionTransfType,
         typename LanczosSolverType,
         typename VectorWithOffsetType>
class TimeVectorsBase  {

public:

	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename TargetParamsType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef WftHelper<ModelType, VectorWithOffsetType, WaveFunctionTransfType> WftHelperType;

	TimeVectorsBase(const ModelType& model,
	                const LeftRightSuperType& lrs,
	                const WaveFunctionTransfType& wft,
	                PsimagLite::String name)
	    : wftHelper_(model, lrs, wft), name_(name), time_(0.0), currentTimeStep_(0)
	{}

	struct ExtraData {

		ExtraData(ProgramGlobals::DirectionEnum dir_,
		          bool allOperatorsApplied_,
		          bool wftAndAdvanceIfNeeded_,
		          VectorSizeType block_,
		          bool isLastCall_)
		    : dir(dir_),
		      allOperatorsApplied(allOperatorsApplied_),
		      wftAndAdvanceIfNeeded(wftAndAdvanceIfNeeded_),
		      block(block_),
		      isLastCall(isLastCall_)
		{}

		ProgramGlobals::DirectionEnum dir;
		bool allOperatorsApplied;
		bool wftAndAdvanceIfNeeded;
		PsimagLite::Vector<SizeType>::Type block;
		bool isLastCall;
	};

	virtual void calcTimeVectors(const VectorSizeType&,
	                             RealType,
	                             const VectorWithOffsetType&,
	                             const ExtraData&)
	{
		err("calcTimeVectors: unimplemented in this base class\n");
	}

	virtual ~TimeVectorsBase() {}

	virtual void timeHasAdvanced() {}

	RealType time() const { return time_; }

	SizeType currentTimeStep() const { return currentTimeStep_; }

	void advanceCurrentTime(RealType tau) { time_ += tau; }

	void advanceCurrentTimeStep() { ++currentTimeStep_; }

	void setCurrentTimeStep(SizeType t) { currentTimeStep_ = t; }

	void setCurrentTime(RealType t) { time_ = t; }

	bool isBase() const { return (name_ == "base"); }

protected:

	const WftHelperType& wftHelper() const { return wftHelper_; }

private:

	WftHelperType wftHelper_;
	PsimagLite::String name_;
	RealType time_;
	SizeType currentTimeStep_;
}; //class TimeVectorsBase
} // namespace Dmrg
/*@}*/
#endif
