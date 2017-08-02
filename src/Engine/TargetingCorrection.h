/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file TargetingCorrection.h
 *
 *  corrects the finite-size algorithm
 *  following PRB 72, 180403(R) (2005)
 *
 */

#ifndef CORRECTION_TARGETTING_H
#define CORRECTION_TARGETTING_H
#include <iostream>
#include "TargetParamsCorrection.h"
#include "TargetingBase.h"
#include <stdexcept>

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingCorrection : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

public:

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType,VectorWithOffsetType_> BaseType;
	typedef typename BaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type TargetVectorType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef TargetParamsCorrection<ModelType> TargetParamsType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename BaseType::InputSimpleOutType InputSimpleOutType;

	enum {DISABLED,ENABLED};

	TargetingCorrection(const LeftRightSuperType& lrs,
	                    const ModelType& model,
	                    const WaveFunctionTransfType& wft,
	                    const SizeType&,
	                    InputValidatorType& io)
	    : BaseType(lrs,model,wft,0),
	      tstStruct_(io,model),
	      progress_("TargetingCorrection")
	{
		this->common().init(&tstStruct_,1);
	}

	RealType normSquared(SizeType i) const
	{
		return PsimagLite::real(this->common().targetVectors()[i]*
		                          this->common().targetVectors()[i]);
	}

	RealType weight(SizeType) const
	{
		assert(this->common().noStageIs(DISABLED));
		return tstStruct_.correctionA();
	}

	RealType gsWeight() const
	{
		return 1;
	}

	void evolve(RealType,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType&,
	            SizeType)
	{
		this->common().cocoon(block1,direction);

		if (direction == ProgramGlobals::INFINITE) return;
		this->common().setAllStagesTo(ENABLED);
		this->common().computeCorrection(direction,block1);
	}

	void save(const typename PsimagLite::Vector<SizeType>::Type& block,
	          PsimagLite::IoSimple::Out& io) const
	{
		PsimagLite::OstringStream msg;
		msg<<"Saving state...";
		progress_.printline(msg,std::cout);

		if (block.size()!=1)
			throw PsimagLite::RuntimeError("GST only supports blocks of size 1\n");

		PsimagLite::String s = "#TCENTRALSITE=" + ttos(block[0]);
		io.printline(s);
		this->common().psi().save(io,"PSI");
	}

	void load(const PsimagLite::String& f)
	{
		this->common().template load<int>(f,0);
	}

	void print(InputSimpleOutType& ioOut) const
	{
		ioOut.print("TARGETSTRUCT",tstStruct_);
		PsimagLite::OstringStream msg;
		msg<<"PSI\n";
		msg<<"CorrectionWeightGroundState=1\n";
		ioOut.print(msg.str());
	}

private:

	TargetParamsType tstStruct_;
	PsimagLite::ProgressIndicator progress_;
};     //class TargetingCorrection

template<typename LanczosSolverType, typename VectorWithOffsetType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingCorrection<LanczosSolverType, VectorWithOffsetType>& tst)
{
	tst.print(os);
	return os;
}
} // namespace Dmrg
/*@}*/
#endif

