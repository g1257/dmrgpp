/*
Copyright (c) 2009, UT-Battelle, LLC
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

/*! \file TargetingGroundState.h
 *
 *  targets the ground state
 *
 */

#ifndef TARGETING_GS_H
#define TARGETING_GS_H
#include <iostream>
#include "String.h"
#include "TargetParamsGroundState.h"
#include "ApplyOperatorLocal.h"
#include <stdexcept>
#include "Tokenizer.h"
#include "TargetingBase.h"
#include "ParametersForSolver.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_,
         typename IoType_>
class TargetingGroundState : public TargetingBase<LanczosSolverTemplate,
                                                  MatrixVectorType_,
                                                  WaveFunctionTransfType_,
                                                  IoType_> {

public:

	typedef TargetingBase<LanczosSolverTemplate,
	                      MatrixVectorType_,
	                      WaveFunctionTransfType_,
                          IoType_> BaseType;
	typedef MatrixVectorType_ MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef LanczosSolverTemplate<ParametersForSolverType,MatrixVectorType,VectorType> LanczosSolverType;
	typedef typename BasisType::BlockType BlockType;
	typedef WaveFunctionTransfType_ WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef VectorType TargetVectorType;
	typedef TargetParamsGroundState<ModelType> TargettingParamsType;
	typedef typename ModelType::InputValidatorType InputValidatorType;

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
		  EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
		  INFINITE=WaveFunctionTransfType::INFINITE};

	TargetingGroundState(const LeftRightSuperType& lrs,
	                     const ModelType& model,
	                     const TargettingParamsType& tstStruct,
	                     const WaveFunctionTransfType& wft,
	                     const SizeType& quantumSector,
	                     InputValidatorType& io)
	    : BaseType(lrs,model,tstStruct,wft,0,0),
	      progress_("TargetingGroundState")
	{}

	RealType weight(SizeType i) const
	{
		throw PsimagLite::RuntimeError("GST: What are you doing here?\n");
		return 0;
	}

	RealType gsWeight() const
	{
		return 1;
	}

	SizeType size() const
	{
		return 0;
	}

	void evolve(RealType Eg,
	            SizeType direction,
	            const BlockType& block1,
	            const BlockType& block2,
	            SizeType loopNumber)
	{
		this->common().cocoon(block1,direction);
	}

	template<typename IoOutputType>
	void save(const typename PsimagLite::Vector<SizeType>::Type& block,IoOutputType& io) const
	{
		PsimagLite::OstringStream msg;
		msg<<"Saving state...";
		progress_.printline(msg,std::cout);

		assert(block.size()>0);
		PsimagLite::String s = "#TCENTRALSITE=" + ttos(block[0]);
		io.printline(s);
		this->common().psi().save(io,"PSI");
	}

	void load(const PsimagLite::String& f)
	{
		this->common().template load<int>(f,0);
	}

private:

	PsimagLite::ProgressIndicator progress_;

};     //class TargetingGroundState

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType,
         typename WaveFunctionTransfType,
         typename IoType_>
std::ostream& operator<<(std::ostream& os,
                         const TargetingGroundState<LanczosSolverTemplate,
                         MatrixVectorType,
                         WaveFunctionTransfType,IoType_>& tst)
{
	os<<"GSTWeightGroundState=1\n";
	return os;
}
} // namespace Dmrg
/*@}*/
#endif

