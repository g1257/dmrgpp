/*
Copyright (c) 2009-2016-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file TargetingRixsStatic.h
 *
 * Implements the targeting required by
 * RIXS Static, restarts from Correction Vector run
 * tv[0] = A^\dagger_{site}|gs>
 * tv[1] = Imaginary of (w*-Htilde+i\eta)^{-1} A^\dagger_{site}|gs>
 * tv[2] = Real of (w*-Htilde+i\eta)^{-1} A^\dagger_{site}|gs>
 * tv[3] = A^\dagger_{sitep}|gs>
 * tv[4] = Imaginary of (w*-Htildep+i\eta)^{-1} A^\dagger_{sitep}|gs>
 * tv[5] = Real of (w*-Htildep+i\eta)^{-1} A^\dagger_{sitep}|gs>
 *
 */

#ifndef TARGETING_RIXS_STATIC_H
#define TARGETING_RIXS_STATIC_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TargetParamsCorrectionVector.h"
#include "VectorWithOffsets.h"
#include "TargetingBase.h"
#include "ParametersForSolver.h"
#include "ParallelTriDiag.h"
#include "FreqEnum.h"
#include "CorrectionVectorSkeleton.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingRixsStatic : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType,VectorWithOffsetType_> BaseType;

public:

	typedef typename BaseType::TargetingCommonType TargetingCommonType;
	typedef typename BaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef TargetParamsCorrectionVector<ModelType> TargetParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename TargetingCommonType::TimeSerializerType TimeSerializerType;
	typedef typename LanczosSolverType::MatrixType LanczosMatrixType;
	typedef CorrectionVectorFunction<LanczosMatrixType,TargetParamsType>
	CorrectionVectorFunctionType;
	typedef ParallelTriDiag<ModelType,LanczosSolverType,VectorWithOffsetType>
	ParallelTriDiagType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType VectorMatrixFieldType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef CorrectionVectorSkeleton<LanczosSolverType,
	VectorWithOffsetType,
	BaseType,
	TargetParamsType> CorrectionVectorSkeletonType;
	typedef typename BasisType::QnType QnType;

	TargetingRixsStatic(const LeftRightSuperType& lrs,
	                    const ModelType& model,
	                    const WaveFunctionTransfType& wft,
	                    const QnType&,
	                    InputValidatorType& ioIn)
	    : BaseType(lrs,model,wft,1),
	      tstStruct_(ioIn,model),
	      ioIn_(ioIn),
	      progress_("TargetingRixsStatic"),
	      gsWeight_(1.0),
	      paramsForSolver_(ioIn,"DynamicDmrg"),
	      skeleton_(ioIn_,tstStruct_,model,lrs,this->common().aoe().energy()),
	      applied_(false),
	      appliedFirst_(false)
	{
		if (!wft.isEnabled())
			err("TargetingRixsStatic needs wft\n");
	}

	SizeType sites() const { return tstStruct_.sites(); }

	SizeType targets() const { return 6; }

	RealType weight(SizeType i) const
	{
		return weight_[i];
	}

	RealType gsWeight() const
	{
		return gsWeight_;
	}

	SizeType size() const
	{
		if (!applied_ && appliedFirst_) {
			return 4;
		}
		return (applied_) ? 6 : 3;
	}

	void evolve(RealType Eg,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType& block2,
	            SizeType loopNumber)
	{
		if (block1.size()!=1 || block2.size()!=1) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "evolve only blocks of one site supported\n";
			throw PsimagLite::RuntimeError(str.c_str());
		}

		SizeType site = block1[0];
		evolve(Eg,direction,site,loopNumber);
		skeleton_.printNormsAndWeights(this->common(), weight_, gsWeight_);
	}

	void write(const VectorSizeType& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, block, prefix);
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix)
	{
		this->common().read(io, prefix);

		TimeSerializerType ts(io, prefix);

		SizeType n = ts.numberOfVectors();
		if (n != 4)
			err("TargetingRixsStatic: number of TVs must be 4\n");

		for (SizeType site = 0; site < 3; ++site) {
			this->common().aoe().targetVectors(site) = ts.vector(site + 1);
		}
	}

private:

	// tv[0] = A_site |gs>
	// tv[1] = imaginary cv for tv[0]
	// tv[2] = real      cv for tv[0]
	// tv[3] = A_sitep |gs>
	// tv[4] = imaginary cv for tv[3]
	// tv[5] = real      cv for tv[3]

	void evolve(RealType,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType site,
	            SizeType loopNumber)
	{
		if (direction == ProgramGlobals::INFINITE) return;

		// see if operator at sitep has been applied and result put into targetVectors[3]
		// if no apply operator at site and add into targetVectors[3]
		// also wft everything

		this->common().aoe().wftSome(site, 0, this->common().aoe().targetVectors().size());

		SizeType max = tstStruct_.sites();

		if (max>2)
			err("You cannot apply more than 2 operators (only SUM is allowed)\n");

		if (!applied_) {
			if (max==1) {
				if (site == tstStruct_.sites(0)) {
					VectorWithOffsetType tmpV1;
					SizeType indexOfOperator = 0;
					this->common().aoe().applyOneOperator(loopNumber,
					                                indexOfOperator,
					                                site,
					                                tmpV1,
					                                this->common().aoe().psi(),
					                                direction,
					                                tstStruct_);
					if (tmpV1.size() > 0) {
						this->common().aoe().targetVectors(3) = tmpV1;
						applied_ = true;
						PsimagLite::OstringStream msg;
						msg<<"Applied operator";
						progress_.printline(msg, std::cout);
					}
				}
			}
			if (max==2) {
				if (site == tstStruct_.sites(0)) {
					VectorWithOffsetType tmpV1;
					SizeType indexOfOperator = 0;
					this->common().aoe().applyOneOperator(loopNumber,
					                                indexOfOperator,
					                                site,
					                                tmpV1,
					                                this->common().aoe().psi(),
					                                direction,
					                                tstStruct_);
					if (tmpV1.size() > 0) {
						this->common().aoe().targetVectors(3) = tmpV1;
						applied_ = false;
						appliedFirst_ = true;
						PsimagLite::OstringStream msg;
						msg<<"Applied first operator";
						progress_.printline(msg, std::cout);
					}
				}
				if (site == tstStruct_.sites(1)) {
					VectorWithOffsetType tmpV2;
					SizeType indexOfOperator = 1;
					this->common().aoe().applyOneOperator(loopNumber,
					                                indexOfOperator,
					                                site,
					                                tmpV2,
					                                this->common().aoe().psi(),
					                                direction,
					                                tstStruct_);
					if (tmpV2.size() > 0) {
						this->common().aoe().targetVectors(3) += tmpV2;
						applied_ = true;
						PsimagLite::OstringStream msg;
						msg<<"Applied second operator";
						progress_.printline(msg, std::cout);
					}
				}
			}
		}

		doCorrectionVector();

		bool doBorderIfBorder = true;
		VectorSizeType block(1, site);
		this->common().cocoon(block, direction, doBorderIfBorder);
	}

	void doCorrectionVector()
	{
		if (!applied_ && appliedFirst_) {
			setWeights(4);
			return;
		}

		if (!applied_) {
			setWeights(3);
			return;
		}

		skeleton_.calcDynVectors(this->common().aoe().targetVectors(3),
		                         this->common().aoe().targetVectors(4),
		                         this->common().aoe().targetVectors(5));
		//		this->common().aoe().targetVectors(4) = this->common().aoe().targetVectors(1);
		//		this->common().aoe().targetVectors(5) = this->common().aoe().targetVectors(2);

		RealType n4 = PsimagLite::real(this->common().aoe().targetVectors(4)*
		                               this->common().aoe().targetVectors(4));
		RealType n5 = PsimagLite::real(this->common().aoe().targetVectors(5)*
		                               this->common().aoe().targetVectors(5));
		std::cout<<"HERE============> n4="<<n4<<" n5="<<n5;
		std::cout<<" "<<this->common().aoe().energy()<<"\n";
		setWeights(6);
	}

	void setWeights(SizeType n)
	{
		gsWeight_ = tstStruct_.gsWeight();

		RealType sum  = n;
		weight_.resize(n, 1);

		for (SizeType r=0;r<weight_.size();r++) weight_[r] = (1.0 - gsWeight_)/sum;
	}

	TargetParamsType tstStruct_;
	InputValidatorType& ioIn_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	typename LanczosSolverType::ParametersSolverType paramsForSolver_;
	CorrectionVectorSkeletonType skeleton_;
	bool applied_;
	bool appliedFirst_;
}; // class TargetingRixsStatic
} // namespace
/*@}*/
#endif // TARGETING_RIXS_STATIC_H
