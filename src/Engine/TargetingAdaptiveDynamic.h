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

/*! \file TargetingAdaptiveDynamic.h
 *
 * Implements the targeting required by
 * arXiv:1012.5543v1
 *
 */

#ifndef TARGETING_ADAPTIVE_DYN_H
#define TARGETING_ADAPTIVE_DYN_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TimeSerializer.h"
#include "TargetParamsAdaptiveDynamic.h"
#include "VectorWithOffsets.h"
#include "ParametersForSolver.h"
#include "TargetingBase.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingAdaptiveDynamic : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {
public:

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType, VectorWithOffsetType_> BaseType;
	typedef typename BaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef TargetParamsAdaptiveDynamic<ModelType> TargetParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename BaseType::IoType IoType;

	enum {DISABLED=BaseType::DISABLED,
		  OPERATOR=BaseType::OPERATOR,
		  WFT_NOADVANCE=BaseType::WFT_NOADVANCE};

	TargetingAdaptiveDynamic(const LeftRightSuperType& lrs,
	                         const ModelType& model,
	                         const WaveFunctionTransfType& wft,
	                         const SizeType&,
	                         InputValidatorType& io)

	    : BaseType(lrs,model,wft,0),
	      tstStruct_(io,model),
	      wft_(wft),
	      lastLanczosVector_(0),
	      dynCounter_(0),
	      progress_("TargetingAdaptiveDynamic"),
	      gsWeight_(1.0),
	      done_(false),
	      weightForContinuedFraction_(0),
	      paramsForSolver_(io,"DynamicDmrg")
	{
		this->common().init(&tstStruct_,2);
		if (!wft.isEnabled()) throw PsimagLite::RuntimeError(" DynamicTargeting "
		                                                     "needs an enabled wft\n");
	}

	RealType weight(SizeType i) const
	{
		assert(this->common().allStages(DISABLED));
		return weight_[i];
	}

	RealType gsWeight() const
	{
		return gsWeight_;
	}

	SizeType size() const
	{
		if (!this->common().allStages(WFT_NOADVANCE)) return 0;
		return lastLanczosVector_;
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
		SizeType numberOfSites = this->lrs().super().block().size();
		if (site>1 && site<numberOfSites-2) return;
		//corner case
		SizeType x = (site==1) ? 0 : numberOfSites-1;
		evolve(Eg,direction,x,loopNumber);
	}

	void print(typename IoType::Out& ioOut) const
	{
		ioOut.print("TARGETSTRUCT",tstStruct_);
		PsimagLite::OstringStream msg;
		msg<<"PSI\n";
		msg<<(*this);
		ioOut.print(msg.str());
	}

	void save(const typename PsimagLite::Vector<SizeType>::Type& block,
	          typename IoType::Out& io) const
	{
		assert(block.size()==1);
		SizeType type = tstStruct_.type();
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;

		if (ab_.size()<2) return;

		paramsForSolver_.Eg = Eg_;
		paramsForSolver_.weight = s2*weightForContinuedFraction_;
		paramsForSolver_.isign = s;
		PostProcType cf(ab_, paramsForSolver_);
		PsimagLite::String str = "#TCENTRALSITE=" + ttos(block[0]);
		io.printline(str);
		this->common().save(block,io,cf,this->common().targetVectors());

		this->common().psi().save(io,"PSI");
	}

	void load(const PsimagLite::String& f)
	{
		this->common().template load<TimeSerializerType>(f);
		lastLanczosVector_ = this->common().targetVectors().size()-1;
		dynCounter_ = 13;
	}

private:

	void evolve(RealType Eg,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType site,
	            SizeType loopNumber)
	{
		Eg_ = Eg;
		VectorWithOffsetType phiNew;
		if (ab_.size()==0)
			this->common().getPhi(phiNew,Eg,direction,site,loopNumber);

		if (!this->common().allStages(WFT_NOADVANCE)) {
			this->common().targetVectors(0) = phiNew;
			return;
		}

		SizeType numberOfSites = this->lrs().super().block().size();
		if (site>0 && site<numberOfSites-1) wftAllDynVectors(site);

		if (!done_) calcDynVectors(site,phiNew);
	}

	void wftAllDynVectors(SizeType site)
	{
		for (SizeType i=0;i<=lastLanczosVector_;i++)
			if (i<2) wftOneDynVector(i,site);
	}

	void wftOneDynVector(SizeType i,SizeType site)
	{
		typename PsimagLite::Vector<SizeType>::Type nk(1,this->model().hilbertSize(site));

		VectorWithOffsetType result;
		result.populateSectors(this->lrs().super());

		// OK, now that we got the partition number right, let's wft:

		// FIXME generalize for su(2)
		wft_.setInitialVector(result,this->common().targetVectors()[i],this->lrs(),nk);
		result.collapseSectors();
		this->common().targetVectors(i) = result;
	}

	void calcDynVectors(SizeType site,const VectorWithOffsetType& phiNew)
	{
		for (SizeType i=0;i<this->common().targetVectors()[0].sectors();i++) {
			VectorType sv;
			SizeType i0 = this->common().targetVectors()[0].sector(i);
			this->common().targetVectors()[0].extract(sv,i0);
			SizeType p = this->lrs().super().findPartitionNumber(this->common().
			                                                     targetVectors()[0].offset(i0));
			if (i==0) {
				if (lastLanczosVector_==0)
					this->common().targetVectors(1) = this->common().targetVectors()[0];
			}
			setLanczosVectors(i0,sv,p,site,phiNew);
		}
		setWeights();
		if (lastLanczosVector_==1 && fabs(weightForContinuedFraction_)<1e-6)
			weightForContinuedFraction_ = PsimagLite::real(this->common().targetVectors()
			                                               [0]*this->common().targetVectors()[0]);
	}

	void setLanczosVectors(SizeType i0,
	                       const VectorType& sv,
	                       SizeType p,
	                       SizeType,
	                       const VectorWithOffsetType& phiNew)
	{
		SizeType threadId = 0;
		RealType fakeTime = 0;
		typename ModelType::ModelHelperType modelHelper(p,this->lrs(),fakeTime,threadId);
		typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
		LanczosMatrixType h(&this->model(),&modelHelper);
		LanczosSolverType lanczosSolver(h,paramsForSolver_);

		RealType a=0,b=0;
		VectorType x(sv.size(),0.0);
		VectorType y = sv;
		if (lastLanczosVector_>0) {
			this->common().targetVectors()[1].extract(y,i0);
			this->common().targetVectors()[0].extract(x,i0);
		}

		if (lastLanczosVector_==0) normalize(y);

		lanczosSolver.oneStepDec(x,y,a,b,lastLanczosVector_);

		if (lastLanczosVector_<2) lastLanczosVector_++;

		//f0 is wft'd, do nothing
		//f1 is wft'd, do nothing
		RealType norm1 = PsimagLite::norm(x);
		if (norm1<1e-6) {
			lanczosSolver.push(ab_,a,b);
			h.matrixVectorProduct(x,y);
			a = PsimagLite::real(x*y);
			lanczosSolver.push(ab_,a,b);
			done_=true;
			return;
		}
		dynCounter_++;
		if (lastLanczosVector_>1 && (dynCounter_%tstStruct_.advanceEach()) != 0)
			return;
		this->common().targetVectors(0).setDataInSector(x,i0);
		this->common().targetVectors(1).setDataInSector(y,i0);
		if ((dynCounter_%tstStruct_.advanceEach()) != 0) return;
		if (ab_.size()>0) {
			lanczosSolver.push(ab_,a,b);
			return;
		}
		// first push:
		VectorType xx(sv.size(),0.0);
		VectorType yy;
		phiNew.extract(yy,i0);
		normalize(yy);
		RealType a1=0,b1=0;
		lanczosSolver.oneStepDec(xx,yy,a1,b1,lastLanczosVector_);
		lanczosSolver.push(ab_,a1,b1);
		if (tstStruct_.advanceEach()<=1) return;
		lanczosSolver.push(ab_,a,b);
	}

	void setWeights()
	{
		RealType sum  = 0;
		weight_.resize(this->common().targetVectors().size());
		for (SizeType r=0;r<weight_.size();r++) {
			weight_[r] = (r>lastLanczosVector_) ? 0 : 1;
			sum += weight_[r];
		}

		gsWeight_ = tstStruct_.gsWeight();
		for (SizeType r=0;r<weight_.size();r++) weight_[r] *= (1-gsWeight_)/sum;
	}

	RealType dynWeightOf(VectorType& v,const VectorType& w) const
	{
		RealType sum = 0;
		for (SizeType i=0;i<v.size();i++) {
			if (i>=w.size()) continue;
			RealType tmp = PsimagLite::real(v[i]*w[i]);
			sum += tmp*tmp;
		}
		return sum;
	}

	void normalize(VectorType& v) const
	{
		RealType x = PsimagLite::norm(v);
		if (fabs(x)<1e-6) return;
		v /= x;
	}

	TargetParamsType tstStruct_;
	const WaveFunctionTransfType& wft_;
	SizeType lastLanczosVector_;
	SizeType dynCounter_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	bool done_;
	RealType Eg_;
	RealType weightForContinuedFraction_;
	TridiagonalMatrixType ab_;
	mutable typename LanczosSolverType::ParametersSolverType paramsForSolver_;
}; // class DynamicTargeting

template<typename LanczosSolverType, typename VectorWithOffsetType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingAdaptiveDynamic<LanczosSolverType,VectorWithOffsetType>&)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // TARGETING_ADAPTIVE_DYN_H

