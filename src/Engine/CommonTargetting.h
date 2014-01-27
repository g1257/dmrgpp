/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file CommonTargetting.h
 *
 * Functionality used by many targetting classes
 *
 */

#ifndef COMMON_TARGETTING_H
#define COMMON_TARGETTING_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "DynamicSerializer.h"
#include "DynamicDmrgParams.h"
#include "VectorWithOffsets.h"
#include "ContinuedFraction.h"
#include <cassert>
#include "ApplyOperatorExpression.h"
#include "TargetHelper.h"

namespace Dmrg {

template<typename TargetHelperType,
         typename VectorWithOffsetType,
         typename LanczosSolverType>
class CommonTargetting  {

public:

	typedef typename TargetHelperType::RealType RealType;
	typedef typename TargetHelperType::ModelType ModelType;
	typedef typename TargetHelperType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef DynamicSerializer<VectorWithOffsetType,PostProcType> DynamicSerializerType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::BasisDataType BasisDataType;
	typedef typename BasisType::BlockType BlockType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef ApplyOperatorExpression<TargetHelperType,
	                                VectorWithOffsetType,
	                                LanczosSolverType> ApplyOperatorExpressionType;
	typedef typename ApplyOperatorExpressionType::VectorSizeType VectorSizeType;
	typedef typename ApplyOperatorExpressionType::ApplyOperatorType ApplyOperatorType;
	typedef typename ApplyOperatorType::BorderEnum BorderEnumType;
	typedef typename TargetHelperType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename TargetHelperType::TargettingParamsType TargettingParamsType;
	typedef typename ApplyOperatorExpressionType::VectorVectorWithOffsetType VectorVectorWithOffsetType;
	typedef typename ApplyOperatorExpressionType::VectorRealType VectorRealType;
	typedef typename ApplyOperatorExpressionType::PairType PairType;

	enum {DISABLED=ApplyOperatorExpressionType::DISABLED,
		  OPERATOR=ApplyOperatorExpressionType::OPERATOR,
		  WFT_NOADVANCE=ApplyOperatorExpressionType::WFT_NOADVANCE};

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
	      EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
	      INFINITE=WaveFunctionTransfType::INFINITE};

	CommonTargetting(const LeftRightSuperType& lrs,
	                 const ModelType& model,
                     const TargettingParamsType& tstStruct,
                     const WaveFunctionTransfType& wft,
	                 SizeType targets)
	    : progress_("CommonTargetting"),
	      targetHelper_(lrs,model,tstStruct,wft),
	      applyOpExpression_(targetHelper_,targets)
	{}

	SizeType getPhi(VectorWithOffsetType& phiNew,
	                RealType Eg,
	                SizeType direction,
	                SizeType site,
	                SizeType loopNumber)
	{
		return applyOpExpression_.getPhi(phiNew,Eg,direction,site,loopNumber);
	}

	VectorWithOffsetType& psi() // <--- FIXME
	{
		return applyOpExpression_.psi();
	}

	const VectorWithOffsetType& psi() const
	{
		return applyOpExpression_.psi();
	}

	RealType normSquared(SizeType i) const
	{
		const VectorWithOffsetType& v = applyOpExpression_.targetVectors()[i];
		// call to mult will conjugate one of the vector
		return std::real(multiply(v,v));
	}

	template<typename IoOutputType>
	void save(const typename PsimagLite::Vector<SizeType>::Type& block,
	          IoOutputType& io,
	          const PostProcType& cf,
	          const typename PsimagLite::Vector<VectorWithOffsetType>::Type& targetVectors) const
	{
		DynamicSerializerType dynS(cf,block[0],targetVectors);
		dynS.save(io);
	}

	template<typename SomeSerializerType>
	void loadTargetVectors(SomeSerializerType& serializer)
	{
		applyOpExpression_.loadTargetVectors(serializer);
	}

	void checkOrder(SizeType i,
	                const typename PsimagLite::Vector<SizeType>::Type& stage) const
	{
		return applyOpExpression_.checkOrder(i,stage);
	}

	bool allStages(SizeType x) const
	{
		return applyOpExpression_.allStages(x);
	}

	bool noStageIs(SizeType x) const
	{
		return applyOpExpression_.noStageIs(x);
	}

	void initialGuess(VectorWithOffsetType& v,
	                  const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		PsimagLite::Vector<SizeType>::Type nk;
		setNk(nk,block);
		setInitialVector(v,applyOpExpression_.psi(),nk);
	}

	void noCocoon(const PsimagLite::String& msg) const
	{
		std::cout<<"-------------&*&*&* In-situ measurements start\n";
		std::cout<<"----- NO IN-SITU MEAS. POSSIBLE, reason="<<msg<<"\n";
		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	// in situ computation:
	void cocoon(SizeType direction,
	            SizeType site,
	            const VectorWithOffsetType& v,
	            const PsimagLite::String& label) const
	{
		std::cout<<"-------------&*&*&* In-situ measurements start\n";

		cocoon_(direction,site,v,label,ApplyOperatorType::BORDER_NO);

		int site2 = findBorderSiteFrom(site,direction);

		if (site2 >= 0) {
			cocoon_(direction,site2,v,label,ApplyOperatorType::BORDER_YES);
		}

		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	void computeCorrection(SizeType direction,
	                       const BlockType& block1)
	{
		const VectorWithOffsetType& psi = applyOpExpression_.psi();
		VectorWithOffsetType& v = applyOpExpression_.targetVectors(0);

		// operators in the one-site basis:
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		SparseMatrixType hmatrix;
		BasisDataType q;

		RealType time = 0;
		targetHelper_.model().setNaturalBasis(creationMatrix,hmatrix,q,block1,time);

		FermionSign fs(targetHelper_.lrs().left(),q.electrons);
		for (SizeType j=0;j<creationMatrix.size();j++) {
			VectorWithOffsetType phiTemp;
			applyOpExpression_.applyOpLocal()(phiTemp,psi,creationMatrix[j],
			              fs,direction,ApplyOperatorType::BORDER_NO);
			if (j==0) v = phiTemp;
			else v += phiTemp;
		}
	}

	void setNk(typename PsimagLite::Vector<SizeType>::Type& nk,
	           const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		for (SizeType i=0;i<block.size();i++)
			nk.push_back(targetHelper_.model().hilbertSize(block[i]));
	}

	RealType setGsWeight(RealType defaultValue) const
	{
		if (!targetHelper_.model().params().gsWeight.first)
			return 	defaultValue;

		return targetHelper_.model().params().gsWeight.second;
	}

	int findFermionSignOfTheOperators() const
	{
		const VectorOperatorType& myoperator = targetHelper_.tstStruct().aOperators();
		int f = 0;

		for (SizeType i = 0; i < myoperator.size(); ++i) {

			RealType norma = norm2(myoperator[i].data);

			if (norma==0) continue;

			if (f == 0) {
				f = myoperator[i].fermionSign;
				continue;
			}

			if (f == myoperator[i].fermionSign) continue;

			PsimagLite::String str("CorrectionVectorTargetting: ");
			str += "inconsistent sign for operators\n";
			throw PsimagLite::RuntimeError(str);

		}

		return f;
	}

	void setAllStagesTo(SizeType x)
	{
		applyOpExpression_.setAllStagesTo(x);
	}

	const RealType& energy() const
	{
		return applyOpExpression_.energy();
	}

	const RealType& currentTime() const
	{
		return applyOpExpression_.currentTime();
	}

	const VectorSizeType& nonZeroQns() const
	{
		 return applyOpExpression_.nonZeroQns();
	}

	const VectorVectorWithOffsetType& targetVectors() const
	{
		return applyOpExpression_.targetVectors();
	}

	VectorWithOffsetType& targetVectors(SizeType i)
	{
		return applyOpExpression_.targetVectors(i);
	}

	void targetVectorsResize(SizeType x)
	{
		applyOpExpression_.targetVectorsResize(x);
	}

	void initTimeVectors(const VectorRealType& times)
	{
		applyOpExpression_.initTimeVectors(times);
	}

	void setTime(RealType t)
	{
		applyOpExpression_.setTime(t);
	}

	void calcTimeVectors(const PairType& startEnd,
	                             RealType Eg,
	                             const VectorWithOffsetType& phi,
	                             SizeType direction,
	                             bool allOperatorsApplied,
	                             const PsimagLite::Vector<SizeType>::Type& block)
	{
		return applyOpExpression_.calcTimeVectors(startEnd,
		                                        Eg,
			                                    phi,
			                                    direction,
			                                    allOperatorsApplied,
			                                    block);
	}

	void cocoonLegacy(SizeType direction,const BlockType& block) const
	{
		const VectorWithOffsetType& psi = applyOpExpression_.psi();
		const VectorWithOffsetType& tv0 = applyOpExpression_.targetVectors()[0];
		const ModelType& model = targetHelper_.model();

		SizeType site = block[0];
		PsimagLite::CrsMatrix<ComplexOrRealType> tmpC(model.naturalOperator("nup",site,0));
		int fermionSign1 = 1;
		const std::pair<SizeType,SizeType> jm1(0,0);
		RealType angularFactor1 = 1.0;
		typename OperatorType::Su2RelatedType su2Related1;
		OperatorType nup(tmpC,fermionSign1,jm1,angularFactor1,su2Related1);

		nup.data = tmpC;
		nup.fermionSign = 1;

		test(psi,psi,direction,"<PSI|nup|PSI>",site,nup,ApplyOperatorType::BORDER_NO);
		PsimagLite::String s = "<P0|nup|P0>";
		test(tv0,tv0,direction,s,site,nup,ApplyOperatorType::BORDER_NO);

		PsimagLite::CrsMatrix<ComplexOrRealType> tmpC2(model.naturalOperator("ndown",site,0));
		OperatorType ndown(tmpC2,fermionSign1,jm1,angularFactor1,su2Related1);
		test(psi,psi,direction,"<PSI|ndown|PSI>",site,ndown,ApplyOperatorType::BORDER_NO);
		s = "<P0|ndown|P0>";
		test(tv0,tv0,direction,s,site,ndown,ApplyOperatorType::BORDER_NO);

		PsimagLite::CrsMatrix<ComplexOrRealType> tmpC3 = (nup.data * ndown.data);
		OperatorType doubleOcc(tmpC3,fermionSign1,jm1,angularFactor1,su2Related1);
		test(psi,psi,direction,"<PSI|doubleOcc|PSI>",site,doubleOcc,ApplyOperatorType::BORDER_NO);
		s = "<P0|doubleOcc|P0>";
		test(tv0,tv0,direction,s,site,doubleOcc,ApplyOperatorType::BORDER_NO);
	}

private:

	void setInitialVector(VectorWithOffsetType& v1,
	                      const VectorWithOffsetType& v2,
	                      const VectorSizeType& nk) const
	{
		const WaveFunctionTransfType& wft = targetHelper_.wft();
		bool noguess = (targetHelper_.model().params().options.find("targetnoguess") !=
		        PsimagLite::String::npos);

		if (noguess)
			wft.createRandomVector(v1);
		else
			wft.setInitialVector(v1,v2,targetHelper_.lrs(),nk);
	}

	int findBorderSiteFrom(SizeType site, SizeType direction) const
	{
		if (site == 1 && direction == EXPAND_ENVIRON) return 0;

		SizeType n = targetHelper_.model().geometry().numberOfSites();
		if (site == n - 2 && direction == EXPAND_SYSTEM) return n - 1;

		return -1;
	}

	void cocoon_(SizeType direction,SizeType site,
	             const VectorWithOffsetType& v,
	             const PsimagLite::String& label,
	             BorderEnumType border) const
	{
		VectorStringType vecStr = getOperatorLabels();

		for (SizeType i=0;i<vecStr.size();i++) {
			const PsimagLite::String& opLabel = vecStr[i];
			OperatorType nup = getOperatorForTest(opLabel,site);

			PsimagLite::String tmpStr = "<"+ label + "|" + opLabel + "|" + label + ">";
			test(v,v,direction,tmpStr,site,nup,border);
		}
	}

	VectorStringType getOperatorLabels() const
	{
		VectorStringType vecStr;
		PsimagLite::tokenizer(targetHelper_.model().params().insitu,vecStr,",");
		return vecStr;
	}

	OperatorType getOperatorForTest(const PsimagLite::String& opLabel,
	                                SizeType site) const
	{
		int fermionSign1 = 1;
		const std::pair<SizeType,SizeType> jm1(0,0);
		RealType angularFactor1 = 1.0;
		typename OperatorType::Su2RelatedType su2Related1;

		OperatorType nup;
		try {
			nup = findOperator(opLabel);
		} catch (std::exception& e) {
			if (opLabel[0] == ':') {
				std::cerr<<e.what();
				throw e;
			}

			PsimagLite::CrsMatrix<ComplexOrRealType> tmpC(targetHelper_.model().naturalOperator(opLabel,site,0));
			nup = OperatorType(tmpC,fermionSign1,jm1,angularFactor1,su2Related1);
		}

		return nup;
	}

	void test(const VectorWithOffsetType& src1,
	          const VectorWithOffsetType& src2,
	          SizeType systemOrEnviron,
	          const PsimagLite::String& label,
	          SizeType site,
	          const OperatorType& A,
	          BorderEnumType border) const
	{
		typename PsimagLite::Vector<SizeType>::Type electrons;
		targetHelper_.model().findElectronsOfOneSite(electrons,site);
		FermionSign fs(targetHelper_.lrs().left(),electrons);
		VectorWithOffsetType dest;
		applyOpExpression_.applyOpLocal()(dest,src1,A,fs,systemOrEnviron,border);

		ComplexOrRealType sum = 0.0;
		for (SizeType ii=0;ii<dest.sectors();ii++) {
			SizeType i = dest.sector(ii);
			SizeType offset1 = dest.offset(i);
			for (SizeType jj=0;jj<src2.sectors();jj++) {
				SizeType j = src2.sector(jj);
				SizeType offset2 = src2.offset(j);
				if (i!=j) continue;
				for (SizeType k=0;k<dest.effectiveSize(i);k++)
					sum+= dest[k+offset1] * std::conj(src2[k+offset2]);
			}
		}
		std::cout<<site<<" "<<sum<<" "<<currentTime();
		std::cout<<" "<<label<<" "<<(src1*src2)<<"\n";
	}

	OperatorType findOperator(const PsimagLite::String& name) const
	{
		if (name.length()<2 || name[0]!=':') {
			PsimagLite::String str("ObserverInterpreter: syntax error for ");
			str += name + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		PsimagLite::String label = name.substr(1,name.length()-1);

		PsimagLite::IoSimple::In io(label);

		CookedOperator<ModelType> cookedOperator(targetHelper_.model());

		return OperatorType(io,cookedOperator,OperatorType::MUST_BE_NONZERO);
	}

	PsimagLite::ProgressIndicator progress_;
	TargetHelperType targetHelper_;
	ApplyOperatorExpressionType applyOpExpression_;
}; // class CommonTargetting

template<typename TargetHelperType,
         typename VectorWithOffsetType,
         typename LanczosSolverType>
std::ostream& operator<<(std::ostream& os,
                         const CommonTargetting<TargetHelperType,
                                                VectorWithOffsetType,
                                                LanczosSolverType>& tst)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // COMMON_TARGETTING_H

