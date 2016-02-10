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

/*! \file TargetingCommon.h
 *
 * Functionality used by many targetting classes
 *
 */

#ifndef TARGETING_COMMON_H
#define TARGETING_COMMON_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "DynamicSerializer.h"
#include "TargetParamsDynamic.h"
#include "VectorWithOffsets.h"
#include "ContinuedFraction.h"
#include <cassert>
#include "ApplyOperatorExpression.h"
#include "IoSimple.h"
#include "Tokenizer.h"

namespace Dmrg {

template<typename TargetHelperType,
         typename VectorWithOffsetType,
         typename LanczosSolverType>
class TargetingCommon  {

	struct NaturalOpStruct {
		NaturalOpStruct(PsimagLite::String label_)
		    : fermionSign(1),dof(0),label(label_)
		{
			SizeType i = 0;
			for (; i < label_.length(); ++i) {
				if (label_[i] == '?') break;
			}

			if (i == label_.length()) return;
			SizeType j = i;
			label = label_.substr(0,j);
			for (; i < label_.length(); ++i) {
				if (label_[i] == '-') break;
			}

			dof = atoi(label_.substr(j+1,i).c_str());
			if (i == label_.length()) return;
			fermionSign = -1;
		}

		int fermionSign;
		SizeType dof;
		PsimagLite::String label;
	}; // struct NaturalOpStruct

public:

	typedef PsimagLite::IoSimple IoType;
	typedef typename IoType::In IoInputType;
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
	typedef typename BasisWithOperatorsType::SymmetryElectronsSzType
	SymmetryElectronsSzType;
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
	typedef typename TargetHelperType::TargetParamsType TargetParamsType;
	typedef typename ApplyOperatorExpressionType::VectorVectorWithOffsetType
	VectorVectorWithOffsetType;
	typedef typename ApplyOperatorExpressionType::VectorRealType VectorRealType;
	typedef typename ApplyOperatorExpressionType::PairType PairType;
	typedef typename ModelType::InputValidatorType InputValidatorType;

	static const SizeType SUM = TargetParamsType::SUM;

	enum {DISABLED=ApplyOperatorExpressionType::DISABLED,
	      OPERATOR=ApplyOperatorExpressionType::OPERATOR,
	      WFT_NOADVANCE=ApplyOperatorExpressionType::WFT_NOADVANCE};

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
	      EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
	      INFINITE=WaveFunctionTransfType::INFINITE};

	TargetingCommon(const LeftRightSuperType& lrs,
	                const ModelType& model,
	                const WaveFunctionTransfType& wft,
	                SizeType indexNoAdvance)
	    : progress_("TargetingCommon"),
	      targetHelper_(lrs,model,wft),
	      applyOpExpression_(targetHelper_,indexNoAdvance),
	      inSitu_(model.geometry().numberOfSites())
	{}

	void init(TargetParamsType* tstStruct, SizeType targets)
	{
		targetHelper_.setTargetStruct(tstStruct);
		applyOpExpression_.init();
		targetVectorsResize(targets);
	}

	SizeType getPhi(VectorWithOffsetType& phiNew,
	                RealType Eg,
	                SizeType direction,
	                SizeType site,
	                SizeType loopNumber)
	{
		return applyOpExpression_.getPhi(phiNew,Eg,direction,site,loopNumber);
	}

	const VectorWithOffsetType& psi() const
	{
		return applyOpExpression_.psi();
	}

	const TargetParamsType& tstStruct() const { return targetHelper_.tstStruct(); }

	template<typename SomeBasisType>
	void setGs(const typename PsimagLite::Vector<VectorType>::Type& v,
	           const SomeBasisType& someBasis)
	{
		applyOpExpression_.psi().set(v,someBasis);
	}

	RealType normSquared(SizeType i) const
	{
		const VectorWithOffsetType& v = applyOpExpression_.targetVectors()[i];
		// call to mult will conjugate one of the vector
		return std::real(v*v);
	}

	void normalizeTimeVectors(SizeType start = 0, SizeType end = 0)
	{
		SizeType total =  applyOpExpression_.targetVectors().size();
		if (end == 0) end = total;
		for (SizeType i = start; i < end; ++i) {
			RealType factor = normSquared(i);
			if (fabs(factor)<1e-10) {
				std::cerr<<"normalizeTimeVectors ";
				std::cerr<<"vector "<<i<<" too small "<<factor<<"\n";
				continue;
			}

			factor = 1.0/sqrt(factor);
			applyOpExpression_.multiplyTimeVector(i,factor);
		}
	}

	void setTime(RealType t)
	{
		applyOpExpression_.setTime(t);
	}

	void timeHasAdvanced() { applyOpExpression_.timeHasAdvanced(); }

	template<typename IoOutputType>
	void save(const typename PsimagLite::Vector<SizeType>::Type& block,
	          IoOutputType& io,
	          const PostProcType& cf,
	          const VectorVectorWithOffsetType& targetVectors) const
	{
		DynamicSerializerType dynS(cf,block[0],targetVectors);
		dynS.save(io);
	}

	template<typename SomeSerializerType>
	void load(const PsimagLite::String& f)
	{
		IoInputType io(f);
		PsimagLite::String loadInto = targetHelper_.model().params().checkpoint.into;
		PsimagLite::String labelForPsi = targetHelper_.model().params().checkpoint.labelForPsi;

		if (loadInto == "All") {
			setAllStagesTo(WFT_NOADVANCE);

			SomeSerializerType ts(io,IoInputType::LAST_INSTANCE);
			for (SizeType i=0;i<targetVectors().size();i++)
				targetVectors(i) = ts.vector(i);

			applyOpExpression_.setTime(ts.time());

			applyOpExpression_.psi().load(io,labelForPsi);

		} else {
			setAllStagesTo(DISABLED);
			io.rewind();
			int site = 0;
			io.readline(site,"#TCENTRALSITE=",IoType::In::LAST_INSTANCE);
			applyOpExpression_.psi().loadOneSector(io,labelForPsi);
		}
	}

	template<typename SomeSerializerType>
	void load(const PsimagLite::String& f,int)
	{
		IoInputType io(f);

		int site=0;
		io.readline(site,"#TCENTRALSITE=",IoType::In::LAST_INSTANCE);
		if (site<0)
			throw PsimagLite::RuntimeError("GST::load(...): site cannot be negative\n");

		applyOpExpression_.psi().load(io,"PSI");
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
	                  const VectorSizeType& block) const
	{
		PsimagLite::Vector<SizeType>::Type nk;
		setNk(nk,block);
		setInitialVector(v,applyOpExpression_.psi(),nk);
	}

	void computeCorrection(SizeType direction,
	                       const BlockType& block1)
	{
		const VectorWithOffsetType& psi = applyOpExpression_.psi();
		VectorWithOffsetType& v = applyOpExpression_.targetVectors(0);

		// operators in the one-site basis:
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		SparseMatrixType hmatrix;
		SymmetryElectronsSzType q;

		RealType time = 0;
		targetHelper_.model().setNaturalBasis(creationMatrix,hmatrix,q,block1,time);

		FermionSign fs(targetHelper_.lrs().left(),q.electrons());
		for (SizeType j=0;j<creationMatrix.size();j++) {
			VectorWithOffsetType phiTemp;
			applyOpExpression_.applyOpLocal()(phiTemp,psi,creationMatrix[j],
			                                  fs,direction,ApplyOperatorType::BORDER_NO);
			if (j==0) v = phiTemp;
			else v += phiTemp;
		}
	}

	int findFermionSignOfTheOperators() const
	{
		const VectorOperatorType& myoperator = targetHelper_.tstStruct().aOperators();
		bool wereSumming = (targetHelper_.tstStruct().concatenation() == SUM);
		int f = 0;

		for (SizeType i = 0; i < myoperator.size(); ++i) {

			RealType norma = norm2(myoperator[i].data);

			if (norma==0 && wereSumming) continue;
			if (isTheIdentity(myoperator[i].data) && !wereSumming) continue;

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

	void initTimeVectors(const VectorRealType& times,InputValidatorType& ioIn)
	{
		applyOpExpression_.initTimeVectors(times,ioIn);
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

	void cocoon(const BlockType& block,SizeType direction) const
	{
		const ModelType& model = targetHelper_.model();
		const VectorVectorWithOffsetType& tv = applyOpExpression_.targetVectors();

		if (model.params().insitu=="") return;

		if (BasisType::useSu2Symmetry()) {
			noCocoon("not when SU(2) symmetry is in use");
			return;
		}

		SizeType max = 1;
		PsimagLite::String magic = "allPvectors";
		if (model.params().options.find(magic) != PsimagLite::String::npos)
			max = tv.size();

		try {
			assert(block.size()>0);
			cocoon(direction,block[0],psi(),"PSI",psi(),"PSI");
			if (tv.size() > 0) {
				for (SizeType i = 0; i < max; ++i)
					cocoon(direction,block[0],tv[i],"P"+ttos(i),tv[i],"P"+ttos(i));
				for (SizeType i = 0; i < max; ++i)
					cocoon(direction,block[0],psi(),"PSI",tv[i],"P"+ttos(i));
			}
		} catch (std::exception& e) {
			noCocoon("unsupported by the model");
		}
	}

	void cocoonLegacy(SizeType direction,const BlockType& block) const
	{
		const VectorWithOffsetType& psi = applyOpExpression_.psi();
		const VectorWithOffsetType& tv0 = applyOpExpression_.targetVectors()[0];
		const ModelType& model = targetHelper_.model();

		SizeType site = block[0];
		SparseMatrixType tmpC(model.naturalOperator("nup",site,0));
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

		SparseMatrixType tmpC2(model.naturalOperator("ndown",site,0));
		OperatorType ndown(tmpC2,fermionSign1,jm1,angularFactor1,su2Related1);
		test(psi,psi,direction,"<PSI|ndown|PSI>",site,ndown,ApplyOperatorType::BORDER_NO);
		s = "<P0|ndown|P0>";
		test(tv0,tv0,direction,s,site,ndown,ApplyOperatorType::BORDER_NO);

		SparseMatrixType tmpC3 = (nup.data * ndown.data);
		OperatorType doubleOcc(tmpC3,fermionSign1,jm1,angularFactor1,su2Related1);
		test(psi,psi,direction,"<PSI|doubleOcc|PSI>",site,doubleOcc,
		     ApplyOperatorType::BORDER_NO);
		s = "<P0|doubleOcc|P0>";
		test(tv0,tv0,direction,s,site,doubleOcc,ApplyOperatorType::BORDER_NO);
	}

	// in situ computation:
	void cocoon(SizeType direction,
	            SizeType site,
	            const VectorWithOffsetType& v1,
	            PsimagLite::String label1,
	            const VectorWithOffsetType& v2,
	            PsimagLite::String label2) const
	{

		std::cout<<"-------------&*&*&* In-situ measurements start\n";
		RealType norm1 = std::norm(v1);
		RealType norm2 = std::norm(v2);
		if (norm1 < 1e-6 || norm2 < 1e-6) {
			std::cout<<"cocoon: At least 1 NORM IS ZERO ";
			std::cout<<label1<<" has norm "<<norm1;
			std::cout<<" "<<label2<<" has norm "<<norm2<<"\n";
			return;
		}

		cocoon_(direction,site,v1,label1,v2,label2,ApplyOperatorType::BORDER_NO);

		int site2 = findBorderSiteFrom(site,direction);

		if (site2 >= 0) {
			cocoon_(direction,site2,v1,label1,v2,label2,ApplyOperatorType::BORDER_YES);
		}

		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	const ComplexOrRealType& inSitu(SizeType site) const
	{
		assert(site < inSitu_.size());
		return inSitu_[site];
	}

private:

	void noCocoon(const PsimagLite::String& msg) const
	{
		std::cout<<"-------------&*&*&* In-situ measurements start\n";
		std::cout<<"----- NO IN-SITU MEAS. POSSIBLE, reason="<<msg<<"\n";
		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	void setNk(typename PsimagLite::Vector<SizeType>::Type& nk,
	           const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		for (SizeType i=0;i<block.size();i++)
			nk.push_back(targetHelper_.model().hilbertSize(block[i]));
	}

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
	             const VectorWithOffsetType& v1,
	             PsimagLite::String label1,
	             const VectorWithOffsetType& v2,
	             PsimagLite::String label2,
	             BorderEnumType border) const
	{
		VectorStringType vecStr = getOperatorLabels();

		for (SizeType i=0;i<vecStr.size();i++) {
			const PsimagLite::String& opLabel = vecStr[i];
			OperatorType nup = getOperatorForTest(opLabel,site);

			PsimagLite::String tmpStr = "<"+ label1 + "|" + opLabel + "|" + label2 + ">";
			test(v1,v2,direction,tmpStr,site,nup,border);
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

			NaturalOpStruct nos(opLabel);
			SparseMatrixType tmpC(targetHelper_.model().naturalOperator(nos.label,
			                                                            site,
			                                                            nos.dof));
			nup = OperatorType(tmpC,nos.fermionSign,jm1,angularFactor1,su2Related1);
		}

		SizeType foundSize = nup.data.row();
		SizeType expectedSize = targetHelper_.model().hilbertSize(site);
		if (foundSize != expectedSize) {
			PsimagLite::String str("getOperatorForTest ");
			str += " Expected size " + ttos(expectedSize);
			str += " but found size " + ttos(foundSize);
			str += " for operator " + opLabel + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		return nup;
	}

	void test(const VectorWithOffsetType& src1,
	          const VectorWithOffsetType& src2,
	          SizeType systemOrEnviron,
	          PsimagLite::String label,
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
			for (SizeType jj=0;jj<src2.sectors();jj++) {
				SizeType j = src2.sector(jj);
				if (i!=j) continue;
				for (SizeType k=0;k<dest.effectiveSize(i);k++)
					sum+= dest.fastAccess(i,k)*
					        std::conj(src2.fastAccess(j,k));
			}
		}

		assert(site < inSitu_.size());
		inSitu_[site] = sum;
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
	mutable VectorType inSitu_;
}; // class TargetingCommon

template<typename TargetHelperType,
         typename VectorWithOffsetType,
         typename LanczosSolverType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingCommon<TargetHelperType,
                         VectorWithOffsetType,
                         LanczosSolverType>& tst)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // TARGETING_COMMON_H

