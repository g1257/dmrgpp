/*
Copyright (c) 2009-2016, UT-Battelle, LLC
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

/*! \file TargetingCommon.h
 *
 * Functionality used by many targetting classes
 *
 */

#ifndef TARGETING_COMMON_H
#define TARGETING_COMMON_H

#include "Braket.h"
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TimeSerializer.h"
#include "TargetParamsDynamic.h"
#include "VectorWithOffsets.h"
#include "ContinuedFraction.h"
#include <cassert>
#include "ApplyOperatorExpression.h"
#include "IoSimple.h"
#include "PsimagLite.h"

namespace Dmrg {

template<typename TargetHelperType,
         typename VectorWithOffsetType,
         typename LanczosSolverType>
class TargetingCommon  {

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
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
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
	typedef Braket<ModelType> BraketType;

	static const SizeType SUM = TargetParamsType::SUM;

	enum {DISABLED=ApplyOperatorExpressionType::DISABLED,
		  OPERATOR=ApplyOperatorExpressionType::OPERATOR,
		  WFT_NOADVANCE=ApplyOperatorExpressionType::WFT_NOADVANCE};

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
	                ProgramGlobals::DirectionEnum direction,
	                SizeType site,
	                SizeType loopNumber)
	{
		return applyOpExpression_.getPhi(phiNew,Eg,direction,site,loopNumber);
	}

	SizeType getPhi(VectorWithOffsetType& phiNew,
	                const VectorWithOffsetType& phiSrc,
	                RealType Eg,
	                ProgramGlobals::DirectionEnum direction,
	                SizeType site,
	                SizeType loopNumber)
	{
		return applyOpExpression_.getPhi(phiNew,phiSrc,Eg,direction,site,loopNumber);
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
		return PsimagLite::real(v*v);
	}

	void normalizeTimeVectors(SizeType start = 0, SizeType end = 0)
	{
		SizeType total =  applyOpExpression_.targetVectors().size();
		if (end == 0) end = total;
		for (SizeType i = start; i < end; ++i) {
			RealType factor = normSquared(i);
			if (fabs(factor) == 0) continue;

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
	void save(const VectorSizeType& block,
	          IoOutputType& io,
	          const PostProcType& cf,
	          const VectorVectorWithOffsetType& targetVectors) const
	{
		cf.save(io);
		SizeType marker = (noStageIs(DISABLED)) ? 1 : 0;

		TimeSerializerType ts(currentTime(),
		                      block[0],
		        applyOpExpression_.targetVectors(),
		        marker);
		ts.save(io);
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
	                  const VectorSizeType& block,
	                  bool noguess) const
	{
		PsimagLite::Vector<SizeType>::Type nk;
		setNk(nk,block);
		setInitialVector(v,applyOpExpression_.psi(), nk, noguess);
	}

	void computeCorrection(ProgramGlobals::DirectionEnum direction,
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

	const VectorVectorWithOffsetType& targetVectors() const
	{
		return applyOpExpression_.targetVectors();
	}

	VectorWithOffsetType& targetVectors(SizeType i)
	{
		return applyOpExpression_.targetVectors(i);
	}

	const VectorWithOffsetType& targetVectors(SizeType i) const
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
	                     ProgramGlobals::DirectionEnum direction,
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

	void applyOneOperator(SizeType loopNumber,
	                      SizeType indexOfOperator,
	                      SizeType site,
	                      VectorWithOffsetType& phiNew,
	                      const VectorWithOffsetType& psiSrc,
	                      SizeType systemOrEnviron)
	{
		applyOpExpression_.applyOneOperator(loopNumber,
		                                    indexOfOperator,
		                                    site,
		                                    phiNew,
		                                    psiSrc,
		                                    systemOrEnviron);
	}

	void applyOneOperator(SizeType loopNumber,
	                      SizeType indexOfOperator,
	                      SizeType site,
	                      VectorWithOffsetType& phiNew,
	                      SizeType systemOrEnviron)
	{
		applyOpExpression_.applyOneOperator(loopNumber,
		                                    indexOfOperator,
		                                    site,
		                                    phiNew,
		                                    applyOpExpression_.psi(),
		                                    systemOrEnviron);
	}

	void wftOneVector(VectorWithOffsetType& phiNew,
	                  const VectorWithOffsetType& src,
	                  SizeType site)
	{
		applyOpExpression_.wftOneVector(phiNew, src, site);
	}

	void wftAll(SizeType site)
	{
		applyOpExpression_.wftAll(site);
	}

	void cocoon(const BlockType& block,
	            ProgramGlobals::DirectionEnum direction) const
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

	void cocoonLegacy(ProgramGlobals::DirectionEnum direction,const BlockType& block) const
	{
		const VectorWithOffsetType& psi = applyOpExpression_.psi();
		const VectorWithOffsetType& tv0 = applyOpExpression_.targetVectors()[0];
		const ModelType& model = targetHelper_.model();

		SizeType site = block[0];
		OperatorType nup = model.naturalOperator("nup",site,0);

		test(psi,psi,direction,"<PSI|nup|PSI>",site,nup,ApplyOperatorType::BORDER_NO);
		PsimagLite::String s = "<P0|nup|P0>";
		test(tv0,tv0,direction,s,site,nup,ApplyOperatorType::BORDER_NO);

		OperatorType ndown = model.naturalOperator("ndown",site,0);
		test(psi,psi,direction,"<PSI|ndown|PSI>",site,ndown,ApplyOperatorType::BORDER_NO);
		s = "<P0|ndown|P0>";
		test(tv0,tv0,direction,s,site,ndown,ApplyOperatorType::BORDER_NO);

		SparseMatrixType tmpC3 = (nup.data * ndown.data);
		OperatorType doubleOcc(tmpC3,
		                       nup.fermionSign,
		                       nup.jm,
		                       nup.angularFactor,
		                       nup.su2Related);
		test(psi,psi,direction,"<PSI|doubleOcc|PSI>",site,doubleOcc,
		     ApplyOperatorType::BORDER_NO);
		s = "<P0|doubleOcc|P0>";
		test(tv0,tv0,direction,s,site,doubleOcc,ApplyOperatorType::BORDER_NO);
	}

	// in situ computation:
	// prints <v1|A|v2>
	void cocoon(ProgramGlobals::DirectionEnum direction,
	            SizeType site,
	            const VectorWithOffsetType& v1,
	            PsimagLite::String label1,
	            const VectorWithOffsetType& v2,
	            PsimagLite::String label2,
	            bool needsPrinting = true) const
	{

		std::cout<<"-------------&*&*&* In-situ measurements start\n";
		RealType norm1 = norm(v1);
		RealType norm2 = norm(v2);
		if (norm1 < 1e-6 || norm2 < 1e-6) {
			std::cout<<"cocoon: At least 1 NORM IS ZERO ";
			std::cout<<label1<<" has norm "<<norm1;
			std::cout<<" "<<label2<<" has norm "<<norm2<<"\n";
			return;
		}

		BorderEnumType border = ApplyOperatorType::BORDER_NO;
		cocoon_(direction,site,v1,label1,v2,label2,border,needsPrinting);

		SizeType numberOfSites = targetHelper_.model().geometry().numberOfSites();

		int site2 = ProgramGlobals::findBorderSiteFrom(site, direction, numberOfSites);

		if (site2 >= 0) {
			border = ApplyOperatorType::BORDER_YES;
			cocoon_(direction,site2,v1,label1,v2,label2,border,needsPrinting);
		}

		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	void calcBracket(ProgramGlobals::DirectionEnum direction,
	                 SizeType site,
	                 const BraketType& braket) const
	{
		if (braket.points() != 1)
			err("Brakets in situ must be one-points\n");

		const VectorWithOffsetType& v1 = getVector(braket.bra());
		const VectorWithOffsetType& v2 = getVector(braket.ket());
		std::cout<<"-------------&*&*&* In-situ measurements start\n";
		RealType norm1 = norm(v1);
		RealType norm2 = norm(v2);
		if (norm1 < 1e-6 || norm2 < 1e-6) {
			std::cout<<"cocoon: At least 1 NORM IS ZERO ";
			std::cout<<braket.bra()<<" has norm "<<norm1;
			std::cout<<" "<<braket.ket()<<" has norm "<<norm2<<"\n";
			return;
		}

		BorderEnumType border = ApplyOperatorType::BORDER_NO;
		test(v1,v2,direction,braket.toString(),site,braket.op(0),border);

		SizeType numberOfSites = targetHelper_.model().geometry().numberOfSites();

		int site2 = ProgramGlobals::findBorderSiteFrom(site, direction, numberOfSites);

		if (site2 >= 0) {
			border = ApplyOperatorType::BORDER_YES;
			test(v1,v2,direction,braket.toString(),site2,braket.op(0),border);
		}

		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	ComplexOrRealType rixsCocoon(ProgramGlobals::DirectionEnum direction,
	                             SizeType site,
	                             SizeType index1,
	                             SizeType index2) const
	{
		ComplexOrRealType value = 0.0;
		VectorStringType vecStr = getOperatorLabels();
		if (vecStr.size() == 0) return value;
		if (vecStr.size() > 1) {
			throw PsimagLite::RuntimeError("rixsCocoon: supports only 1 insitu operator\n");
		}

		SizeType numberOfSites = targetHelper_.model().geometry().numberOfSites();
		BorderEnumType border = (site == 0 || site == numberOfSites - 1) ?
		            ApplyOperatorType::BORDER_YES : ApplyOperatorType::BORDER_NO;

		const VectorWithOffsetType& v1 =  applyOpExpression_.targetVectors(index1);
		const VectorWithOffsetType& v2 =  applyOpExpression_.targetVectors(index2);

		assert(vecStr.size() == 1);
		for (SizeType i=0;i<vecStr.size();i++) {
			PsimagLite::String opLabel = vecStr[i];

			BraketType Braket(targetHelper_.model(),"<gs|"+opLabel+"[" + ttos(site) + "]|gs>");

			OperatorType A = Braket.op(0);

			value = test_(v1,v2,direction,site,A,border);
		}

		return value;
	}


	const ComplexOrRealType& inSitu(SizeType site) const
	{
		assert(site < inSitu_.size());
		return inSitu_[site];
	}

private:

	void setQuantumNumbers(const VectorWithOffsetType& v)
	{
		applyOpExpression_.setQuantumNumbers(v);
	}

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
	                      const VectorSizeType& nk,
	                      bool noguess) const
	{
		const WaveFunctionTransfType& wft = targetHelper_.wft();

		if (noguess)
			wft.createRandomVector(v1);
		else
			wft.setInitialVector(v1,v2,targetHelper_.lrs(),nk);
	}

	// prints <v1|A|v2>
	void cocoon_(ProgramGlobals::DirectionEnum direction,
	             SizeType site,
	             const VectorWithOffsetType& v1,
	             PsimagLite::String label1,
	             const VectorWithOffsetType& v2,
	             PsimagLite::String label2,
	             BorderEnumType border,
	             bool wantsPrinting) const
	{
		VectorStringType vecStr = getOperatorLabels();

		for (SizeType i=0;i<vecStr.size();i++) {
			PsimagLite::String opLabel = braketIfNeeded(vecStr[i],
			                                            site,
			                                            label1,
			                                            label2);

			BraketType Braket(targetHelper_.model(), opLabel);

			OperatorType nup = Braket.op(0);

			if (wantsPrinting) test(v1,v2,direction,opLabel,site,nup,border);
			else test_(v1,v2,direction,site,nup,border);
		}
	}

	PsimagLite::String braketIfNeeded(PsimagLite::String opLabel,
	                                  SizeType site,
	                                  PsimagLite::String label1,
	                                  PsimagLite::String label2) const
	{
		if (label1 == "PSI") label1 = "gs";
		if (label2 == "PSI") label2 = "gs";
		if (opLabel.length() == 0 || opLabel[0] == '<') return opLabel;
		return "<"+label1 + "|" + opLabel + "|" + label2 + ">";
	}

	VectorStringType getOperatorLabels() const
	{
		VectorStringType vecStr;
		PsimagLite::split(vecStr, targetHelper_.model().params().insitu, ",");
		return vecStr;
	}

	const VectorWithOffsetType& getVector(PsimagLite::String braOrKet) const
	{
		if (braOrKet == "gs")
			return applyOpExpression_.psi();

		int ind = BraketType::getPtype(braOrKet);
		if (ind <= 0)
			err("Malformed braket " + braOrKet + "\n");

		return applyOpExpression_.targetVectors(ind - 1);
	}

	// prints <src2|A|src1>
	void test(const VectorWithOffsetType& src1,
	          const VectorWithOffsetType& src2,
	          SizeType systemOrEnviron,
	          PsimagLite::String label,
	          SizeType site,
	          const OperatorType& A,
	          BorderEnumType border) const
	{
		ComplexOrRealType sum = test_(src1,src2,systemOrEnviron,site,A,border);
		std::cout<<site<<" "<<sum<<" "<<currentTime();
		std::cout<<" "<<label<<" "<<(src1*src2)<<"\n";
	}

	// returns <src2|A|src1>
	ComplexOrRealType test_(const VectorWithOffsetType& src1,
	                        const VectorWithOffsetType& src2,
	                        SizeType systemOrEnviron,
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
					        PsimagLite::conj(src2.fastAccess(j,k));
			}
		}

		assert(site < inSitu_.size());
		inSitu_[site] = sum;
		return sum;
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

