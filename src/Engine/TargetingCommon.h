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

/*! \file TargetingCommon.h
 *
 * Functionality used by many targeting classes
 *
 */

#ifndef TARGETING_COMMON_H
#define TARGETING_COMMON_H

#include "Braket.h"
#include "ProgressIndicator.h"
#include "TimeSerializer.h"
#include "TargetParamsDynamic.h"
#include "VectorWithOffsets.h"
#include "ContinuedFraction.h"
#include <cassert>
#include "ApplyOperatorExpression.h"
#include "Io/IoSelector.h"
#include "PsimagLite.h"
#include "GetBraOrKet.h"
#include "RestartStruct.h"
#include "SdhsReinterpret.h"

namespace Dmrg {

template<typename LambdaType>
class LambdaForTests {

public:

	void operator()(const PsimagLite::GetBraOrKet& bra,
	                const PsimagLite::GetBraOrKet& ket,
	                LambdaType* t) const
	{
		if (!t) return;
		t->operator()(bra, ket);
	}
};

template<>
class LambdaForTests<int> {

public:

	void operator()(const PsimagLite::GetBraOrKet&,
	                const PsimagLite::GetBraOrKet&,
	                int*) const
	{}
};

template<typename TargetHelperType,
         typename VectorWithOffsetType_,
         typename LanczosSolverType_>
class TargetingCommon  {

public:

	enum SetTvsEnum { NO_TVS = false, READ_AND_SET_TVS = true};

	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef LanczosSolverType_ LanczosSolverType;
	typedef PsimagLite::IoSelector IoType;
	typedef typename IoType::In IoInputType;
	typedef typename TargetHelperType::RealType RealType;
	typedef typename TargetHelperType::ModelType ModelType;
	typedef typename TargetHelperType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::VectorQnType VectorQnType;
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
	typedef typename ApplyOperatorExpressionType::TargetParamsType TargetParamsType;
	typedef typename ApplyOperatorExpressionType::VectorVectorWithOffsetType
	VectorVectorWithOffsetType;
	typedef typename ApplyOperatorExpressionType::VectorRealType VectorRealType;
	typedef typename ApplyOperatorExpressionType::PairType PairType;
	typedef typename ApplyOperatorExpressionType::TimeVectorsBaseType TimeVectorsBaseType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef Braket<ModelType> BraketType;
	typedef FermionSign FermionSignType;
	typedef typename ApplyOperatorExpressionType::StageEnumType StageEnumType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef RestartStruct RestartStructType;
	typedef typename ApplyOperatorExpressionType::VectorVectorVectorWithOffsetType
	VectorVectorVectorWithOffsetType;
	typedef SdhsReinterpret<BraketType> SdhsReinterpretType;

	enum class OpLabelCategory { DRESSED, BARE };

	TargetingCommon(const LeftRightSuperType& lrs,
	                const ModelType& model,
	                const WaveFunctionTransfType& wft,
	                SizeType indexNoAdvance)
	    : cocoonType_(OpLabelCategory::DRESSED),
	      progress_("TargetingCommon"),
	      targetHelper_(lrs, model, wft),
	      aoe_(targetHelper_, indexNoAdvance),
	      inSitu_(model.superGeometry().numberOfSites())
	{
		PsimagLite::split(meas_, model.params().insitu, ",");
		SizeType n = meas_.size();
		for (SizeType i = 0; i < n; ++i) {
			const bool isDressed = isOpLabelDressed(meas_[i]);

			// check this early that what's passed makes sense
			if (isDressed)
				BraketType(model, meas_[i]);

			OpLabelCategory cocoonExpected = (isDressed) ? OpLabelCategory::DRESSED
			                                             : OpLabelCategory::BARE;
			if (i == 0) {
				cocoonType_ = cocoonExpected;
				continue;
			}

			if (cocoonType_ != cocoonExpected)
				err("FATAL: If one label is dressed (bare) then all must be dressed (bare)\n");
		}
	}

	void postCtor(SizeType tstSites, SizeType targets)
	{
		aoe_.postCtor(tstSites);
		aoe_.targetVectorsResize(targets);
	}

	// START read/write

	void write(PsimagLite::IoSelector::Out& io,
	           const VectorSizeType& block,
	           PsimagLite::String prefix) const
	{
		if (block.size() != 1)
			err(PsimagLite::String(__FILE__) + " write() only supports blocks.size=1\n");

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Saving state...";
		progress_.printline(msgg, std::cout);

		io.write(block[0], prefix + "/TargetCentralSite");
		aoe_.writePsi(io, prefix);
	}

	void writeNGSTs(PsimagLite::IoSelector::Out& io,
	                PsimagLite::String prefix,
	                const VectorSizeType& block,
	                PsimagLite::String name,
	                const PostProcType& cf) const
	{
		cf.write(io, prefix);
		writeNGSTs(io, prefix, block, name);
	}

	void writeNGSTs(PsimagLite::IoSelector::Out& io,
	                PsimagLite::String prefix,
	                const VectorSizeType& block,
	                PsimagLite::String name) const
	{
		SizeType site = block[0];
		TimeSerializerType ts(currentTimeStep(),
		                      time(),
		                      site,
		                      aoe_,
		                      name);
		ts.write(io, prefix);
	}

	void read(IoInputType& io, PsimagLite::String prefix)
	{
		const PsimagLite::String label = targetHelper_.model().params().
		        checkpoint.labelForEnergy();
		prefix += "/";
		aoe_.loadEnergy(io, label);
		aoe_.readPsi(io, prefix);
	}

	void readGSandNGSTs(IoInputType& io,
	                    PsimagLite::String prefix,
	                    PsimagLite::String name)
	{
		read(io, prefix);

		TimeSerializerType* ts = 0;

		try {
			ts = new TimeSerializerType(io, prefix);
		} catch (...) {
			return;
		}

		const typename TimeSerializerType::VectorStageEnumType& stages = ts->stages();
		const SizeType rstages = stages.size();        // read stages
		const SizeType dstages = aoe_.stages().size(); // destination stages

		const RestartStructType& checkpoint = targetHelper_.model().params().checkpoint;

		const SizeType dstagesOrZero = (checkpoint.mapStages()) ? dstages : 0;

		if (dstagesOrZero > 0 && dstagesOrZero != rstages)
			err("TargetingCommon::readGSandNGSTs: stages cannot be set from previous\n");

		for (SizeType i = 0; i < dstagesOrZero; ++i)
			aoe_.setStage(i, stages[i]);

		SizeType rtvs = ts->numberOfVectors();        // read tvs
		SizeType dtvs = aoe_.tvs();  // destination tvs

		int tvForPsi = checkpoint.sourceTvForPsi();
		if (tvForPsi >= 0) {
			SizeType tvForPsiUnsigned = tvForPsi;
			if (tvForPsiUnsigned >= rtvs)
				err("TargetingCommon::readGSandNGSTs: sourceTvForPsi >= " + ttos(rtvs) + "\n");

			std::cout<<"FIXME TODO WARNING: Need better spec for TvForPsi\n";
			std::cerr<<"FIXME TODO WARNING: Need better spec for TvForPsi\n";

			aoe_.setOnlyOnePsi(ts->vector(tvForPsiUnsigned));
		}

		ApplyOperatorExpressionType& aoeNonConst = const_cast<ApplyOperatorExpressionType&>(aoe_);
		for (SizeType i = 0; i < dtvs; ++i) {
			const int j = checkpoint.mappingTvs(i);
			if (j < 0) continue;
			const SizeType jj = j;
			if (jj >= rtvs) {
				err("TargetingCommon::readGSandNGSTs: tvs mapping failed " +
				    ttos(j) + " >= " + ttos(rtvs) + "\n");
			}

			aoeNonConst.targetVectorsNonConst(i) = ts->vector(j);
		}

		// FIXME TODO check that the NGST name changes instead
		bool sameNgst = isThisNgstSameAsPrevious(name, ts->name(), dtvs, rtvs);
		SizeType cTimeStep = (sameNgst) ? ts->currentTimeStep() : 0;
		setCurrentTimeStep(cTimeStep);

		RealType timeReal = (sameNgst) ? ts->time() : 0;
		setCurrentTime(timeReal);

		delete ts;
		ts = 0;
	}

	// END read/write

	const ApplyOperatorExpressionType& aoe() const { return aoe_; }

	ApplyOperatorExpressionType& aoeNonConst() { return aoe_; }

	RealType time() const
	{
		return (aoe_.hasTimeVectors()) ? aoe_.timeVectors().time() : 0;
	}

	// START Cocoons

	template<typename SomeLambdaType = int>
	void cocoon(const BlockType& block,
	            ProgramGlobals::DirectionEnum direction,
	            bool doBorderIfBorder,
	            SomeLambdaType* someLambda = nullptr) const
	{
		if (aoe_.noStageIs(StageEnumType::DISABLED))
			std::cout<<"ALL OPERATORS HAVE BEEN APPLIED\n";
		else
			std::cout<<"NOT ALL OPERATORS APPLIED YET\n";

		if (cocoonType_ == OpLabelCategory::BARE)
			err("TargetingCommon.h: cocoon(): Bare specs only for RIXS\n");

		SizeType n = meas_.size();
		assert(block.size()>0);
		SizeType site = block[0];
		SizeType numberOfSites = targetHelper_.model().superGeometry().numberOfSites();
		BorderEnumType border = ApplyOperatorType::BORDER_NO;
		if (site == 0 &&
		        direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON)
			border = ApplyOperatorType::BORDER_YES;
		if (site == numberOfSites - 1 &&
		        direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			border = ApplyOperatorType::BORDER_YES;

		const SizeType expectedSize = targetHelper_.model().hilbertSize(site);

		LambdaForTests<SomeLambdaType> lambdaForTests;

		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String opLabel = meas_[i];

			BraketType braket(targetHelper_.model(), opLabel);

			SdhsReinterpretType sdhs(braket, {site});

			if (sdhs.forbidden()) continue;

			const typename BraketType::AlgebraType& nup = sdhs.op(0);

			const SizeType foundSize = nup.getStorage().rows();
			if (foundSize != expectedSize) continue;

			const VectorWithOffsetType* v1 = getVector(braket.bra());
			const VectorWithOffsetType* v2 = getVector(braket.ket());

			if (!v1 || !v2) continue;

			if (v1->size() == 0 || v2->size() == 0) continue;

			lambdaForTests(braket.bra(), braket.ket(), someLambda);
			test(*v1, *v2, direction, opLabel, site, nup, border);
			// don't repeat for border because this is called twice if needed
		}

		if (!doBorderIfBorder) return;

		SizeType site2 = numberOfSites;

		if (site == 1 && direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON)
			site2 = 0;
		if (site == numberOfSites - 2 && direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			site2 = numberOfSites - 1;
		if (site2 == numberOfSites) return;

		VectorSizeType block1(1, site2);
		cocoon(block1, direction, false);
	}

	// FIXME TODO REMOVE
	ComplexOrRealType rixsCocoon(ProgramGlobals::DirectionEnum direction,
	                             SizeType site,
	                             SizeType index1,
	                             SizeType index2,
	                             bool needsShift) const
	{
		const ModelType& model = targetHelper_.model();
		SizeType h = model.hilbertSize(site);
		typename OperatorType::Su2RelatedType su2Related1;
		SparseMatrixType idSparse;
		idSparse.makeDiagonal(h, 1.0);
		OperatorType id(idSparse,
		                ProgramGlobals::FermionOrBosonEnum::BOSON,
		                PairType(0, 0),
		                1.0,
		                su2Related1);
		ComplexOrRealType value = 0.0;
		SizeType n = meas_.size();
		if (n== 0) return value;
		if (n > 1)
			err("rixsCocoon: supports only 1 insitu operator\n");

		if (cocoonType_ != OpLabelCategory::BARE)
			err("rixsCocoon: supports only bare operators\n");

		SizeType numberOfSites = targetHelper_.model().superGeometry().numberOfSites();
		BorderEnumType border = (site == 0 || site == numberOfSites - 1) ?
		            ApplyOperatorType::BORDER_YES : ApplyOperatorType::BORDER_NO;

		const VectorWithOffsetType& v1 =  aoe_.targetVectors(index1);
		const VectorWithOffsetType& v2 =  aoe_.targetVectors(index2);

		assert(n == 1);
		for (SizeType i = 0; i < n; ++i) {
			const PsimagLite::String opLabel = meas_[i];

			BraketType Braket(targetHelper_.model(),
			                  "<gs|"+opLabel+"[" + ttos(site) + "]|gs>");
			if (needsShift) {
				const typename BraketType::AlgebraType& A = Braket.op(0);
				value = test_(v1, v2, direction, site, A, border);
			} else {
				value = test_(v1, v2, direction, site, id, border);
			}
		}

		return value;
	}

	// END Cocoons

	void setAllStagesTo(StageEnumType x)
	{
		SizeType n = aoe_.stages().size();
		for (SizeType i = 0; i < n; ++i)
			aoe_.setStage(i, x);
	}

	RealType normSquared(SizeType i) const
	{
		const VectorWithOffsetType& v = aoe_.targetVectors(i);
		if (v.size() == 0) return 0;
		// call to mult will conjugate one of the vector
		return PsimagLite::real(v*v);
	}

	void normalizeTimeVectors(SizeType start = 0, SizeType end = 0)
	{
		SizeType total =  aoe_.tvs();
		if (end == 0) end = total;
		for (SizeType i = start; i < end; ++i) {
			RealType factor = normSquared(i);
			if (fabs(factor) == 0) continue;

			factor = 1.0/sqrt(factor);
			aoe_.multiplyTimeVector(i,factor);
		}
	}

	void initialGuess(VectorWithOffsetType& v,
	                  const VectorWithOffsetType& psi,
	                  const VectorSizeType& block,
	                  bool noguess) const
	{
		PsimagLite::Vector<SizeType>::Type nk;
		setNk(nk, block);
		const WaveFunctionTransfType& wft = targetHelper_.wft();

		if (noguess)
			wft.createRandomVector(v);
		else
			wft.setInitialVector(v, psi, targetHelper_.lrs(), nk);
	}

	void computeCorrection(ProgramGlobals::DirectionEnum direction,
	                       const BlockType& block1)
	{
		const VectorWithOffsetType& psi = aoe_.ensureOnlyOnePsi("computeCorrection");
		ApplyOperatorExpressionType& aoeNonConst = const_cast<ApplyOperatorExpressionType&>(aoe_);
		VectorWithOffsetType& v = aoeNonConst.targetVectorsNonConst(0);

		// operators in the one-site basis:
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		VectorQnType q;
		targetHelper_.model().setOperatorMatrices(creationMatrix, q, block1);

		typename BasisWithOperatorsType::VectorBoolType signs(q.size());
		for (SizeType i = 0; i < q.size(); ++i) signs[i] = q[i].oddElectrons;

		assert(0 < block1.size());
		const SizeType splitSize = targetHelper_.model().hilbertSize(block1[0]);

		FermionSign fs(targetHelper_.lrs().left(), signs);
		for (SizeType j=0;j<creationMatrix.size();j++) {
			VectorWithOffsetType phiTemp;
			const OperatorType& cm = creationMatrix[j];
			aoe_.applyOpLocal()(phiTemp,
			                    psi,
			                    cm,
			                    fs,
			                    splitSize,
			                    direction,
			                    ApplyOperatorType::BORDER_NO);
			if (j==0) v = phiTemp;
			else v += phiTemp;
		}
	}

	int findFermionSignOfTheOperators(typename TargetParamsType::ConcatEnum concat,
	                                  const VectorOperatorType& myoperator) const
	{
		bool wereSumming = (concat == TargetParamsType::ConcatEnum::SUM);
		ProgramGlobals::FermionOrBosonEnum forB = ProgramGlobals::FermionOrBosonEnum::BOSON;

		for (SizeType i = 0; i < myoperator.size(); ++i) {

			RealType norma = norm2(myoperator[i].getCRS());

			if (norma==0 && wereSumming) continue;
			if (isTheIdentity(myoperator[i].getCRS()) && !wereSumming) continue;

			if (i == 0) {
				forB = myoperator[i].fermionOrBoson();
				continue;
			}

			if (forB == myoperator[i].fermionOrBoson()) continue;

			PsimagLite::String str("CorrectionVectorTargeting: ");
			str += "inconsistent sign for operators\n";
			throw PsimagLite::RuntimeError(str);

		}

		return (forB == ProgramGlobals::FermionOrBosonEnum::FERMION) ? -1 : 1;
	}

	const typename VectorWithOffsetType_::value_type& inSitu(SizeType site) const
	{
		assert(site < inSitu_.size());
		return inSitu_[site];
	}

	void calcBracket(ProgramGlobals::DirectionEnum direction,
	                 SizeType site,
	                 const BraketType& braket) const
	{
		if (braket.points() != 1)
			err("Brakets in situ must be one-points\n");

		const VectorWithOffsetType* v1 = getVector(braket.bra());
		const VectorWithOffsetType* v2 = getVector(braket.ket());

		if (!v1 || !v2) return;

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
		// v1 == bra; v2 = ket
		test(v1,v2,direction,braket.toString(),site,braket.op(0),border);

		SizeType numberOfSites = targetHelper_.model().geometry().numberOfSites();

		int site2 = ProgramGlobals::findBorderSiteFrom(site, direction, numberOfSites);

		if (site2 >= 0) {
			border = ApplyOperatorType::BORDER_YES;
			// v1 == bra; v2 = ket
			test(*v1, *v2, direction, braket.toString(), site2, braket.op(0), border);
		}

		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	bool withLegacyBugs() const { return targetHelper_.withLegacyBugs(); }

	// returns <src2|A|src1>
	template<typename SomeAlgebraType>
	ComplexOrRealType testRealWork(const VectorWithOffsetType& src1,
	                               const VectorWithOffsetType& src2,
	                               const ProgramGlobals::DirectionEnum systemOrEnviron,
	                               SizeType site,
	                               const SomeAlgebraType& A,
	                               BorderEnumType border) const
	{
		const SizeType splitSize = targetHelper_.model().hilbertSize(site);
		typename PsimagLite::Vector<bool>::Type oddElectrons;
		targetHelper_.model().findOddElectronsOfOneSite(oddElectrons, site);
		FermionSign fs(targetHelper_.lrs().left(), oddElectrons);
		VectorWithOffsetType dest;
		aoe_.applyOpLocal()(dest, src1, A, fs, splitSize, systemOrEnviron, border);

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

	void printNormsAndWeights(const RealType& gsWeight,
	                          const VectorRealType& weights) const
	{
		if (aoe_.allStages(StageEnumType::DISABLED))
			return;

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"gsWeight="<<gsWeight<<" weights= ";
		for (SizeType i = 0; i < weights.size(); ++i)
			msg<<weights[i]<<" ";
		progress_.printline(msgg, std::cout);

		PsimagLite::OstringStream msgg2(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
		msg2<<"gsNorms=";

		const SizeType nsectors = aoe_.psiConst().size();

		for (SizeType sectorIndex = 0; sectorIndex < nsectors; ++sectorIndex) {
			const SizeType nexcited = aoe_.psiConst()[sectorIndex].size();
			for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex) {
				if (nexcited > 1)
					msg2<<"e"<<excitedIndex<<";";
				if (nsectors > 1)
					msg2<<"s"<<sectorIndex<<";";
				msg2<<norm(*(aoe_.psiConst()[sectorIndex][excitedIndex]))<<" ";
			}
		}

		msg2<<" norms= ";
		for (SizeType i = 0; i < weights.size(); i++)
			msg2<<normSquared(i)<<" ";
		progress_.printline(msgg2, std::cout);

		if (currentTimeStep() == 0) return;

		PsimagLite::OstringStream msgg3(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg3 = msgg3();
		msg3<<"CurrentTime="<<time();
		msg3<<" CurrentTimeStep="<<currentTimeStep();
		progress_.printline(msgg3, std::cout);
	}

private:

	SizeType currentTimeStep() const
	{
		return (aoe_.hasTimeVectors()) ? aoe_.timeVectors().currentTimeStep() : 0;
	}

	void setCurrentTimeStep(SizeType ts)
	{
		TimeVectorsBaseType* ptr = const_cast<TimeVectorsBaseType*>(&aoe_.timeVectors());
		ptr->setCurrentTimeStep(ts);
	}

	void setCurrentTime(RealType t)
	{
		TimeVectorsBaseType* ptr = const_cast<TimeVectorsBaseType*>(&aoe_.timeVectors());
		ptr->setCurrentTime(t);
	}

	bool isThisNgstSameAsPrevious(PsimagLite::String nameThis,
	                              PsimagLite::String namePrev,
	                              SizeType tvsThis,
	                              SizeType tvsPrev) const
	{
		if (nameThis != "LEGACY" && namePrev != "LEGACY")
			return (nameThis == namePrev);

		return (tvsThis == tvsPrev);
	}

	void setQuantumNumbers(const VectorWithOffsetType& v)
	{
		aoe_.setQuantumNumbers(v);
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
		SizeType n = meas_.size();
		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String opLabel = braketTheBare(meas_[i],
			                                           site,
			                                           label1,
			                                           label2);

			BraketType Braket(targetHelper_.model(), opLabel);

			const typename BraketType::AlgebraType& nup = Braket.op(0);

			if (wantsPrinting) test(v1,v2,direction,opLabel,site,nup,border);
			else test_(v1,v2,direction,site,nup,border);
		}
	}

	static PsimagLite::String braketTheBare(PsimagLite::String opLabel,
	                                        SizeType,
	                                        PsimagLite::String label1,
	                                        PsimagLite::String label2)
	{
		if (label1 == "PSI") label1 = "gs";
		if (label2 == "PSI") label2 = "gs";
		if (opLabel.length() == 0 || opLabel[0] == '<')
			err("Expecting bare opspec, got empty or dressed: " + opLabel + "\n");

		return "<"+label1 + "|" + opLabel + "|" + label2 + ">";
	}

	static bool isOpLabelDressed(PsimagLite::String opLabel)
	{
		SizeType n = opLabel.length();
		if (n < 2) return false;
		assert(n > 1);
		return (opLabel[0] == '<' && opLabel[n - 1] == '>');
	}

	const VectorWithOffsetType* getVector(const PsimagLite::GetBraOrKet& getBraOrKet) const
	{
		if (getBraOrKet.isPvector()) {
			const SizeType pIndex = getBraOrKet.pIndex();
			return (pIndex >= aoe_.tvs()) ? nullptr
			                              : &(aoe_.targetVectors(pIndex));
		}

		const SizeType sectorIndex = getBraOrKet.sectorIndex();
		if (aoe_.psiConst().size() <= sectorIndex) {
			std::cout<<"getVector: sectorIndex = "<<sectorIndex<<">=";
			std::cout<<aoe_.psiConst().size()<<" (skipping)\n";
			return nullptr;
		}

		const SizeType levelIndex = getBraOrKet.levelIndex();
		if (aoe_.psiConst()[sectorIndex].size() <= levelIndex)
			return nullptr;

		return aoe_.psiConst()[sectorIndex][levelIndex];
	}

	// prints <src2|A|src1>
	template<typename SomeAlgebraType>
	void test(const VectorWithOffsetType& src1,
	          const VectorWithOffsetType& src2,
	          const ProgramGlobals::DirectionEnum systemOrEnviron,
	          PsimagLite::String label,
	          SizeType site,
	          const SomeAlgebraType& A,
	          BorderEnumType border) const
	{
		ComplexOrRealType sum = test_(src1,src2,systemOrEnviron,site,A,border);
		std::cout<<site<<" "<<sum<<" "<<time();
		std::cout<<" "<<label<<" "<<(src1*src2)<<"\n";
	}

	// returns <src2|A|src1>; but if !withLegacyBugs returns <src1|A|src2>
	template<typename SomeAlgebraType>
	ComplexOrRealType test_(const VectorWithOffsetType& src1,
	                        const VectorWithOffsetType& src2,
	                        const ProgramGlobals::DirectionEnum systemOrEnviron,
	                        SizeType site,
	                        const SomeAlgebraType& A,
	                        BorderEnumType border) const
	{
		if (targetHelper_.withLegacyBugs())
			return testRealWork(src1, src2, systemOrEnviron, site, A, border);
		else
			return testRealWork(src2, src1, systemOrEnviron, site, A, border);
	}

	TargetingCommon(const TargetingCommon&) = delete;

	TargetingCommon& operator=(const TargetingCommon&) = delete;

	OpLabelCategory cocoonType_;
	VectorStringType meas_;
	PsimagLite::ProgressIndicator progress_;
	TargetHelperType targetHelper_;
	ApplyOperatorExpressionType aoe_;
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

