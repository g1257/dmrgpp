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

/*! \file TimeVectorsChebyshev.h
 *
 *
 */

#ifndef TIME_VECTORS_CHEBYSHEV
#define TIME_VECTORS_CHEBYSHEV
#include <iostream>
#include <vector>
#include "TimeVectorsBase.h"
#include "ParallelTriDiag.h"
#include "NoPthreadsNg.h"
#include "Parallelizer.h"
#include "ScaledHamiltonian.h"
#include "Sort.h"
#include "ExpressionCalculator.h"
#include "PredicateAwesome.h"

namespace Dmrg {

template<typename TargetParamsType,
         typename ModelType,
         typename WaveFunctionTransfType,
         typename LanczosSolverType,
         typename VectorWithOffsetType>
class TimeVectorsChebyshev : public  TimeVectorsBase<TargetParamsType,
        ModelType,
        WaveFunctionTransfType,
        LanczosSolverType,
        VectorWithOffsetType> {

	typedef TimeVectorsBase<TargetParamsType,
	ModelType,
	WaveFunctionTransfType,
	LanczosSolverType,
	VectorWithOffsetType> BaseType;
	typedef typename BaseType::PairType PairType;
	typedef typename TargetParamsType::RealType RealType;
	typedef typename BaseType::VectorRealType VectorRealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType
	BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Matrix<ComplexOrRealType> MatrixRealType;
	typedef ParallelTriDiag<ModelType,LanczosSolverType,VectorWithOffsetType>
	ParallelTriDiagType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType
	MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::TargetVectorType TargetVectorType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType
	VectorMatrixFieldType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type
	VectorVectorWithOffsetType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef typename LanczosSolverType::MatrixType MatrixLanczosType;
	typedef ScaledHamiltonian<MatrixLanczosType, TargetParamsType> ScaledMatrixType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::ExpressionCalculator<int> ExpressionCalculatorType;
	typedef PsimagLite::PrepassData<int> PrepassDataType;
	typedef PsimagLite::ExpressionPrepass<PrepassDataType> ExpressionPrepassType;
	typedef PsimagLite::PredicateAwesome<> PredicateAwesomeType;

public:

	TimeVectorsChebyshev(const SizeType& currentTimeStep,
	                     const TargetParamsType& tstStruct,
	                     const VectorRealType& times,
	                     VectorVectorWithOffsetType& targetVectors,
	                     const ModelType& model,
	                     const WaveFunctionTransfType& wft,
	                     const LeftRightSuperType& lrs,
	                     InputValidatorType& ioIn)
	    : BaseType(model, lrs, wft),
	      currentTimeStep_(currentTimeStep),
	      tstStruct_(tstStruct),
	      times_(times),
	      targetVectors_(targetVectors),
	      model_(model),
	      wft_(wft),
	      lrs_(lrs),
	      ioIn_(ioIn),
	      timeHasAdvanced_(false),
	      correctVectorsAwesomePred_("0==1") // never correct
	{
		try {
			ioIn_.readline(correctVectorsAwesomePred_, "ChebyshevCorrectVector=");
		} catch (std::exception&) {}
	}

	virtual void calcTimeVectors(const VectorSizeType& indices,
	                             RealType Eg,
	                             const VectorWithOffsetType& phi,
	                             const typename BaseType::ExtraData& extra)
	{
		if (extra.wftAndAdvanceIfNeeded) {
			const SizeType noAdvance = indices[0];
			VectorWithOffsetType phiNew;
			if (targetVectors_[noAdvance].size() > 0) {
				BaseType::wftHelper().wftOneVector(phiNew,
				                                   targetVectors_[noAdvance],
				                                   extra.block[0]);

				targetVectors_[noAdvance] = phiNew;
			}
		}

		SizeType startOfWft = 1;
		if (currentTimeStep_ == 0) {
			SizeType indexOf1 = indices[startOfWft];
			assert(indexOf1 < targetVectors_.size());
			VectorWithOffsetType& tv1 =
			        const_cast<VectorWithOffsetType&>(targetVectors_[indexOf1]);
			tv1  = phi;
			startOfWft = 2;
		}

		// WFT 1 if !time advanced
		// WFT 2 if time advanced
		assert(0 < extra.block.size());
		SizeType n = indices.size();

		for (SizeType i = startOfWft; i < n; ++i) {
			SizeType ii = indices[i];
			BaseType::wftHelper().wftSome(targetVectors_, extra.block[0], ii, ii + 1);
		}

		assert(n > 0);
		if (currentTimeStep_ == 0 && tstStruct_.noOperator() && tstStruct_.skipTimeZero()) {
			for (SizeType i = 0; i < n; ++i) {
				SizeType ii = indices[i];
				targetVectors_[ii] = phi;
			}

			return;
		}

		const VectorWithOffsetType* ptr0 = &(targetVectors_[indices[0]]);
		const VectorWithOffsetType* ptr1 = &phi;
		if (ptr0 != ptr1)
			targetVectors_[indices[0]] = phi;

		if (times_.size() == 1 && fabs(times_[0])<1e-10) return;

		if (timeHasAdvanced_) {
			for (SizeType i = 0; i < n - 1; ++i) {
				SizeType ii = indices[i];
				SizeType jj = indices[i + 1];
				targetVectors_[ii] = targetVectors_[jj];
			}
		}

		for (SizeType i = 2; i < n; ++i) {
			SizeType ii = indices[i];
			assert(ii < targetVectors_.size());
			assert(ii != 1);
			targetVectors_[ii] = phi;
			SizeType prev = indices[i - 1];
			SizeType prevMinus2 = indices[i - 2];
			calcTargetVector(targetVectors_[ii], phi, prev, prevMinus2, Eg);
		}

		assert(extra.block.size() > 0);
		const SizeType site = extra.block[0];
		const SizeType nsites = model_.geometry().numberOfSites();
		const SizeType center = nsites/2;
		PredicateAwesomeType pred(correctVectorsAwesomePred_);
		const bool flagcorrection = pred.isTrue("%s", site, "%c", center, "%n", nsites);
		if (flagcorrection)
			correctVectors(indices, Eg);

		timeHasAdvanced_ = false;
	}

	void timeHasAdvanced()
	{
		timeHasAdvanced_ = true;
	}

	RealType time() const { return currentTimeStep_*tstStruct_.tau(); }

private:

	void calcTargetVector(VectorWithOffsetType& v,
	                      const VectorWithOffsetType& phi,
	                      SizeType prev,
	                      SizeType prevMinus2,
	                      RealType Eg)

	{
		for (SizeType ii=0;ii<phi.sectors();ii++) {
			SizeType i0 = phi.sector(ii);
			TargetVectorType r;
			calcTargetVector(r, phi, prev, prevMinus2, i0, Eg);
			v.setDataInSector(r,i0);
		}
	}

	void calcTargetVector(TargetVectorType& r,
	                      const VectorWithOffsetType& phi,
	                      SizeType prev,
	                      SizeType prevMinus2,
	                      SizeType i0,
	                      RealType Eg)
	{
		SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
		typename ModelHelperType::Aux aux(p, lrs_);
		typename ModelType::HamiltonianConnectionType hc(lrs_,
		                                                 model_.geometry(),
		                                                 ModelType::modelLinks(),
		                                                 time(),
		                                                 0);
		MatrixLanczosType lanczosHelper(model_, hc, aux);

		ProgramGlobals::VerboseEnum verbose = (model_.params().options.isSet("VerboseCheby"))
		        ? ProgramGlobals::VerboseEnum::YES
		        : ProgramGlobals::VerboseEnum::NO;

		// defining Hprime matrix:
		ScaledMatrixType lanczosHelper2(lanczosHelper, tstStruct_, Eg, verbose);

		SizeType total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		r.resize(total);
		if (currentTimeStep_ == 0) {
			phi.extract(phi2,i0);
			lanczosHelper2.matrixVectorProduct(r,phi2); // applying Hprime
		} else {
			TargetVectorType x2(total);
			VectorWithOffsetType x = 2.0*targetVectors_[prev];
			x.extract(x2,i0);
			targetVectors_[prevMinus2].extract(phi2,i0);
			lanczosHelper2.matrixVectorProduct(r,x2); // applying Hprime
			r += (-1.0)*phi2;
		}
	}

	void correctVectors(const VectorSizeType& indices, RealType Eg)
	{
		// take the first vector and compute V and weights
		const SizeType n = indices.size();
		if (n != 3) {
			PsimagLite::String msg = "TimeVectorsChebyshev:: correctVectors " +
			        PsimagLite::String(" indices.size() = " + ttos(indices.size()) +
			                           " != 3 (ignoring)\n");
			std::cout<<msg;
			std::cerr<<msg;
			return;
		}

		const SizeType indexToUse = 0;
		const VectorWithOffsetType& phi = targetVectors_[indices[indexToUse]];

		if (phi.sectors() == 0) {
			PsimagLite::String msg = "TimeVectorsChebyshev:: correctVectors " +
			        PsimagLite::String(" called too early maybe? (ignoring)\n");
			std::cout<<msg;
			std::cerr<<msg;
			return;
		}

		if (phi.sectors() != 1)
			err("Cheby: correctVectors: does NOT work for VectorWithOffsets\n");


		VectorMatrixFieldType V(phi.sectors());
		VectorType weights;

		VectorSizeType permutation;
		computeAuxForCorrection(V, weights, permutation, phi, phi.sector(0), Eg);
		if (weights.size() == 0) return;
		assert(V.size() == 1);

		// correct all vectors based on the same V and weights
		const SizeType nvectors = indices.size(); // = 3
		for (SizeType i = 0; i < nvectors; ++i)
			correctVectors(targetVectors_[indices[i]], V[0], weights, permutation);
	}

	void correctVectors(VectorWithOffsetType& phi,
	                    const MatrixRealType& Vmatrix,
	                    const VectorType& weights,
	                    const VectorSizeType& permutation)

	{
		for (SizeType ii = 0; ii < phi.sectors(); ++ii) {
			SizeType i0 = phi.sector(ii);
			TargetVectorType r;
			phi.extract(r, i0);
			correctVector(r, Vmatrix, weights, permutation);
			phi.setDataInSector(r, i0);
		}
	}

	void computeAuxForCorrection(VectorMatrixFieldType& V,
	                             VectorType& weights,
	                             VectorSizeType& permutation,
	                             const VectorWithOffsetType& phi,
	                             SizeType i0,
	                             RealType Eg) const
	{
		const SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
		typename ModelHelperType::Aux aux(p, lrs_);
		typename ModelType::HamiltonianConnectionType hc(lrs_,
		                                                 model_.geometry(),
		                                                 ModelType::modelLinks(),
		                                                 time(),
		                                                 0);
		MatrixLanczosType lanczosHelper(model_, hc, aux);
		ProgramGlobals::VerboseEnum verbose = (model_.params().options.isSet("VerboseCheby"))
		        ? ProgramGlobals::VerboseEnum::YES
		        : ProgramGlobals::VerboseEnum::NO;
		// defining Hprime matrix:
		ScaledMatrixType lanczosHelper2(lanczosHelper, tstStruct_, Eg, verbose);

		const RealType fakeTime = 0;
		typedef PsimagLite::NoPthreadsNg<ParallelTriDiagType> ParallelizerType;
		ParallelizerType threadedTriDiag(PsimagLite::CodeSectionParams(1));


		VectorMatrixFieldType T(phi.sectors());
		VectorSizeType steps(phi.sectors());

		ParallelTriDiagType helperTriDiag(phi,
		                                  T,
		                                  V,
		                                  steps,
		                                  lrs_,
		                                  fakeTime,
		                                  model_,
		                                  ioIn_);

		threadedTriDiag.loopCreate(helperTriDiag);

		VectorVectorRealType eigs(phi.sectors());

		for (SizeType ii = 0;ii < phi.sectors(); ++ii)
			PsimagLite::diag(T[ii], eigs[ii], 'V');

		VectorRealType& veigs = eigs[0];
		lanczosHelper2.scale(veigs);

		for (SizeType i = 0; i < veigs.size(); ++i)
			veigs[i] = -fabs(veigs[i]);

		const MatrixRealType& Vmatrix = V[0];

		const SizeType small = Vmatrix.cols();
		const SizeType big = Vmatrix.rows();
		PsimagLite::Sort<VectorRealType> sort;
		permutation.resize(small);
		sort.sort(veigs, permutation);

		SizeType counter = 0;
		for (; counter < small; ++counter)
			if (veigs[counter] > -1.0) break;

		if (counter == 0) return;

		const SizeType nbad = counter;
		weights.resize(nbad);

		VectorType sv;
		phi.extract(sv, i0);

		for (SizeType alpha = 0; alpha < nbad; ++alpha) {
			ComplexOrRealType sum = 0.0;
			for (SizeType j = 0; j < big; ++j)
				sum += PsimagLite::conj(Vmatrix(j, permutation[alpha]))*sv[j];
			weights[alpha] = sum;
		}
	}

	void correctVector(VectorType& r,
	                   const MatrixRealType& Vmatrix,
	                   const VectorType& weights,
	                   const VectorSizeType& permutation) const
	{

		const SizeType big = Vmatrix.rows();
		const SizeType nbad = weights.size();
		ComplexOrRealType sum = 0;
		ComplexOrRealType sumOld = 0;

		for (SizeType i = 0; i < big; ++i) {
			sumOld += r[i]*PsimagLite::conj(r[i]);
			for (SizeType alpha = 0; alpha < nbad; ++alpha)
				r[i] -= weights[alpha] * Vmatrix(i, permutation[alpha]);
			sum += r[i]*PsimagLite::conj(r[i]);
		}

		const RealType factor = sqrt(PsimagLite::real(sumOld))/sqrt(PsimagLite::real(sum));
		for (SizeType i = 0; i < big; ++i)
			r[i] *= factor;
	}

	const SizeType& currentTimeStep_;
	const TargetParamsType& tstStruct_;
	const VectorRealType& times_;
	VectorVectorWithOffsetType& targetVectors_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
	InputValidatorType& ioIn_;
	bool timeHasAdvanced_;
	PsimagLite::String correctVectorsAwesomePred_;
}; //class TimeVectorsChebyshev
} // namespace Dmrg
/*@}*/
#endif

