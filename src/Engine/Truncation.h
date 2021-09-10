/*
Copyright (c) 2009-2014-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

/*! \file Truncation.h
 *
 * Implements the truncation for the DMRG algorithm
 *
 */

#ifndef DMRG_TRUNCATION_H
#define DMRG_TRUNCATION_H

#include "DensityMatrixLocal.h"
#include "DensityMatrixSvd.h"
#include "Sort.h"
#include "Concurrency.h"
#include "Io/IoNg.h"
#include "Profiling.h"
#include "PredicateAwesome.h"
#include "OutputFileOrNot.h"

namespace Dmrg {

template<typename ParametersType,
         typename TargetingType>
class Truncation  {

	typedef typename TargetingType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType VectorSizeType;
	typedef typename BasisWithOperatorsType::PairSizeSizeType PairSizeSizeType;
	typedef typename LeftRightSuperType::ProgressIndicatorType ProgressIndicatorType;
	typedef typename TargetingType::SparseMatrixType SparseMatrixType;
	typedef typename TargetingType::RealType RealType;
	typedef typename TargetingType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef DensityMatrixLocal<TargetingType> DensityMatrixLocalType;
	typedef DensityMatrixSvd<TargetingType> DensityMatrixSvdType;
	typedef DensityMatrixBase<TargetingType> DensityMatrixBaseType;
	typedef typename DensityMatrixBaseType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename TargetingType::ModelType ModelType;
	typedef typename ModelType::SuperGeometryType SuperGeometryType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	typedef typename DensityMatrixBaseType::Params ParamsDensityMatrixType;
	typedef BlockDiagonalMatrixType TransformType;
	typedef PsimagLite::IoNg::Out IoOutType;

	struct TruncationCache {

		BlockDiagonalMatrixType transform;
		VectorRealType eigs;
		typename PsimagLite::Vector<SizeType>::Type removedIndices;

	}; // TruncationCache

	Truncation(const LeftRightSuperType& lrs,
	           WaveFunctionTransfType& waveFunctionTransformation,
	           const ParametersType& parameters,
	           const SuperGeometryType& geometry,
	           OutputFileOrNot& ioOut)
	    : lrs_(lrs),
	      waveFunctionTransformation_(waveFunctionTransformation),
	      parameters_(parameters),
	      superGeometry_(geometry),
	      ioOut_(ioOut),
	      progress_("Truncation"),
	      error_(0.0)
	{
		firstCall_ = true;
		if (parameters_.truncationControl.first < 0) return;
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"has tolerance= "<<parameters_.truncationControl.first;
		msg<<" minimum m= "<<parameters_.truncationControl.second;
		progress_.printline(msgg, std::cout);
	}

	void changeBasisFinite(BasisWithOperatorsType& pS,
	                       BasisWithOperatorsType& pE,
	                       const TargetingType& target,
	                       SizeType keptStates,
	                       ProgramGlobals::DirectionEnum direction)
	{
		PsimagLite::Profiling profiling("TruncationChangeBasis", std::cout);
		DensityMatrixBaseType* dmS = 0;

		if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			changeBasis(pS,target,keptStates,direction, &dmS);
			assert(dmS);
			truncateBasis(pS,lrs_.right(), *dmS, direction);
		} else {
			changeBasis(pE,target,keptStates,direction, &dmS);
			assert(dmS);
			truncateBasis(pE,lrs_.left(), *dmS, direction);
		}

		delete dmS;
		dmS = 0;
	}

	const TransformType& transform(ProgramGlobals::DirectionEnum direction) const
	{
		return (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? leftCache_.transform
		                                                                   : rightCache_.transform;
	}

	const RealType& error() const { return error_; }

	void changeBasisInfinite(BasisWithOperatorsType& sBasis,
	                         BasisWithOperatorsType& eBasis,
	                         const TargetingType& target,
	                         SizeType keptStates)
	{
		PsimagLite::Profiling profiling("TruncationChangeBasis", std::cout);

		DensityMatrixBaseType* dmS = 0;
		changeBasis(sBasis,
		            target,
		            keptStates,
		            ProgramGlobals::DirectionEnum::EXPAND_SYSTEM,
		            &dmS);
		assert(dmS);
		truncateBasis(sBasis,
		              lrs_.right(),
		              *dmS,
		              ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		delete dmS;
		dmS = 0;

		DensityMatrixBaseType* dmE = 0;
		changeBasis(eBasis,
		            target,
		            keptStates,
		            ProgramGlobals::DirectionEnum::EXPAND_ENVIRON,
		            &dmE);
		assert(dmE);
		truncateBasis(eBasis,
		              lrs_.left(),
		              *dmE,
		              ProgramGlobals::DirectionEnum::EXPAND_ENVIRON);
		delete dmE;
		dmE = 0;
	}

private:

	void changeBasis(BasisWithOperatorsType& rSprime,
	                 const TargetingType& target,
	                 SizeType keptStates,
	                 ProgramGlobals::DirectionEnum direction,
	                 DensityMatrixBaseType** dm)
	{
		/* PSIDOC Truncation
			Let us define the density matrices for system:
			\begin{equation}
			(\hat{\rho}_S)_{\alpha,\alpha'} = \sum_{\beta\in\mathcal{V}(E')}
			\psi_{\alpha',\beta}^*\psi_{\alpha,\beta}
			\label{eq:rhoSystem}
			\end{equation}
			in $\mathcal{V}(S')$,
			and environment:
			\begin{equation}
			(\hat{\rho}_E )_{\beta,\beta'}= \sum_{\alpha\in \mathcal{V}(S')}
			\psi_{\alpha,\beta'}^*\psi_{\alpha,\beta}
			\label{eq:rhoEnviron}
			\end{equation}
			in $\mathcal{V}(E')$.
			*/

		const ProgramGlobals::DirectionEnum expandSys = ProgramGlobals::DirectionEnum::EXPAND_SYSTEM;
		const BasisWithOperatorsType& pBasis = (direction == expandSys) ?
		            lrs_.left() : lrs_.right();

		bool debug = false;
		bool useSvd = !parameters_.options.isSet("truncationNoSvd");
		bool enablePersistentSvd = parameters_.options.isSet("EnablePersistentSvd");
		bool serialSvd = parameters_.options.isSet("SerialSvd");
		ParamsDensityMatrixType p(useSvd, direction, debug, enablePersistentSvd, serialSvd);
		TruncationCache& cache = (direction == expandSys) ? leftCache_ :
		                                                    rightCache_;

		if (BasisType::useSu2Symmetry()) {
			err("SU(2) no longer supported\n");
		} else if (p.useSvd) {
			*dm = new DensityMatrixSvdType(target,lrs_,p);
		} else {
			*dm = new DensityMatrixLocalType(target,lrs_,p);
		}

		assert(*dm);
		DensityMatrixBaseType* dmS = *dm;
		assert(dmS);

		dmS->diag(cache.eigs,'V');

		typename PsimagLite::Vector<SizeType>::Type perm(cache.eigs.size());
		PsimagLite::Sort<VectorRealType> sort;
		sort.sort(cache.eigs, perm);

		printSumAndCheckEigs(cache.eigs);

		updateKeptStates(keptStates, cache.eigs);

		cache.transform = dmS->operator()();
		if (parameters_.options.isSet("nodmrgtransform")) {
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"SolverOptions=nodmrgtransform, setting transform to identity";
			progress_.printline(msgg, std::cout);
			cache.transform.setTo(1.0);
		}

		rSprime = pBasis;
		rSprime.changeBasis(cache.removedIndices, perm, keptStates, parameters_);
	}

	void truncateBasis(BasisWithOperatorsType& rPrime,
	                   const BasisWithOperatorsType& oppoBasis,
	                   const DensityMatrixBaseType& dms,
	                   ProgramGlobals::DirectionEnum direction)
	{
		bool expandSys = (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		const BasisWithOperatorsType& basis = (expandSys) ? lrs_.left() : lrs_.right();
		size_t mostRecent = superGeometry_.hollowOutRadius( maxOpsPerSiteLeft());
		size_t numOfOp = basis.numberOfLocalOperators();
		PairSizeSizeType startEnd(0, numOfOp);
		if (startEnd.second > mostRecent) {
			if (expandSys) startEnd.first = startEnd.second - mostRecent;
			else startEnd.second = mostRecent;
		}

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		TruncationCache& cache = (expandSys) ? leftCache_ : rightCache_;

		cache.transform.truncate(cache.removedIndices);

		rPrime.truncateBasis(cache.transform,
		                     cache.eigs,
		                     cache.removedIndices,
		                     startEnd,
		                     parameters_.gemmRnb,
		                     PsimagLite::Concurrency::codeSectionParams.npthreadsLevelTwo,
		                     parameters_.opOnSiteThreshold);

		LeftRightSuperType* lrs = 0;
		if (expandSys)
			lrs = new LeftRightSuperType(rPrime,
			                             const_cast<BasisWithOperatorsType&>(oppoBasis),
			                             const_cast<BasisType&>(lrs_.super()));
		else
			lrs = new LeftRightSuperType(const_cast<BasisWithOperatorsType&>(oppoBasis),
			                             rPrime,
			                             const_cast<BasisType&>(lrs_.super()));

		bool twoSiteDmrg = waveFunctionTransformation_.options().twoSiteDmrg;
		bool wftInPatches = (waveFunctionTransformation_.options().accel ==
		                     WaveFunctionTransfType::WftOptionsType::ACCEL_PATCHES);
		const LeftRightSuperType& lrsForWft = (twoSiteDmrg || wftInPatches) ? lrs_ : *lrs;
		waveFunctionTransformation_.push(cache.transform,
		                                 direction,
		                                 lrsForWft,
		                                 dms.vts(),
		                                 dms.s(),
		                                 dms.qns());

		msg<<"new size of basis="<<rPrime.size();
		msg<<" transform is "<<cache.transform.rows()<<" x "<<cache.transform.cols();
		msg<<" with "<<cache.transform.blocks()<<" symmetry blocks";
		progress_.printline(msgg, std::cout);

		delete lrs;
		lrs = 0;
	}

	// there is no right
	SizeType maxOpsPerSiteLeft() const
	{
		const SizeType n = lrs_.left().block().size();
		SizeType max = lrs_.left().operatorsPerSite(0);
		for (SizeType i = 1; i < n; ++i) {
			if (max < lrs_.left().operatorsPerSite(i))
				max = lrs_.left().operatorsPerSite(i);
		}

		return max;
	}

	void printSumAndCheckEigs(const VectorRealType& eigs) const
	{
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();

		RealType sum = checkAndSum(eigs);

		msg<<"Sum of Singular Values of the Density Matrix: "<<sum;
		progress_.printline(msgg, std::cout);

		if (fabs(sum - 1) > 1e-3)
			err("printSumAndCheckEigs: DM eigs don't amount to one\n");
	}

	void updateKeptStates(SizeType& keptStates,
	                      const VectorRealType& eigs)
	{
		dumpEigs(eigs);

		SizeType newKeptStates = computeKeptStates(keptStates, eigs);
		SizeType statesToRemove = 0;
		if (eigs.size()>=newKeptStates)
			statesToRemove = eigs.size() - newKeptStates;
		RealType discWeight = sumUpTo(eigs, statesToRemove);
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		if (newKeptStates != keptStates) {
			// we report that the "m" value has been changed and...
			msg<<"Reducing kept states to "<<newKeptStates<<" from "<<keptStates;
			// ... we change it:
			keptStates = newKeptStates;
		} else {
			// we report that the "m" value remains the same
			msg<<"Not changing kept states="<<keptStates;
		}

		progress_.printline(msgg, std::cout);

		error_ = discWeight;

		// we report the discarded weight
		PsimagLite::OstringStream msgg2(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
		msg2<<"Discarded weight (Truncation error): "<< discWeight;
		progress_.printline(msgg2, std::cout);

		if (parameters_.options.isSet("calcAndPrintEntropies"))
			calcAndPrintEntropies(eigs);
	}

	void calcAndPrintEntropies(const VectorRealType& eigs)
	{
		RealType rntentropy = entropy(eigs, 1.0);
		const RealType reyniIndex = 2.0;
		RealType r2p0entropy = entropy(eigs, reyniIndex);
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		// von-neumann entaglement entropy; and 2nd order Reyni entropy
		msg<<"EntropyVonNeumann= "<<rntentropy;
		msg<<"; ReyniIndex= "<<reyniIndex<<" EntropyReyni="<<r2p0entropy;
		progress_.printline(msgg, std::cout);
	}

	RealType entropy(const VectorRealType& eigs,
	                 const RealType reyniIndex) const
	{
		RealType ent = 0;
		RealType val = 0;

		// von-neumann entanglement entropy
		if (reyniIndex == 1.0) {
			for (SizeType i = 0; i < eigs.size(); ++i) {
				val = eigs[i];
				if (val <= 1e-12) continue;
				ent += -1.0*val*log(val);
			}

			return ent;
		}

		// Reyni entropy of index n
		val = 0.0;
		for (SizeType i = 0;i < eigs.size(); ++i) {
			if (eigs[i] <= 1e-12) continue;
			val += pow(eigs[i], reyniIndex);
		}

		ent = 1.0/(1.0 - reyniIndex);
		return ent*log(val);
	}

	/* PSIDOC RemovalOfStates
		Let $m_S$ (here given by \verb!keptStates_! be a fixed number that
		corresponds to the number of states in $\mathcal{V}(S')$ that we want to keep.
		Consider the first $m_S$ eigenvectors $w^S$,
		 and let us call the Hilbert space spanned by them, $\mathcal{V}_R(S')$,
		 the DMRG-reduced Hilbert space on
		block $S'$. If $m_S\ge\#\mathcal{V}(S')$ then we keep all eigenvectors
		and there is effectively no truncation.
		We truncate the matrices $(H^{S' {\rm new\,\,basis}})$
		(and other operators as necessary)
		such that they now act on this truncated Hilbert space, $\mathcal{V}_R(S')$.
		We proceed in the same manner for the environment.
		!PTEX-END */
	//! eigenvalues are ordered in increasing order
	SizeType computeKeptStates(SizeType& keptStates,
	                           const VectorRealType& eigs) const
	{
		if (parameters_.truncationControl.first < 0) return keptStates;
		int start = eigs.size() - keptStates;
		if (start<0) start = 0;
		int maxToRemove = eigs.size()-parameters_.keptStatesInfinite;
		if (maxToRemove<0) maxToRemove = 0;
		SizeType total = parameters_.keptStatesInfinite;
		RealType discWeight=sumUpTo(eigs,start);
		// maybe we should use int instead of SizeType here!!!

		for (int i=start;i<maxToRemove;i++) {
			// calculate the discarded weight if we keep i states.
			discWeight += fabs(eigs[i]);
			// if the discarded weight
			// gets larger than the tolerance, we break the loop.
			if (discWeight > parameters_.truncationControl.first) {
				total = eigs.size() - i;
				discWeight -= fabs(eigs[i]);
				break;
			}
		}

		// if total is too big we keep it unchanged
		if (total>=keptStates)
			return keptStates;

		if (total < parameters_.truncationControl.second)
			return parameters_.truncationControl.second;

		return total;
	}

	RealType sumUpTo(const VectorRealType& eigs,
	                 SizeType x) const
	{
		RealType discWeight = 0;
		for (SizeType i = 0; i < x; ++i)
			discWeight += fabs(eigs[i]);
		return discWeight;
	}

	RealType checkAndSum(const VectorRealType& eigs) const
	{
		SizeType x = eigs.size();
		RealType sum = 0;
		for (SizeType i = 0; i < x; ++i) {
			const RealType val = eigs[i];
			if (val < 0)
				err("checkAndSum: Density Matrix eigenvalue is less than zero\n");
			if (val > 1)
				err("checkAndSum: Density Matrix eigenvalue is greater than one\n");
			sum += eigs[i];
		}

		return sum;
	}

	void dumpEigs(const VectorRealType& eigs)
	{
		SizeType last = lrs_.left().block().size();
		assert(lrs_.left().block().size() > 0);
		SizeType index = lrs_.left().block()[last - 1];

		PsimagLite::String predicate = parameters_.saveDensityMatrixEigenvalues;
		const SizeType center = superGeometry_.numberOfSites()/2;
		PsimagLite::PredicateAwesome<>::replaceAll(predicate, "c", ttos(center));
		PsimagLite::PredicateAwesome<> pAwesome(predicate);

		if (!pAwesome.isTrue("s", index))
			return;

		PsimagLite::String label("DensityMatrixEigenvalues");

		if (firstCall_) {
			ioOut_.createGroup(label);
			SizeType n = superGeometry_.numberOfSites();
			counterVector_.resize(n, 0);
			ioOut_.write(n, label + "/Size");
			for (SizeType i = 0; i < n; ++i)
				ioOut_.createGroup(label + "/" + ttos(i));
			firstCall_ = false;
		}

		assert(index < counterVector_.size());
		SizeType counter = counterVector_[index];

		ioOut_.write(eigs, label + "/" + ttos(index) + "/" + ttos(counter));


		ioOut_.write(counter + 1,
		             label + "/" + ttos(index) + "/Size",
		             (counter == 0) ?  IoOutType::Serializer::NO_OVERWRITE :
		                               IoOutType::Serializer::ALLOW_OVERWRITE);

		++counterVector_[index];
	}

	const LeftRightSuperType& lrs_;
	WaveFunctionTransfType& waveFunctionTransformation_;
	const ParametersType& parameters_;
	const SuperGeometryType& superGeometry_;
	OutputFileOrNot& ioOut_;
	ProgressIndicatorType progress_;
	RealType error_;
	TruncationCache leftCache_;
	TruncationCache rightCache_;
	VectorSizeType counterVector_;
	static bool firstCall_;
}; // class Truncation

template<typename T1, typename T2>
bool Truncation<T1, T2>::firstCall_ = true;

} // namespace
/*@}*/
#endif // DMRG_TRUNCATION_H

