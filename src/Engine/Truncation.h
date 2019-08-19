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
#include "DensityMatrixSu2.h"
#include "Sort.h"
#include "Concurrency.h"
#include "Io/IoNg.h"
#include "Profiling.h"

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
	typedef DensityMatrixSu2<TargetingType> DensityMatrixSu2Type;
	typedef DensityMatrixBase<TargetingType> DensityMatrixBaseType;
	typedef typename DensityMatrixBaseType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename TargetingType::ModelType ModelType;
	typedef typename ModelType::GeometryType GeometryType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
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

	Truncation(ReflectionSymmetryType& reflectionOperator,
	           WaveFunctionTransfType& waveFunctionTransformation,
	           const ParametersType& parameters,
	           const GeometryType& geometry,
	           IoOutType& ioOut)
	    : reflectionOperator_(reflectionOperator),
	      lrs_(reflectionOperator_.leftRightSuper()),
	      waveFunctionTransformation_(waveFunctionTransformation),
	      parameters_(parameters),
	      geometry_(geometry),
	      ioOut_(ioOut),
	      progress_("Truncation"),
	      error_(0.0)
	{
		if (parameters_.truncationControl.first < 0) return;
		PsimagLite::OstringStream msg;
		msg<<"has tolerance= "<<parameters_.truncationControl.first;
		msg<<" minimum m= "<<parameters_.truncationControl.second;
		progress_.printline(msg,std::cout);
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
		            keptStates, ProgramGlobals::DirectionEnum::EXPAND_SYSTEM,
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
		bool useSvd = (parameters_.options.find("truncationNoSvd") == PsimagLite::String::npos);
		bool enablePersistentSvd = (parameters_.options.find("EnablePersistentSvd") !=
		        PsimagLite::String::npos);
		ParamsDensityMatrixType p(useSvd, direction, debug, enablePersistentSvd);
		TruncationCache& cache = (direction == expandSys) ? leftCache_ :
		                                                    rightCache_;

		if (BasisType::useSu2Symmetry()) {
			if (p.useSvd) {
				std::cerr<<"WARNING: SVD for truncation NOT supported with SU(2)\n";
				p.useSvd = false;
			}

			*dm = new DensityMatrixSu2Type(target,lrs_,p);
		} else if (p.useSvd) {
			*dm = new DensityMatrixSvdType(target,lrs_,p);
		} else {
			*dm = new DensityMatrixLocalType(target,lrs_,p);
		}

		assert(*dm);
		DensityMatrixBaseType* dmS = *dm;
		assert(dmS);

		dmS->diag(cache.eigs,'V');

		updateKeptStates(keptStates, cache.eigs);

		cache.transform = dmS->operator()();
		if (parameters_.options.find("nodmrgtransform") != PsimagLite::String::npos) {
			PsimagLite::OstringStream msg;
			msg<<"SolverOptions=nodmrgtransform, setting transform to identity";
			progress_.printline(msg,std::cout);
			cache.transform.setTo(1.0);
		}

		rSprime = pBasis;
		rSprime.changeBasis(cache.removedIndices,cache.eigs,keptStates,parameters_);
	}

	void truncateBasis(BasisWithOperatorsType& rPrime,
	                   const BasisWithOperatorsType& oppoBasis,
	                   const DensityMatrixBaseType& dms,
	                   ProgramGlobals::DirectionEnum direction)
	{
		bool expandSys = (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		const BasisWithOperatorsType& basis = (expandSys) ? lrs_.left() : lrs_.right();
		SizeType site = 0; // FIXME for model Immm
		size_t mostRecent = lrs_.left().operatorsPerSite(site)*geometry_.maxConnections();
		size_t numOfOp = basis.numberOfOperators();
		PairSizeSizeType startEnd(0, numOfOp);
		if (startEnd.second > mostRecent) {
			if (expandSys) startEnd.first = startEnd.second - mostRecent;
			else startEnd.second = mostRecent;
		}

		PsimagLite::OstringStream msg;
		TruncationCache& cache = (expandSys) ? leftCache_ : rightCache_;

		cache.transform.truncate(cache.removedIndices);

		const bool blasIsThreadSafe = (parameters_.options.find("blasNotThreadSafe") ==
		        PsimagLite::String::npos);
		rPrime.truncateBasis(cache.transform,
		                     cache.eigs,
		                     cache.removedIndices,
		                     startEnd,
		                     blasIsThreadSafe);
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
		progress_.printline(msg,std::cout);

		delete lrs;
		lrs = 0;
	}

	void updateKeptStates(SizeType& keptStates,
	                      const VectorRealType& eigs2)
	{
		VectorRealType eigs = eigs2;
		typename PsimagLite::Vector<SizeType>::Type perm(eigs.size());
		PsimagLite::Sort<VectorRealType> sort;
		sort.sort(eigs,perm);
		dumpEigs(eigs);

		SizeType newKeptStates = computeKeptStates(keptStates,eigs);
		SizeType statesToRemove = 0;
		if (eigs.size()>=newKeptStates)
			statesToRemove = eigs.size()-newKeptStates;
		RealType discWeight = sumUpTo(eigs,statesToRemove);
		PsimagLite::OstringStream msg;
		if (newKeptStates != keptStates) {
			// we report that the "m" value has been changed and...
			msg<<"Reducing kept states to "<<newKeptStates<<" from "<<keptStates;
			// ... we change it:
			keptStates = newKeptStates;
		} else {
			// we report that the "m" value remains the same
			msg<<"Not changing kept states="<<keptStates;
		}

		progress_.printline(msg,std::cout);

		error_ = discWeight;

		// we report the discarded weight
		PsimagLite::OstringStream msg2;
		msg2<<"Discarded weight (Truncation error): "<< discWeight;
		progress_.printline(msg2,std::cout);

		if (parameters_.options.find("calcAndPrintEntropies") != PsimagLite::String::npos)
			calcAndPrintEntropies(eigs);
	}

	void calcAndPrintEntropies(const VectorRealType& eigs)
	{
		RealType rntentropy = entropy(eigs, 1.0);
		const RealType reyniIndex = 2.0;
		RealType r2p0entropy = entropy(eigs, reyniIndex);
		PsimagLite::OstringStream msg;
		// von-neumann entaglement entropy; and 2nd order Reyni entropy
		msg<<"EntropyVonNeumann= "<<rntentropy;
		msg<<"; ReyniIndex= "<<reyniIndex<<" EntropyReyni="<<r2p0entropy;
		progress_.printline(msg, std::cout);
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
		for (SizeType i=0;i<x;i++)
			discWeight += fabs(eigs[i]);
		return discWeight;
	}

	void dumpEigs(const VectorRealType& eigs)
	{
		if (parameters_.options.find("saveDensityMatrixEigenvalues") == PsimagLite::String::npos)
			return;

		PsimagLite::String label("DensityMatrixEigenvalues");

		static bool firstCall = true;
		if (firstCall) {
			ioOut_.createGroup(label);
			SizeType n = geometry_.numberOfSites();
			counterVector_.resize(n, 0);
			ioOut_.write(n, label + "/Size");
			for (SizeType i = 0; i < n; ++i)
				ioOut_.createGroup(label + "/" + ttos(i));
			firstCall = false;
		}

		SizeType last = lrs_.left().block().size();
		assert(lrs_.left().block().size() > 0);
		SizeType index = lrs_.left().block()[last - 1];

		assert(index < counterVector_.size());
		SizeType counter = counterVector_[index];

		ioOut_.write(eigs, label + "/" + ttos(index) + "/" + ttos(counter));


		ioOut_.write(counter + 1,
		             label + "/" + ttos(index) + "/Size",
		             (counter == 0) ?  IoOutType::Serializer::NO_OVERWRITE :
		                               IoOutType::Serializer::ALLOW_OVERWRITE);

		++counterVector_[index];
	}

	ReflectionSymmetryType& reflectionOperator_;
	const LeftRightSuperType& lrs_;
	WaveFunctionTransfType& waveFunctionTransformation_;
	const ParametersType& parameters_;
	const GeometryType& geometry_;
	IoOutType& ioOut_;
	ProgressIndicatorType progress_;
	RealType error_;
	TruncationCache leftCache_;
	TruncationCache rightCache_;
	VectorSizeType counterVector_;
}; // class Truncation

} // namespace
/*@}*/
#endif // DMRG_TRUNCATION_H

