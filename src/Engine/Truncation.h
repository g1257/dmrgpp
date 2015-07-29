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

/*! \file Truncation.h
 *
 * Implements the truncation for the DMRG algorithm
 *
 */

#ifndef DMRG_TRUNCATION_H
#define DMRG_TRUNCATION_H

#include "DensityMatrix.h"
#include "Sort.h"
#include "Concurrency.h"

namespace Dmrg {

template<typename LeftRightSuperType,
         typename ParametersType,
         typename TargettingType>
class Truncation  {
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename LeftRightSuperType::ProgressIndicatorType ProgressIndicatorType;
	typedef typename TargettingType::RealType RealType;
	typedef typename TargettingType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef DensityMatrix<BasisType,BasisWithOperatorsType,TargettingType> DensityMatrixType;
	typedef typename TargettingType::ModelType ModelType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
		  EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM};

public:

	typedef typename LeftRightSuperType::SparseMatrixType TransformType;
	typedef typename DensityMatrixType::BlockMatrixType BlockMatrixType;

	struct TruncationCache {
		TruncationCache()
		    : transform(0,0)
		{}

		BlockMatrixType transform;
		typename PsimagLite::Vector<RealType>::Type eigs;
		typename PsimagLite::Vector<SizeType>::Type removedIndices;
	}; // TruncationCache

	Truncation(ReflectionSymmetryType& reflectionOperator,
	           WaveFunctionTransfType& waveFunctionTransformation,
	           const ParametersType& parameters,
	           bool verbose)
	    : reflectionOperator_(reflectionOperator),
	      lrs_(reflectionOperator_.leftRightSuper()),
	      waveFunctionTransformation_(waveFunctionTransformation),
	      parameters_(parameters),
	      verbose_(verbose),
	      progress_("Truncation"),
	      error_(0.0)
	{
		if (parameters_.tolerance<0) return;
		PsimagLite::OstringStream msg;
		msg<<"has tolerance= "<<parameters_.tolerance;
		progress_.printline(msg,std::cout);
	}

	void operator()(BasisWithOperatorsType& pS,
	                BasisWithOperatorsType& pE,
	                const TargettingType& target,
	                SizeType keptStates,
	                SizeType direction)
	{
		if (direction==EXPAND_SYSTEM) {
			progress_.print("for Environment\n",std::cout);
			changeBasis(pS,target,keptStates,direction);
			truncateBasisSystem(pS,lrs_.right());
		} else {
			progress_.print("for System\n",std::cout);
			changeBasis(pE,target,keptStates,direction);
			truncateBasisEnviron(pE,lrs_.left());
		}
	}

	const TransformType& transform() const
	{
		return ftransform_;
	}

	const RealType& error() const { return error_; }

	void changeBasis(BasisWithOperatorsType& sBasis,
	                 BasisWithOperatorsType& eBasis,
	                 const TargettingType& target,
	                 SizeType keptStates)
	{
		changeBasis(sBasis,target,keptStates,EXPAND_SYSTEM);
		changeBasis(eBasis,target,keptStates,EXPAND_ENVIRON);

		truncateBasisSystem(sBasis,lrs_.right());
		truncateBasisEnviron(eBasis,lrs_.left());

	}

private:

	void changeBasis(BasisWithOperatorsType& rSprime,
	                 const TargettingType& target,
	                 SizeType keptStates,
	                 SizeType direction)
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

		const BasisWithOperatorsType& pBasis = (direction==EXPAND_SYSTEM) ?
		            lrs_.left() : lrs_.right();
		const BasisWithOperatorsType& pBasisSummed = (direction==EXPAND_SYSTEM) ?
		            lrs_.right() : lrs_.left();

		DensityMatrixType dmS(target,pBasis,pBasisSummed,lrs_.super(),direction);
		dmS.check(direction);

		if (verbose_ && PsimagLite::Concurrency::root()) {
			std::cerr<<"Trying to diagonalize density-matrix with size=";
			std::cerr<<dmS.rank()<<"\n";
		}

		/* PSIDOC DiagOfDensityMatrix
			We then diagonalize $\hat{\rho}_S$, and obtain its eigenvalues and
		    eigenvectors,
			$w^S_{\alpha,\alpha'}$ in $\mathcal{V}(S')$ ordered in decreasing
		    eigenvalue order.
			We change basis for the operator $H^{S'}$ (and other operators as necessary),
		    as follows:
			\begin{equation}
			(H^{S' {\rm new\,\,basis}})_{\alpha,\alpha'}=(w^S)^{-1}_{\alpha,\gamma}
		    (H^{ S'})_{\gamma,\gamma'}w^S_{\gamma',\alpha'}.
			\label{eq:transformation}
			\end{equation}
		*/

		TruncationCache& cache = (direction==EXPAND_SYSTEM) ? leftCache_ : rightCache_;

		dmS.diag(cache.eigs,'V');
		dmS.check2(direction);
		updateKeptStates(keptStates,cache.eigs);

		//! transform basis: dmS^\dagger * operator matrix * dms
		cache.transform = dmS();
		if (parameters_.options.find("nodmrgtransform") != PsimagLite::String::npos) {
			cache.transform.setDiagonal();
		}

		rSprime = pBasis;
		rSprime.changeBasis(cache.removedIndices,cache.eigs,keptStates,parameters_);

		PsimagLite::OstringStream msg2;
		msg2<<"done with entanglement";
		progress_.printline(msg2,std::cout);

	}

	void truncateBasisSystem(BasisWithOperatorsType& rSprime,
	                         const BasisWithOperatorsType& eBasis)
	{

		PsimagLite::OstringStream msg;
		TruncationCache& cache = leftCache_;

		rSprime.truncateBasis(ftransform_,cache.transform,
		                      cache.eigs,cache.removedIndices);
		LeftRightSuperType lrs(rSprime,(BasisWithOperatorsType&) eBasis,
		                       (BasisType&)lrs_.super());
		bool twoSiteDmrg = (parameters_.options.find("twositedmrg")!=PsimagLite::String::npos);
		const LeftRightSuperType& lrsForWft = (twoSiteDmrg) ? lrs_ : lrs;
		waveFunctionTransformation_.push(ftransform_,EXPAND_SYSTEM,lrsForWft);

		msg<<"new size of basis="<<rSprime.size();
		progress_.printline(msg,std::cout);
	}

	void truncateBasisEnviron(BasisWithOperatorsType& rEprime,
	                          const BasisWithOperatorsType& sBasis)
	{

		PsimagLite::OstringStream msg;
		TruncationCache& cache = rightCache_;

		rEprime.truncateBasis(ftransform_,cache.transform,
		                      cache.eigs,cache.removedIndices);
		LeftRightSuperType lrs((BasisWithOperatorsType&) sBasis,
		                       rEprime,(BasisType&)lrs_.super());
		bool twoSiteDmrg = (parameters_.options.find("twositedmrg")!=PsimagLite::String::npos);
		const LeftRightSuperType& lrsForWft = (twoSiteDmrg) ? lrs_ : lrs;
		waveFunctionTransformation_.push(ftransform_,EXPAND_ENVIRON,lrsForWft);
		msg<<"new size of basis="<<rEprime.size();
		progress_.printline(msg,std::cout);
	}

	void updateKeptStates(SizeType& keptStates,
	                      const typename PsimagLite::Vector<RealType>::Type& eigs2)
	{
		typename PsimagLite::Vector<RealType>::Type eigs = eigs2;
		typename PsimagLite::Vector<SizeType>::Type perm(eigs.size());
		PsimagLite::Sort<typename PsimagLite::Vector<RealType>::Type> sort;
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
		error_ = discWeight;
		progress_.printline(msg,std::cout);
		// we report the discarded weight
		msg<<"Discarded weight (Truncation error): "<< discWeight ;
		progress_.printline(msg,std::cout);
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
	                           const typename PsimagLite::Vector<RealType>::Type& eigs) const
	{
		if (parameters_.tolerance<0) return keptStates;
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
			if (discWeight>parameters_.tolerance) {
				total = eigs.size() - i;
				discWeight -= fabs(eigs[i]);
				break;
			}
		}

		// if total is too small or too big we keep it unchanged
		if (total>=keptStates || total<parameters_.keptStatesInfinite)
			return keptStates;

		return total;
	}

	RealType sumUpTo(const typename PsimagLite::Vector<RealType>::Type& eigs,
	                 SizeType x) const
	{
		RealType discWeight = 0;
		for (SizeType i=0;i<x;i++)
			discWeight += fabs(eigs[i]);
		return discWeight;
	}

	void dumpEigs(const typename PsimagLite::Vector<RealType>::Type& eigs) const
	{
		static SizeType counter = 0;
		if (parameters_.fileForDensityMatrixEigs=="") return;
		PsimagLite::String file(parameters_.fileForDensityMatrixEigs);
		file += ttos(counter);
		PsimagLite::IoSimple::Out io(file);
		io<<eigs;
		counter++;
	}

	ReflectionSymmetryType& reflectionOperator_;
	const LeftRightSuperType& lrs_;
	WaveFunctionTransfType& waveFunctionTransformation_;
	const ParametersType& parameters_;
	bool verbose_;
	ProgressIndicatorType progress_;
	RealType error_;
	TransformType ftransform_;
	TruncationCache leftCache_,rightCache_;

}; // class Truncation

} // namespace
/*@}*/
#endif // DMRG_TRUNCATION_H

