
/*
Copyright (c) 2009-2011, UT-Battelle, LLC
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
// END LICENSE BLOCK
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
		typedef typename TargettingType::ConcurrencyType ConcurrencyType;
		typedef DensityMatrix<RealType,BasisType,BasisWithOperatorsType,TargettingType> DensityMatrixType;
		typedef typename TargettingType::ModelType ModelType;
		typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;

		enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
		EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM};
		
	public:

		typedef typename DensityMatrixType::BuildingBlockType TransformType;

		Truncation(ReflectionSymmetryType& reflectionOperator,
		           WaveFunctionTransfType& waveFunctionTransformation,
		           ConcurrencyType& concurrency,
		           const ParametersType& parameters,
		           size_t maxConnections,
		           bool verbose)
		: reflectionOperator_(reflectionOperator),
		  lrs_(reflectionOperator_.leftRightSuper()),
		  waveFunctionTransformation_(waveFunctionTransformation),
		  concurrency_(concurrency),
		  parameters_(parameters),
		  maxConnections_(maxConnections),
		  verbose_(verbose),
		  progress_("Truncation",0),
		  keptStates_(0),
		  error_(0.0)
		{
			if (parameters_.tolerance<0) return;
			std::ostringstream msg;
			msg<<"has tolerance= "<<parameters_.tolerance;
			progress_.printline(msg,std::cout);
		}

		void operator()(BasisWithOperatorsType& pS,
		                BasisWithOperatorsType& pE,
		                const TargettingType& target,
		                size_t keptStates,
		                size_t direction)
		{
// 			TransformType transform;
			if (direction==EXPAND_SYSTEM) {
				progress_.print("for Environment\n",std::cout);
				this->operator()(pS,target,keptStates,direction);
			} else {
				progress_.print("for System\n",std::cout);
				this->operator()(pE,target,keptStates,direction);
				
			}
		}

		//! Truncate basis
		void operator()(BasisWithOperatorsType &rSprime,
						const TargettingType& target,
				  size_t keptStates,
				  size_t direction) 
		{
			keptStates_ = keptStates;
			//size_t opPerSite = lrs_.left().numberOfOperatorsPerSite();
// 			size_t mostRecent = lrs_.left().numberOfOperatorsPerSite() * maxConnections_;
			if (direction==EXPAND_SYSTEM) {
// 				size_t numOfOp = lrs_.left().numberOfOperators();
// 				std::pair<size_t,size_t> startEnd(0,numOfOp);
// 				if (startEnd.second>mostRecent) // && numOfOp>4*opPerSite)
// 					startEnd.first = startEnd.second - mostRecent;
				changeAndTruncateBasis(rSprime,target,lrs_.left(),lrs_.right(),
									   direction); //startEnd);
			} else {
// 				size_t numOfOp = lrs_.right().numberOfOperators();
// 				std::pair<size_t,size_t> startEnd(0,numOfOp);
// 				if (startEnd.second>mostRecent) // && numOfOp>4*opPerSite)
// 					startEnd.second = mostRecent;
				changeAndTruncateBasis(rSprime,target,lrs_.right(),lrs_.left(),
									   direction); //startEnd);
			}
		}

		const TransformType& transform() const { return ftransform_; }
		
		const RealType& error() const { return error_; }

	private:

		void changeAndTruncateBasis(BasisWithOperatorsType &rSprime,
		                            const TargettingType& target,
		                            BasisWithOperatorsType const &pBasis,
		                            BasisWithOperatorsType const &pBasisSummed,
		                            size_t direction)
// 		                            const std::pair<size_t,size_t>& startEnd)
		{
			/** !PTEX-START Truncation
			Let us define the density matrices for system:
			\begin{equation}
			(\hat{\rho}_S)_{\alpha,\alpha'} = \sum_{\beta\in\mathcal{V}(E')}\psi_{\alpha',\beta}^*\psi_{\alpha,\beta}
			\label{eq:rhoSystem}
			\end{equation}
			in $\mathcal{V}(S')$,
			and environment:
			\begin{equation}
			(\hat{\rho}_E )_{\beta,\beta'}= \sum_{\alpha\in \mathcal{V}(S')}\psi_{\alpha,\beta'}^*\psi_{\alpha,\beta}
			\label{eq:rhoEnviron}
			\end{equation}
			in $\mathcal{V}(E')$.
			!PTEX-END */
			DensityMatrixType dmS(target,pBasis,pBasisSummed,lrs_.super(),direction);
			dmS.check(direction);
			
			if (verbose_ && concurrency_.root()) {
				std::cerr<<"Trying to diagonalize density-matrix with size=";
				std::cerr<<dmS.rank()<<"\n";
			}
			std::vector<RealType> eigs;
			/** !PTEX-START DiagOfDensityMatrix
			We then diagonalize $\hat{\rho}_S$, and obtain its eigenvalues and eigenvectors, 
			$w^S_{\alpha,\alpha'}$ in $\mathcal{V}(S')$ ordered in decreasing eigenvalue order.
			We change basis for the operator $H^{S'}$ (and other operators as necessary), as follows:
			\begin{equation}
			(H^{S' {\rm new\,\,basis}})_{\alpha,\alpha'}=(w^S)^{-1}_{\alpha,\gamma} (H^{ S'})_{\gamma,\gamma'}w^S_{\gamma',\alpha'}.
			\label{eq:transformation}
			\end{equation}
			!PTEX-END */
			dmS.diag(eigs,'V',concurrency_);
			dmS.check2(direction);
			
			updateKeptStates(eigs);
			
			if (verbose_ && concurrency_.root())
				std::cerr<<"Done with density-matrix diag.\n";
			
			//! transform basis: dmS^\dagger * operator matrix * dms
			rSprime = pBasis;
			if (verbose_ && concurrency_.root())
				std::cerr<<"About to changeBasis...\n";

			error_ = rSprime.changeBasis(ftransform_,dmS(),eigs,keptStates_,
			                             parameters_,concurrency_); //startEnd);
			reflectionOperator_.changeBasis(ftransform_,direction);

			if (direction == EXPAND_SYSTEM) {
				LeftRightSuperType lrs(rSprime,
										(BasisWithOperatorsType&) pBasisSummed,
										(BasisType&)lrs_.super());
				waveFunctionTransformation_.push(ftransform_,direction,lrs);
			} else {
				LeftRightSuperType lrs((BasisWithOperatorsType&) pBasisSummed,
										rSprime,
										(BasisType&)lrs_.super());
				waveFunctionTransformation_.push(ftransform_,direction,lrs);
			}
			std::ostringstream msg;
			msg<<"new size of basis="<<rSprime.size();
			progress_.printline(msg,std::cout);
		}

		void updateKeptStates(const std::vector<RealType>& eigs2)
		{
// 			std::cerr<<eigs;
// 			std::cerr<<"-----------------\n";
			std::vector<RealType> eigs = eigs2;
			std::vector<size_t> perm(eigs.size());
			Sort<std::vector<RealType> > sort;
			sort.sort(eigs,perm);
			
			size_t newKeptStates = computeKeptStates(eigs);
			size_t statesToRemove = 0;
			if (eigs.size()>=newKeptStates)
				statesToRemove = eigs.size()-newKeptStates;
			RealType discWeight = sumUpTo(eigs,statesToRemove);
// 			std::cerr<<"newKeptstates="<<newKeptStates<<"\n";
			std::ostringstream msg;
			if (newKeptStates != keptStates_) {
				// we report that the "m" value has been changed and...
				msg<<"Reducing kept states to "<<newKeptStates<<" from "<<keptStates_;
				// ... we change it:
				keptStates_ = newKeptStates;
			} else {
				// we report that the "m" value remains the same
				msg<<"Not changing kept states="<<keptStates_;
			}
			progress_.printline(msg,std::cout);
			// we report the discarded weight
			msg<<"Discarded weight (Truncation error): "<< discWeight ; 
			progress_.printline(msg,std::cout);
			
		}

		/** !PTEX-START RemovalOfStates
		Let $m_S$ (here given by \verb!keptStates_! be a fixed number that 
		corresponds to the number of states in $\mathcal{V}(S')$ that we want to keep. 
		Consider the first $m_S$ eigenvectors $w^S$, 
		 and let us call the Hilbert space spanned by them, $\mathcal{V}_R(S')$, 
		 the DMRG-reduced Hilbert space on 
		block $S'$. If $m_S\ge\#\mathcal{V}(S')$ then we keep all eigenvectors 
		and there is effectively no truncation.
		We truncate the matrices $(H^{S' {\rm new\,\,basis}})$ (and other operators as necessary)
		such that they now act on this truncated Hilbert space, $\mathcal{V}_R(S')$.
		We proceed in the same manner for the environment.
		!PTEX-END */
		//! eigenvalues are ordered in increasing order
		size_t computeKeptStates(const std::vector<RealType>& eigs) const
		{
			if (parameters_.tolerance<0) return keptStates_;
			int start = eigs.size() - keptStates_;
			if (start<0) start = 0;
			int maxToRemove = eigs.size()-parameters_.keptStatesInfinite;
			if (maxToRemove<0) maxToRemove = 0;
			size_t total = parameters_.keptStatesInfinite;
			RealType discWeight=sumUpTo(eigs,start);
			// maybe we should use int instead of size_t here!!!
			
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
			if (total>=keptStates_ || total<parameters_.keptStatesInfinite)
				return keptStates_ ; 

			return total;
		}

		RealType sumUpTo(const std::vector<RealType>& eigs,size_t x) const
		{
			RealType discWeight = 0;
			for (size_t i=0;i<x;i++)
				discWeight += fabs(eigs[i]); 
			return discWeight;
		}

		ReflectionSymmetryType& reflectionOperator_;
		const LeftRightSuperType& lrs_;
		WaveFunctionTransfType& waveFunctionTransformation_;
		ConcurrencyType& concurrency_;
		const ParametersType& parameters_;
		size_t maxConnections_;
		bool verbose_;
		ProgressIndicatorType progress_;
		size_t keptStates_;
		RealType error_;
		TransformType ftransform_;
	}; // class Truncation
	
} // namespace
/*@}*/
#endif // DMRG_TRUNCATION_H

