
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
		
		enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
		EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM};
		
	public:

		typedef typename DensityMatrixType::BuildingBlockType TransformType;

		Truncation(LeftRightSuperType& lrs,
		           WaveFunctionTransfType& waveFunctionTransformation,
		           ConcurrencyType& concurrency,
		           const ParametersType& parameters,
		           bool verbose)
		: lrs_(lrs),
		  waveFunctionTransformation_(waveFunctionTransformation),
		  concurrency_(concurrency),
		  parameters_(parameters),
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
// 			TransformType& ftransform,
		                     size_t keptStates,
		                     size_t direction)
		{
			keptStates_ = keptStates;
			if (direction==EXPAND_SYSTEM) {
				changeAndTruncateBasis(rSprime,
				                       target,
				                       lrs_.left(),
				                       lrs_.right(),
				                       direction);
			} else {
				changeAndTruncateBasis(rSprime,
				                       target,
				                       lrs_.right(),
				                       lrs_.left(),
				                       direction);
			}
		}
		
		const TransformType& transform() const { return ftransform_; }
		
		const RealType& error() const { return error_; }
		
	private:

		//! Truncate basis 
		void changeAndTruncateBasis(BasisWithOperatorsType &rSprime,
		                            const TargettingType& target,
		                            BasisWithOperatorsType const &pBasis,
		                            BasisWithOperatorsType const &pBasisSummed,
		                            size_t direction)
		{
			DensityMatrixType dmS(target,pBasis,pBasisSummed,lrs_.super(),direction);
			dmS.check(direction);
			
			if (verbose_ && concurrency_.root()) {
				std::cerr<<"Trying to diagonalize density-matrix with size=";
				std::cerr<<dmS.rank()<<"\n";
			}
			std::vector<RealType> eigs;
			dmS.diag(eigs,'V',concurrency_);
			dmS.check2(direction);
			
			updateKeptStates(eigs);
			
			if (verbose_ && concurrency_.root())
				std::cerr<<"Done with density-matrix diag.\n";
			
			//! transform basis: dmS^\dagger * operator matrix * dms
			rSprime = pBasis;
			if (verbose_ && concurrency_.root())
				std::cerr<<"About to changeBasis...\n";

			error_ = rSprime.changeBasis(ftransform_,
			                                        dmS(),
			                                        eigs,
			                                        keptStates_,
			                                        parameters_,
			                                        concurrency_);
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
		
		void updateKeptStates(const std::vector<RealType>& eigs)
		{
			if (parameters_.tolerance<0) return;
			size_t min = parameters_.keptStatesInfinite;
			size_t total = min;
			for (size_t i=min;i<eigs.size();i++) {
				if (fabs(eigs[i])<parameters_.tolerance) {
					total = i;
					break;
				}
			}
			if (total>=keptStates_ || total<min)  return;
			std::ostringstream msg;
			msg<<"Updating kept states to "<<total<<" from "<<keptStates_;
			progress_.printline(msg,std::cout);
			
			keptStates_ = total;
		}
		
		const LeftRightSuperType& lrs_;
		WaveFunctionTransfType& waveFunctionTransformation_;
		ConcurrencyType& concurrency_;
		const ParametersType& parameters_;
		bool verbose_;
		ProgressIndicatorType progress_;
		size_t keptStates_;
		RealType error_;
		TransformType ftransform_;
	}; // class Truncation
	
} // namespace
/*@}*/
#endif // DMRG_TRUNCATION_H

