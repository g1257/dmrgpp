// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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

/*! \file MettsStochastics.h
 *
 *  Stochastics (random choices) needed for the METTS algorithm 
 *
 */
 
#ifndef METTS_STOCHASTICS_H
#define METTS_STOCHASTICS_H
#include <iostream>
#include <vector>
#include "Random48.h"

namespace Dmrg {
	template<typename ModelType>
	class MettsStochastics  {
	public:
		typedef std::pair<size_t,size_t> PairType;
		typedef typename ModelType::RealType RealType;

		MettsStochastics(const ModelType& model)
		: model_(model), random48_(34328811),
		  pureStates_(model_.geometry().numberOfSites()),
		  addedSites_(0)
		{
			size_t addedBlockSize = 1;
			model_.setNaturalBasis(basisOfOneSite_,quantumNumbsOneSite_,addedBlockSize);
			
			for (size_t i=0;i<pureStates_.size();i++) 
				pureStates_[i] = size_t(random48_()*basisOfOneSite_.size());
		}

		size_t chooseRandomState(size_t i) const
		{
			return basisOfOneSite_[pureStates_[i]];
		}
		
		void update(size_t qn,const PairType& sites)
		{
			addedSites_.push_back(sites.first);
			addedSites_.push_back(sites.second);
			
			// fix target quantum number
			size_t symm = getSymmetry();
			
			while(symm!=qn) {
				raiseOrLowerSymm(sites.first,(symm<qn));
				symm = getSymmetry();
				if (symm==qn) break;
				raiseOrLowerSymm(sites.second,(symm<qn));
				symm = getSymmetry();
			}
		}

	private:
		// assumes states in basisOfOneSite_ are ordered in increasing
		// symmetry
		void raiseOrLowerSymm(size_t site,bool raiseSymm)
		{
			if (raiseSymm && pureStates_[site]<basisOfOneSite_.size()-1) {
				pureStates_[site]++;
				return;
			}
			if (pureStates_[site]>0) pureStates_[site]--;
		}
		
		// assumes local symmetry througout
		size_t getSymmetry() const
		{
			size_t sum = 0;
			for (size_t i=0;i<addedSites_.size();i++) {
				sum += quantumNumbsOneSite_[pureStates_[addedSites_[i]]];
			}
			return sum;
		}
		
		const ModelType& model_;
		PsimagLite::Random48<RealType> random48_;
		std::vector<size_t> pureStates_;
		std::vector<size_t> addedSites_;
		std::vector<size_t> quantumNumbsOneSite_;
		typename ModelType::HilbertBasisType basisOfOneSite_;
	};  //class MettsStochastics
} // namespace Dmrg
/*@}*/
#endif //METTS_STOCHASTICS_H
