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

/*! \file MettsStochastics.h
 *
 *  Stochastics (random choices) needed for the METTS algorithm 
 *
 */
 
#ifndef METTS_STOCHASTICS_H
#define METTS_STOCHASTICS_H
#include <iostream>
#include <vector>
#include <stdexcept>
#include "ProgressIndicator.h"
#include "Random48.h"
#include <algorithm>
#include "TypeToString.h"

namespace Dmrg {
	template<typename ModelType>
	class MettsStochastics  {
	public:
		typedef std::pair<size_t,size_t> PairType;
		typedef typename ModelType::RealType RealType;
		typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
		typedef typename PsimagLite::Random48<RealType> RngType;
		typedef typename RngType::LongType LongType;
		typedef typename ModelType::HilbertBasisType HilbertBasisType;

		MettsStochastics(const ModelType& model,const LongType& seed)
		: model_(model),
		  progress_("MettsStochastics",0),
		  random48_(seed),
		  addedSites_(0)
		{}

		size_t hilbertSize(size_t site) const { return model_.hilbertSize(site); }

		size_t chooseRandomState(size_t site) const
		{
			std::vector<size_t> quantumNumbsOneSite;
			HilbertBasisType basisOfOneSite;
			basisForOneSite(quantumNumbsOneSite,basisOfOneSite,site);
			return basisOfOneSite[pureStates_[site]];
		}

		// recommendation for understanding this
		// think of this first as having basisOfOneSite[i] = i
		size_t chooseRandomState(const std::vector<RealType>& probs,size_t site) const
		{
			std::vector<size_t> quantumNumbsOneSite;
			HilbertBasisType basisOfOneSite;
			basisForOneSite(quantumNumbsOneSite,basisOfOneSite,site);

			RealType r = random48_();
			RealType s1 = 0;
			RealType s2 = 0;
			for (size_t i=0;i<probs.size();++i) {
				s2 = s1 + probs[i];
				if (s1<r && r<=s2) return basisOfOneSite[i];
				s1 = s2;
			}
			std::string s(__FILE__);
			s += std::string(" ") + ttos(__LINE__) + " " + __FUNCTION__ +
			     " Probabilities don't amount to 1\n";
			throw std::runtime_error(s.c_str());
		}

		void update(size_t qn,const PairType& sites)
		{
			std::vector<size_t> currentSites;
			currentSites.push_back(sites.first);
			currentSites.push_back(sites.second);

			if (addedSites_.size()==0) {
				if (sites.first!=1)
					throw std::runtime_error("MettsStochastics::update(...): must start from 0\n");
				pureStates_.resize(sites.second+2);
				for (size_t i=0;i<pureStates_.size();i++)
					pureStates_[i] = size_t(random48_()*model_.hilbertSize(i));

				addedSites_.push_back(sites.first-1);
				addedSites_.push_back(sites.second+1);
				currentSites.push_back(sites.first-1);
				currentSites.push_back(sites.second+1);
			}

			if (sites.first!=sites.second) {
				addedSites_.push_back(sites.first);
				addedSites_.push_back(sites.second);
				qnVsSize_.resize(addedSites_.size()+1,0);
				qnVsSize_[addedSites_.size()]=qn;
				getStochasticsUpToThisPoint(qn,currentSites);
				return; // INFINITE
			}

			//FINITE: all this is bogus because we're not using
			// pureStates_ anymore in the finite phase
			if (std::find(addedSites_.begin(),addedSites_.end(),
				sites.first) != addedSites_.end()) {
				//getStochasticsForLattice();
				addedSites_.clear();
			}
			addedSites_.push_back(sites.first);
		}

		void setCollapseBasis(std::vector<RealType>& collapseBasisWeights,size_t site) const
		{
			size_t nk = model_.hilbertSize(site);
			for (size_t alpha=0;alpha<nk;alpha++) {
				RealType randomNumber = random48_();
				collapseBasisWeights[alpha] = randomNumber;
			}
			RealType norm1  = 1.0/PsimagLite::norm(collapseBasisWeights);
			collapseBasisWeights *= norm1;
		}


	private:

		void basisForOneSite(std::vector<size_t>& quantumNumbsOneSite,HilbertBasisType& basisOfOneSite,size_t site) const
		{
			std::vector<size_t> block(1,site);
			model_.setNaturalBasis(basisOfOneSite,quantumNumbsOneSite,block);
		}
// 		void getStochasticsForLattice()
// 		{
// 			for (size_t i=0;i<pureStates_.size();i++) 
// 				pureStates_[i] = size_t(random48_()*basisOfOneSite_.size());
// 			size_t sys = addedSites_.size()/2;
// 			size_t env = sys;
// 			addedSites_.clear();
// 			addedSites_.push_back(0);
// 			addedSites_.push_back(2*env-1);
// 			for (size_t i=1;i<sys;i++) {
// 				addedSites_.push_back(i);
// 				addedSites_.push_back(2*env-i-1);
// 				std::vector<size_t> sites;
// 				sites.push_back(i);
// 				sites.push_back(2*env-i-1);
// 				getStochasticsUpToThisPoint(qnVsSize_[addedSites_.size()],sites);
// 			}
// 		}

		void getStochasticsUpToThisPoint(size_t qn,
		                                 const std::vector<size_t>& currentSites) 
		{
			// fix target quantum number
			size_t symm = getSymmetry();
			size_t counter = 0;
			while(symm!=qn) {
				counter++;
				if (counter>1e6) {
					std::string s(__FILE__);
					s += std::string(" ") + ttos(__LINE__) + std::string(" ") + 
					std::string(__FUNCTION__);
					s += std::string(" too many iterations\n");
					throw std::runtime_error(s.c_str());
				}
				for (size_t i=0;i<currentSites.size();i++) {
					size_t thisSite = currentSites[i];
					if (i==1 && currentSites[0]==currentSites[1]) break;
					raiseOrLowerSymm(thisSite,(symm<qn));
					symm = getSymmetry();
					if (symm==qn) break;
				}
// 				raiseOrLowerSymm(sites.first,(symm<qn));
// 				symm = getSymmetry();
// 				if (sites.second==sites.first) continue;
// 				if (symm==qn) break;
// 				raiseOrLowerSymm(sites.second,(symm<qn));
// 				symm = getSymmetry();
			}
			std::ostringstream msg;
			msg<<"targetQn="<<qn<<" sites="<<addedSites_.size()<<" Pure=";
			for (size_t i=0;i<pureStates_.size();i++)
				msg<<pureStates_[i]<<" ";
			progress_.printline(msg,std::cout);
		}

		// assumes states in basisOfOneSite_ are ordered in increasing
		// symmetry
		void raiseOrLowerSymm(size_t site,bool raiseSymm)
		{
			if (raiseSymm) {
				if (pureStates_[site]<model_.hilbertSize(site)-1)
					pureStates_[site]++;
				return;
			}
			
			if (pureStates_[site]>0) pureStates_[site]--;
		}

		// assumes local symmetry througout
		size_t getSymmetry() const
		{
			std::vector<size_t> quantumNumbsOneSite;
			HilbertBasisType basisOfOneSite;
			
			size_t sum = 0;
			for (size_t i=0;i<addedSites_.size();i++) {
				basisForOneSite(quantumNumbsOneSite,basisOfOneSite,addedSites_[i]);
				sum += quantumNumbsOneSite[pureStates_[addedSites_[i]]];
			}
			return sum;
		}

		const ModelType& model_;
		PsimagLite::ProgressIndicator progress_;
		RngType random48_;
		std::vector<size_t> pureStates_;
		std::vector<size_t> addedSites_;
// 		std::vector<size_t> quantumNumbsOneSite_;
// 		typename ModelType::HilbertBasisType basisOfOneSite_;
		std::vector<size_t> qnVsSize_;
	};  //class MettsStochastics
} // namespace Dmrg
/*@}*/
#endif //METTS_STOCHASTICS_H
