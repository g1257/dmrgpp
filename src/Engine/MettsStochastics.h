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
#include <algorithm>
#include "TypeToString.h"
#include "Utils.h"

namespace Dmrg {
	template<typename ModelType,typename RngType_>
	class MettsStochastics  {

	public:

		typedef std::pair<size_t,size_t> PairType;
		typedef typename ModelType::RealType RealType;
		typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
		typedef typename ModelType::HilbertBasisType HilbertBasisType;
		typedef RngType_ RngType;
		typedef typename RngType::LongType LongType;

		MettsStochastics(const ModelType& model,int long long seed)
		: model_(model),
		  rng_(seed),
		  progress_("MettsStochastics",0),
		  addedSites_(0)
		{}

		size_t hilbertSize(size_t site) const { return model_.hilbertSize(site); }

		size_t chooseRandomState(size_t site) const
		{
			typename PsimagLite::Vector<size_t>::Type quantumNumbsOneSite;
			HilbertBasisType basisOfOneSite;
			basisForOneSite(quantumNumbsOneSite,basisOfOneSite,site);
			size_t tmp = size_t(rng_()*model_.hilbertSize(site));
			assert(tmp<basisOfOneSite.size());
			return basisOfOneSite[tmp];
		}

		// recommendation for understanding this
		// think of this first as having basisOfOneSite[i] = i
		size_t chooseRandomState(const typename PsimagLite::Vector<RealType>::Type& probs,size_t site) const
		{
			typename PsimagLite::Vector<size_t>::Type quantumNumbsOneSite;
			HilbertBasisType basisOfOneSite;
			basisForOneSite(quantumNumbsOneSite,basisOfOneSite,site);

			std::cerr<<"basisOfOneSite ";
			for (size_t i=0;i<basisOfOneSite.size();i++)
				std::cerr<<basisOfOneSite[i]<<" ";
			std::cerr<<"\n";

			RealType r = rng_();
			std::cout<<"RANDOM="<<r<<"\n";
			RealType s1 = 0;
			RealType s2 = 0;
			for (size_t i=0;i<probs.size();++i) {
				s2 = s1 + probs[i];
				if (s1<r && r<=s2) return basisOfOneSite[i];
				s1 = s2;
			}
			PsimagLite::String s(__FILE__);
			s += PsimagLite::String(" ") + ttos(__LINE__) + " " + __FUNCTION__ +
			     " Probabilities don't amount to 1\n";
			throw PsimagLite::RuntimeError(s.c_str());
		}

		// call only from INFINITE
		void update(size_t qn,const typename PsimagLite::Vector<size_t>::Type& block1,const typename PsimagLite::Vector<size_t>::Type& block2,size_t seed)
		{
			if (addedSites_.size()==0) {
				pureStates_.resize(block2[block2.size()-1]+2*block2.size());
				initialSetOfPures(seed);
				for (size_t i=0;i<block1.size();i++)
					for (size_t j=0;j<block1.size();j++)
						addedSites_.push_back(block1[i]+j-block1.size());
				for (size_t i=0;i<block2.size();i++)
					for (size_t j=0;j<block2.size();j++)
						addedSites_.push_back(block2[i]+j+block2.size());
			}

			for (size_t i=0;i<block1.size();i++)
				addedSites_.push_back(block1[i]);
			for (size_t i=0;i<block2.size();i++)
				addedSites_.push_back(block2[i]);
			qnVsSize_.resize(addedSites_.size()+1,0);
			qnVsSize_[addedSites_.size()]=qn;
		}

		void setCollapseBasis(typename PsimagLite::Vector<RealType>::Type& collapseBasisWeights,size_t site) const
		{
			size_t nk = model_.hilbertSize(site);
			for (size_t alpha=0;alpha<nk;alpha++) {
				RealType randomNumber = rng_();
				collapseBasisWeights[alpha] = randomNumber;
			}
			RealType norm1  = 1.0/PsimagLite::norm(collapseBasisWeights);
			collapseBasisWeights *= norm1;
		}

	private:

		void basisForOneSite(typename PsimagLite::Vector<size_t>::Type& quantumNumbsOneSite,
		                     HilbertBasisType& basisOfOneSite,size_t site) const
		{
			typename PsimagLite::Vector<size_t>::Type block(1,site);
			model_.setNaturalBasis(basisOfOneSite,quantumNumbsOneSite,block);
		}

		void initialSetOfPures(LongType seed)
		{
			for (size_t i=0;i<pureStates_.size();i++)
				pureStates_[i] = size_t(rng_()*model_.hilbertSize(i));
		}

		const ModelType& model_;
		mutable RngType rng_;
		PsimagLite::ProgressIndicator progress_;
		typename PsimagLite::Vector<size_t>::Type pureStates_;
		typename PsimagLite::Vector<size_t>::Type addedSites_;
		typename PsimagLite::Vector<size_t>::Type qnVsSize_;
	};  //class MettsStochastics
} // namespace Dmrg
/*@}*/
#endif //METTS_STOCHASTICS_H
