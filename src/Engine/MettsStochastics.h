/*
Copyright (c) 2009-2012, 2013, UT-Battelle, LLC
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
template<typename ModelType_,typename RngType_>
class MettsStochastics  {

public:

	typedef std::pair<SizeType,SizeType> PairType;
	typedef ModelType_ ModelType;
	typedef typename ModelType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelType::HilbertBasisType HilbertBasisType;
	typedef RngType_ RngType;
	typedef typename RngType::LongType LongType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

	MettsStochastics(const ModelType& model,
	                 int long seed,
	                 const VectorSizeType& pure)
	    : model_(model),
	      rng_(seed),
	      pure_(pure),
	      progress_("MettsStochastics"),
	      addedSites_(0)
	{}

	const ModelType& model() const { return model_; }

	SizeType chooseRandomState(SizeType site) const
	{
		if (site < pure_.size()) return pure_[site];

		return SizeType(rng_()*model_.hilbertSize(site));
	}

	SizeType chooseRandomState(const VectorRealType& probs) const
	{
		RealType r = rng_();
		RealType s1 = 0;
		RealType s2 = 0;
		for (SizeType i=0;i<probs.size();++i) {
			s2 = s1 + probs[i];
			if (s1<r && r<=s2) return i;
			s1 = s2;
		}

		PsimagLite::String s(__FILE__);
		s += PsimagLite::String(" ") + ttos(__LINE__) + " " + __FUNCTION__ +
		        " Probabilities don't amount to 1\n";
		throw PsimagLite::RuntimeError(s.c_str());
	}

	// call only from INFINITE
	void update(const QnType& qn,
	            const typename PsimagLite::Vector<SizeType>::Type& block1,
	            const typename PsimagLite::Vector<SizeType>::Type& block2,
	            SizeType seed)
	{
		if (addedSites_.size()==0) {
			pureStates_.resize(block2[block2.size()-1]+block2.size()+1);
			initialSetOfPures(seed);
			for (SizeType i=0;i<block1.size();i++)
				for (SizeType j=0;j<block1.size();j++)
					addedSites_.push_back(block1[i]+j-block1.size());
			for (SizeType i=0;i<block2.size();i++)
				for (SizeType j=0;j<block2.size();j++)
					addedSites_.push_back(block2[i]+j+block2.size());
		}

		for (SizeType i=0;i<block1.size();i++)
			addedSites_.push_back(block1[i]);
		for (SizeType i=0;i<block2.size();i++)
			addedSites_.push_back(block2[i]);

		qnVsSize_.resize(addedSites_.size() + 1, QnType::zero());
		qnVsSize_[addedSites_.size()]=qn;
	}

	void setCollapseBasis(typename PsimagLite::Vector<RealType>::Type& collapseBasisWeights,
	                      SizeType site) const
	{
		SizeType nk = model_.hilbertSize(site);
		for (SizeType alpha=0;alpha<nk;alpha++) {
			RealType randomNumber = rng_();
			collapseBasisWeights[alpha] = randomNumber;
		}
		RealType norm1  = 1.0/PsimagLite::norm(collapseBasisWeights);
		collapseBasisWeights *= norm1;
	}

private:

	void initialSetOfPures(LongType)
	{
		for (SizeType i=0;i<pureStates_.size();i++)
			pureStates_[i] = SizeType(rng_()*model_.hilbertSize(i));
	}

	const ModelType& model_;
	mutable RngType rng_;
	const VectorSizeType& pure_;
	PsimagLite::ProgressIndicator progress_;
	typename PsimagLite::Vector<SizeType>::Type pureStates_;
	typename PsimagLite::Vector<SizeType>::Type addedSites_;
	VectorQnType qnVsSize_;
};  //class MettsStochastics
} // namespace Dmrg
/*@}*/
#endif //METTS_STOCHASTICS_H

