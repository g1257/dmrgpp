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

/*! \file MettsCollapse.h
 *
 *  Stochastics (random choices) needed for the METTS algorithm
 *
 */

#ifndef METTS_COLLAPSE_H
#define METTS_COLLAPSE_H
#include <iostream>
#include <vector>
#include <cassert>
#include "ProgressIndicator.h"
#include "PackIndices.h"
#include "Matrix.h"
#include "ProgramGlobals.h"

namespace Dmrg {
template<typename VectorWithOffsetType,
         typename MettsStochasticsType,
         typename TargettingParamsType>
class MettsCollapse  {

	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename MettsStochasticsType::PairType PairType;
	typedef typename MettsStochasticsType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename MettsStochasticsType::RealType RealType;
	typedef typename MettsStochasticsType::RngType RngType;
	typedef typename MettsStochasticsType::ModelType ModelType;
	typedef PsimagLite::Matrix<RealType> MatrixType;

	enum {EXPAND_ENVIRON=ProgramGlobals::EXPAND_ENVIRON,
		  EXPAND_SYSTEM=ProgramGlobals::EXPAND_SYSTEM};

	static const bool COLLAPSE_INTO_RANDOM_BASIS = false;

public:

	typedef PsimagLite::PackIndices PackIndicesType;

	MettsCollapse(const MettsStochasticsType& mettsStochastics,
	              const LeftRightSuperType& lrs,
	              const TargettingParamsType& targetParams)
	    : mettsStochastics_(mettsStochastics),
	      lrs_(lrs),
	      rng_(targetParams.rngSeed),
	      targetParams_(targetParams),
	      progress_("MettsCollapse"),
	      prevDirection_(ProgramGlobals::INFINITE),
	      collapseBasis_(0,0)
	{}

	bool operator()(VectorWithOffsetType& c,
	                const VectorWithOffsetType& eToTheBetaH,
	                typename PsimagLite::Vector<SizeType>::Type& block,
	                SizeType direction)
	{
		assert(direction!=ProgramGlobals::INFINITE);

		if (targetParams_.collapse.find("every")!=PsimagLite::String::npos ||
		    collapseBasis_.n_row() == 0) setCollapseBasis(block);

		internalAction(c,eToTheBetaH,block,direction,false);

		if (atBorder(direction,block)) {
			typename PsimagLite::Vector<SizeType>::Type block2;
			setBlockToBorder(block2,block);
			internalAction(c,eToTheBetaH,block2,direction,true);
		}

		for (SizeType i=0;i<block.size();i++)
			sitesSeen_.push_back(block[i]);

		if (direction==prevDirection_) return false;

		prevDirection_ = direction;

		assert(block.size()>0);
		bool allSitesSeen = checkSites(block[block.size()-1]);
		if (!allSitesSeen) return false;

		sitesSeen_.clear();

		return true;
	}

	void setNk(typename PsimagLite::Vector<SizeType>::Type& nk,
	           const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		for (SizeType i=0;i<block.size();i++)
			nk.push_back(mettsStochastics_.model().hilbertSize(block[i]));
	}

	SizeType volumeOf(const typename PsimagLite::Vector<SizeType>::Type& v) const
	{
		assert(v.size()>0);
		SizeType ret = v[0];
		for (SizeType i=1;i<v.size();i++) ret *= v[i];
		return ret;
	}

	SizeType volumeOf(const typename PsimagLite::Vector<SizeType>::Type& alphaFixed,
	                  const typename PsimagLite::Vector<SizeType>::Type& nk) const
	{
		assert(alphaFixed.size()>0);
		assert(alphaFixed.size()==nk.size());
		SizeType sum = alphaFixed[0];
		for (SizeType i=1;i<alphaFixed.size();i++)
			sum += alphaFixed[i]*nk[i-1];
		return sum;
	}

private:

	void internalAction(VectorWithOffsetType& dest2,
	                    const VectorWithOffsetType& src2,
	                    const typename PsimagLite::Vector<SizeType>::Type& block,
	                    SizeType direction,
	                    bool border) const
	{
		if (dest2.size()==0) {
			dest2 =  src2;
		}

		typename PsimagLite::Vector<SizeType>::Type nk;
		setNk(nk,block);
		SizeType volumeOfNk = volumeOf(nk);

		typename PsimagLite::Vector<RealType>::Type p(volumeOfNk,0);
		probability(p,dest2,direction,volumeOfNk,border);
		RealType sum = 0;
		for (SizeType i=0;i<p.size();i++)
			sum += p[i];
		assert(fabs(sum-1.0)<1e-6);

		SizeType volumeOfIndexFixed = mettsStochastics_.chooseRandomState(p);

		for (SizeType i=0;i<block.size();i++)
			std::cerr<<"SITES="<<block[i]<<" ";
		std::cerr<<" PROBS=";
		for (SizeType i=0;i<p.size();i++) std::cerr<<p[i]<<" ";
		std::cerr<<" CHOSEN="<<volumeOfIndexFixed<<" BORDER="<<border<<"\n";
		// m1 == indexFixed in FIXME write paper reference here

		VectorWithOffsetType dest;
		collapseVector(dest,dest2,direction,volumeOfIndexFixed,volumeOfNk,border);
		RealType x = std::norm(dest);

		if (x<1e-6 && dest.sectors()==0) {
			std::cout<<"norm of dest2= "<<std::norm(dest2)<<"\n";
			dest2 = dest;
			return;
		}

		assert(x>1e-6);
		dest2 = (1.0/x) * dest;
		assert(dest2.size()==src2.size());
	}

	void collapseVector(VectorWithOffsetType& dest2, // <<---- CPS
	                    const VectorWithOffsetType& src, // <--- MPS
	                    SizeType direction,
	                    SizeType indexFixed, // <--- m1
	                    SizeType nk,  // <-- size of the Hilbert sp. of one site
	                    bool border) const
	{
		VectorWithOffsetType dest = dest2;
		dest.populateSectors(lrs_.super());
		for (SizeType ii=0;ii<dest.sectors();ii++) {
			SizeType i0 = dest.sector(ii);
			collapseVector(dest,src,direction,i0,indexFixed,nk,border);
		}
		dest.collapseSectors();

		dest2 =  dest;
		std::cerr<<" Norm of the collapsed="<<std::norm(dest2)<<"\n";
	}

	void collapseVector(VectorWithOffsetType& w, // <<---- CPS
	                    const VectorWithOffsetType& v, // <--- MPS
	                    SizeType direction,
	                    SizeType m, // <-- non-zero sector
	                    SizeType indexFixed, // <--- m1
	                    SizeType nk,// <-- size of the Hilbert sp. of one site
	                    bool border) const
	{
		if (direction==EXPAND_SYSTEM)  collapseVectorLeft(w,v,m,indexFixed,nk,border);
		else  collapseVectorRight(w,v,m,indexFixed,nk,border);
	}

	void compare1(const VectorType& v1,const VectorType& v2) const
	{
		assert(v1.size()==v2.size());
		SizeType count = 0;
		for (SizeType i=0;i<v1.size();i++) {
			for (SizeType j=0;j<v2.size();j++) {
				if (fabs(v1[i]-v2[i])>1e-3) count++;
			}
		}
		RealType fraction = static_cast<RealType>(count)/v1.size();
		std::cout<<__FILE__<<" "<<__LINE__<<" count="<<fraction<<"%\n";
	}

	void collapseVectorLeft(VectorWithOffsetType& w, // <<---- CPS
	                        const VectorWithOffsetType& v, // <--- MPS
	                        SizeType m, // <-- non-zero sector
	                        SizeType indexFixed, // <--- m1
	                        SizeType nk, // <-- size of the Hilbert sp. of one site
	                        bool border) const
	{
		if (border)
			collapseLeftBorder(w,v,m,indexFixed,nk);
		else
			collapseVectorLeft(w,v,m,indexFixed,nk);
	}

	void collapseVectorLeft(VectorWithOffsetType& w, // <<---- CPS
	                        const VectorWithOffsetType& v, // <--- MPS
	                        SizeType m, // <-- non-zero sector
	                        SizeType indexFixed, // <--- m1
	                        SizeType nk) const // <-- size of the Hilbert sp. of one site
	{
		SizeType offset = lrs_.super().partition(m);
		int total = lrs_.super().partition(m+1) - offset;

		SizeType ns = lrs_.left().size();
		PackIndicesType packSuper(ns);
		PackIndicesType packLeft(ns/nk);
		for (SizeType i=0;i<SizeType(total);i++) {
			SizeType alpha,beta;
			packSuper.unpack(alpha,beta,lrs_.super().permutation(i+offset));

			SizeType alpha0,alpha1;
			packLeft.unpack(alpha0,alpha1,lrs_.left().permutation(alpha));

			for (SizeType alpha1Prime=0;alpha1Prime<nk;alpha1Prime++) {
				SizeType alphaPrime = packLeft.pack(alpha0,
				                                    alpha1Prime,
				                                    lrs_.left().permutationInverse());
				SizeType iprime = packSuper.pack(alphaPrime,
				                                 beta,
				                                 lrs_.super().permutationInverse());
				w[i+offset] += v[iprime]*collapseBasis_(alpha1Prime,indexFixed) *
				        collapseBasis_(alpha1,indexFixed);

			}
		}
	}

	void collapseLeftBorder(VectorWithOffsetType& w, // <<---- CPS
	                        const VectorWithOffsetType& v, // <--- MPS
	                        SizeType m, // <-- non-zero sector
	                        SizeType indexFixed, // <--- m1
	                        SizeType nk) const // <-- size of the Hilbert sp. of one site
	{
		assert(lrs_.right().size()==nk);
		SizeType offset = lrs_.super().partition(m);
		int total = lrs_.super().partition(m+1) - offset;

		SizeType ns = lrs_.left().size();
		PackIndicesType packSuper(ns);

		for (SizeType i=0;i<SizeType(total);i++) {
			SizeType alpha,beta;
			packSuper.unpack(alpha,beta,lrs_.super().permutation(i+offset));

			for (SizeType betaPrime=0;betaPrime<nk;betaPrime++) {
				SizeType iprime = packSuper.pack(alpha,
				                                 betaPrime,
				                                 lrs_.super().permutationInverse());
				w[i+offset] += v[iprime]*collapseBasis_(betaPrime,indexFixed) *
				        collapseBasis_(beta,indexFixed);
			}
		}
	}

	void collapseVectorRight(VectorWithOffsetType& w, // <<---- CPS
	                         const VectorWithOffsetType& v, // <--- MPS
	                         SizeType m, // <-- non-zero sector
	                         SizeType indexFixed, // <--- m1
	                         SizeType nk, // <-- size of the Hilbert sp. of one site
	                         bool border) const
	{
		if (border)
			collapseRightBorder(w,v,m,indexFixed,nk);
		else
			collapseVectorRight(w,v,m,indexFixed,nk);
	}

	void collapseVectorRight(VectorWithOffsetType& w, // <<---- CPS
	                         const VectorWithOffsetType& v, // <--- MPS
	                         SizeType m, // <-- non-zero sector
	                         SizeType indexFixed, // <--- m1
	                         SizeType nk) const // <-- size of the Hilbert sp. of one site
	{
		SizeType offset = lrs_.super().partition(m);
		int total = lrs_.super().partition(m+1) - offset;

		SizeType ns = lrs_.left().size();
		PackIndicesType packSuper(ns);
		PackIndicesType packRight(nk);
		for (SizeType i=0;i<SizeType(total);i++) {
			SizeType alpha,beta;
			packSuper.unpack(alpha,beta,lrs_.super().permutation(i+offset));

			SizeType beta0,beta1;
			packRight.unpack(beta0,beta1,lrs_.right().permutation(beta));

			for (SizeType beta0Prime=0;beta0Prime<nk;beta0Prime++) {
				SizeType betaPrime =  packRight.pack(beta0Prime,
				                                     beta1,
				                                     lrs_.right().permutationInverse());
				SizeType iprime = packSuper.pack(alpha,
				                                 betaPrime,
				                                 lrs_.super().permutationInverse());
				w[i+offset] += v[iprime]*collapseBasis_(beta0Prime,indexFixed) *
				        collapseBasis_(beta0,indexFixed);
			}
		}
	}

	void collapseRightBorder(VectorWithOffsetType& w, // <<---- CPS
	                         const VectorWithOffsetType& v, // <--- MPS
	                         SizeType m, // <-- non-zero sector
	                         SizeType indexFixed, // <--- m1
	                         SizeType nk) const // <-- size of the Hilbert sp. of one site
	{
		assert(lrs_.left().size()==nk);
		SizeType offset = lrs_.super().partition(m);
		int total = lrs_.super().partition(m+1) - offset;

		SizeType ns = lrs_.left().size();
		PackIndicesType packSuper(ns);

		for (SizeType i=0;i<SizeType(total);i++) {
			SizeType alpha,beta;
			packSuper.unpack(alpha,beta,lrs_.super().permutation(i+offset));

			for (SizeType alphaPrime=0;alphaPrime<nk;alphaPrime++) {
				SizeType iprime = packSuper.pack(alphaPrime,
				                                 beta,
				                                 lrs_.super().permutationInverse());
				w[i+offset] += v[iprime]*collapseBasis_(alphaPrime,indexFixed) *
				        collapseBasis_(alpha,indexFixed);
			}
		}
	}

	// p[m] = norm2 of the collapsed_m
	void probability(typename PsimagLite::Vector<RealType>::Type& p,
	                 const VectorWithOffsetType& src,
	                 SizeType direction,
	                 SizeType volumeOfNk,
	                 bool border) const
	{
		RealType tmp = std::norm(src);
		if (fabs(tmp-1.0)>1e-3)
			std::cerr<<"probability "<<tmp<<"\n";

		RealType sum = 0;
		for (SizeType alpha=0;alpha<volumeOfNk;alpha++) {
			VectorWithOffsetType dest;
			collapseVector(dest,src,direction,alpha,volumeOfNk,border);
			RealType x = std::norm(dest);
			sum += x*x;
			p[alpha] = x*x;
		}
		if (fabs(sum-1.0)>1e-3)
			std::cerr<<"probability sum="<<sum<<"\n";
		assert(fabs(sum)>1e-6);
		for (SizeType alpha=0;alpha<p.size();alpha++) {
			p[alpha] /= sum;
			assert(p[alpha]>=0 && p[alpha]<=1);
		}
	}

	bool checkSites(SizeType site) const
	{
		std::cerr<<"MettsCollapse: SITES SEEN ";
		for (SizeType i=0;i<sitesSeen_.size();++i)
			std::cerr<<sitesSeen_[i]<<" ";
		std::cerr<<"\n";

		SizeType sitesPerBlock = mettsStochastics_.model().params().sitesPerBlock;
		for (SizeType i=sitesPerBlock;i<site+1;i++) {
			bool seen = (std::find(sitesSeen_.begin(),sitesSeen_.end(),i) != sitesSeen_.end());
			if (!seen) return false;
		}

		return true;
	}

	void setCollapseBasis(const typename PsimagLite::Vector<SizeType>::Type& block)
	{
		SizeType nk = 1;
		for (SizeType i=0;i<block.size();i++)
			nk *= mettsStochastics_.model().hilbertSize(block[i]);

		collapseBasis_.resize(nk,nk);
		for (SizeType i=0;i<nk;i++)
			for (SizeType j=0;j<nk;j++)
				collapseBasis_(i,j) = (i==j) ? 1.0 : 0.0;

		assert(block.size()>0);
		SizeType site = block[0];

		if (targetParams_.collapse.find("random")!=PsimagLite::String::npos)
			rotationNd(collapseBasis_,
			           mettsStochastics_.model().hilbertSize(site),
			           block.size());

		if (targetParams_.collapse.find("particle")!=PsimagLite::String::npos)
			particleCollapse(collapseBasis_);

		std::cout<<"Collapse basis:\n";
		std::cout<<collapseBasis_;
		checkBasis();
	}

	void particleCollapse(MatrixType& m) const
	{
		if (m.n_row()!=4 || m.n_col()!=4)
			throw PsimagLite::RuntimeError("particleCollapse: only for 4 states\n");

		MatrixType m2(m.n_row(),m.n_col());
		RealType theta = M_PI*rng_();
		rotation2d(m2,1,2,theta);
		m=m2;
	}

	void rotationNd(MatrixType& m,
	                SizeType oneSiteHilbertSize,
	                SizeType blockSize) const
	{
		if (blockSize == 1)
			return rotationNd(m,oneSiteHilbertSize);

		MatrixType mold;
		for (SizeType i=0;i<blockSize;i++) {
			MatrixType aux1(oneSiteHilbertSize,oneSiteHilbertSize);
			rotationNd(aux1,oneSiteHilbertSize);
			if (i == 0) {
				m = aux1;
			} else {
				m.reset(0,0);
				outerProduct(m,mold,aux1);
			}
			mold = m;
		}
	}

	void rotationNd(MatrixType& m,SizeType hilbertSize) const
	{
		for (SizeType i=0;i<hilbertSize;i++) {
			MatrixType aux1(m.n_row(),m.n_col());
			RealType theta = M_PI*rng_();
			SizeType i1 = i;
			SizeType i2 = i+1;
			if (i==hilbertSize-1) i2=0;
			rotation2d(aux1,i1,i2,theta);
			if (i==0) m = aux1;
			else m = (m*aux1);
		}
	}

	void rotation2d(MatrixType& m,SizeType x,SizeType y,const RealType& theta) const
	{
		std::cout<<"Theta="<<theta<<"\n";
		for (SizeType i=0;i<m.n_row();i++) m(i,i) = 1.0;
		m(x,x) = m(y,y) = cos(theta);
		m(x,y) = sin(theta);
		m(y,x) = -sin(theta);
	}

	bool atBorder(SizeType direction,
	              const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		typename PsimagLite::Vector<SizeType>::Type nk;
		setNk(nk,block);
		SizeType volumeOfNk = volumeOf(nk);
		bool b1 = (direction==EXPAND_SYSTEM && lrs_.right().size()==volumeOfNk);
		bool b2 = (direction==EXPAND_ENVIRON && lrs_.left().size()==volumeOfNk);
		return (b1 || b2);
	}

	void setBlockToBorder(typename PsimagLite::Vector<SizeType>::Type& block2,
	                      const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		block2 = block;
		assert(block.size()>0);
		SizeType site = block[block.size()-1];
		bool leftCorner = (site+2==lrs_.super().block().size()) ? false : true;
		int offset = (leftCorner) ? -block.size() : block.size();
		for (SizeType i=0;i<block2.size();i++) {
			block2[i] = block[i] + offset;
		}
	}

	void checkBasis() const
	{
		MatrixType transpose;
		transposeConjugate(transpose,collapseBasis_);
		MatrixType shouldBeI = collapseBasis_ * transpose;
		if (isTheIdentity(shouldBeI)) return;
		throw PsimagLite::RuntimeError("checkBasis\n");
	}

	const MettsStochasticsType& mettsStochastics_;
	const LeftRightSuperType& lrs_;
	mutable RngType rng_;
	const TargettingParamsType& targetParams_;
	PsimagLite::ProgressIndicator progress_;
	SizeType prevDirection_;
	MatrixType collapseBasis_;
	typename PsimagLite::Vector<SizeType>::Type sitesSeen_;
};  //class MettsCollapse
} // namespace Dmrg
/*@}*/
#endif //METTS_COLLAPSE_H

