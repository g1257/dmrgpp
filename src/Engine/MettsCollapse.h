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
#include "Random48.h"

namespace Dmrg {
	template<typename VectorWithOffsetType,typename MettsStochasticsType>
	class MettsCollapse  {

		typedef typename VectorWithOffsetType::VectorType VectorType;
		typedef typename MettsStochasticsType::PairType PairType;
		typedef typename MettsStochasticsType::LeftRightSuperType LeftRightSuperType;
		typedef typename MettsStochasticsType::RealType RealType;
		typedef typename MettsStochasticsType::RngType RngType;
		typedef PsimagLite::Matrix<RealType> MatrixType;

		enum {EXPAND_ENVIRON=ProgramGlobals::EXPAND_ENVIRON,
		      EXPAND_SYSTEM=ProgramGlobals::EXPAND_SYSTEM};
		
		static const bool COLLAPSE_INTO_RANDOM_BASIS = false;

	public:

		typedef PsimagLite::PackIndices PackIndicesType;

		MettsCollapse(const MettsStochasticsType& mettsStochastics,
			      const LeftRightSuperType& lrs,
			      RngType& rng)
		: mettsStochastics_(mettsStochastics),
		  lrs_(lrs),
		  rng_(rng),
		  progress_("MettsCollapse",0),
		  prevDirection_(ProgramGlobals::INFINITE),
		  collapseBasis_(mettsStochastics_.hilbertSize(0),mettsStochastics_.hilbertSize(0))
		{
			setCollapseBasis();
		}

		bool operator()(VectorWithOffsetType& c,
		                const VectorWithOffsetType& eToTheBetaH,
		                size_t site,
		                size_t direction)
		{

			internalAction(c,eToTheBetaH,site,direction);
			sitesSeen_.push_back(site);
			
			if (direction==prevDirection_) return false;
			prevDirection_ = direction;

			bool allSitesSeen = checkSites(site);
			if (!allSitesSeen) return false;

			sitesSeen_.clear();

			return true;
		}

	private:

		void internalAction(VectorWithOffsetType& dest2,
		                    const VectorWithOffsetType& src2,
		                    size_t site,
		                    size_t direction) const
		{
			size_t nk = mettsStochastics_.hilbertSize(site);
			if (dest2.size()==0) {
				dest2 =  src2;
			}

			std::vector<RealType> p(nk,0);
			probability(p,src2,direction,nk);
			RealType sum = 0;
			for (size_t i=0;i<p.size();i++)
				sum += p[i];
			assert(fabs(sum-1.0)<1e-6);

			VectorWithOffsetType dest;
			size_t indexFixed = mettsStochastics_.chooseRandomState(p,site);
			std::cerr<<"SITE="<<site<<" PROBS=";
			for (size_t i=0;i<p.size();i++) std::cerr<<p[i]<<" ";
			std::cerr<<" CHOSEN="<<indexFixed<<"\n";
			 // m1 == indexFixed in FIXME write paper reference here

			collapseVector(dest,dest2,direction,indexFixed,nk,true);
			RealType x = std::norm(dest);
			assert(x>1e-6);
			dest2 = (1.0/x) * dest;
			assert(dest2.size()==src2.size());
		}

		void collapseVector(VectorWithOffsetType& dest, // <<---- CPS
		                    const VectorWithOffsetType& src, // <--- MPS
		                    size_t direction,
		                    size_t indexFixed, // <--- m1
		                    size_t nk,  // <-- size of the Hilbert sp. of one site
		                    bool option = false) const
		{
			assert(src.sectors()==1);
			dest = src;
			for (size_t ii=0;ii<src.sectors();ii++) {
				size_t i0 = src.sector(ii);
				VectorType vdest,vsrc;
				src.extract(vsrc,i0);
				collapseVector(vdest,vsrc,direction,i0,indexFixed,nk,option);
				dest.setDataInSector(vdest,i0);
			}
			//assert(std::norm(dest)>1e-6);
		}

		void collapseVector(VectorType& w, // <<---- CPS
		                    const VectorType& v, // <--- MPS
		                    size_t direction,
		                    size_t m, // <-- non-zero sector
		                    size_t indexFixed, // <--- m1
		                    size_t nk,// <-- size of the Hilbert sp. of one site
		                    bool option) const
		{
			if (direction==EXPAND_SYSTEM)  collapseVectorLeft(w,v,m,indexFixed,nk);
			else  collapseVectorRight(w,v,m,indexFixed,nk);
			//if (option && direction!=prevDirection_) compare1(w,v);
//			assert(fabs(PsimagLite::norm(w))>1e-6);
		}

		void compare1(const VectorType& v1,const VectorType& v2) const
		{
			assert(v1.size()==v2.size());
			size_t count = 0;
			for (size_t i=0;i<v1.size();i++) {
				for (size_t j=0;j<v2.size();j++) {
					if (fabs(v1[i]-v2[i])>1e-3) count++;
				}
			}
			RealType fraction = static_cast<RealType>(count)/v1.size();
			std::cout<<__FILE__<<" "<<__LINE__<<" count="<<fraction<<"%\n";
		}

		void collapseVectorLeft(VectorType& w, // <<---- CPS
					const VectorType& v, // <--- MPS
					size_t m, // <-- non-zero sector
					size_t indexFixed, // <--- m1
					size_t nk) const // <-- size of the Hilbert sp. of one site
		{
			size_t offset = lrs_.super().partition(m);
			int total = lrs_.super().partition(m+1) - offset;

			size_t ns = lrs_.left().size();
			PackIndicesType packSuper(ns);
			PackIndicesType packLeft(ns/nk);
			w.resize(total);
			for (size_t i=0;i<size_t(total);i++) {
				w[i] = 0;
				size_t alpha,beta;
				packSuper.unpack(alpha,beta,lrs_.super().permutation(i+offset));

				size_t alpha0,alpha1;
				packLeft.unpack(alpha0,alpha1,lrs_.left().permutation(alpha));

				RealType sum = 0;
				for (size_t alpha1Prime=0;alpha1Prime<nk;alpha1Prime++) {
					size_t alphaPrime = packLeft.pack(alpha0,alpha1Prime,lrs_.left().permutationInverse());
					size_t iprime = packSuper.pack(alphaPrime,beta,lrs_.super().permutationInverse());
					if (iprime<offset || iprime>=offset+total) continue;
					sum += v[iprime-offset]*collapseBasis_(alpha1Prime,indexFixed);
				}
				w[i] = sum * collapseBasis_(alpha1,indexFixed);
			}
		}

		void collapseVectorRight(VectorType& w, // <<---- CPS
		                         const VectorType& v, // <--- MPS
		                         size_t m, // <-- non-zero sector
		                         size_t indexFixed, // <--- m1
		                         size_t nk) const // <-- size of the Hilbert sp. of one site
		{
			size_t offset = lrs_.super().partition(m);
			int total = lrs_.super().partition(m+1) - offset;

			size_t ns = lrs_.left().size();
			PackIndicesType packSuper(ns);
			PackIndicesType packRight(nk);
			w.resize(total);
			for (size_t i=0;i<size_t(total);i++) {
				w[i] = 0;
				size_t alpha,beta;
				packSuper.unpack(alpha,beta,lrs_.super().permutation(i+offset));

				size_t beta0,beta1;
				packRight.unpack(beta0,beta1,lrs_.right().permutation(beta));

				RealType sum = 0;
				for (size_t beta0Prime=0;beta0Prime<nk;beta0Prime++) {
					size_t betaPrime =  packRight.pack(beta0Prime,beta1,lrs_.right().permutationInverse());
					size_t iprime = packSuper.pack(alpha,betaPrime,lrs_.super().permutationInverse());
					if (iprime<offset || iprime>=offset+total) continue;
					sum += v[iprime-offset]*collapseBasis_(beta0Prime,indexFixed);
				}

				w[i] = sum * collapseBasis_(beta0,indexFixed);
			}
		}

		// p[m] = norm2 of the collapsed_m
		void probability(std::vector<RealType>& p,
		                 const VectorWithOffsetType& src,
		                 size_t direction,
		                 size_t nk) const // <-- size of the Hilbert sp. of one site
		{
			RealType sum = 0;
			for (size_t alpha=0;alpha<nk;alpha++) {
				VectorWithOffsetType dest;
				collapseVector(dest,src,direction,alpha,nk);
				RealType x = std::norm(dest);
				sum += x*x;
				p[alpha] = x*x;
			}
			assert(fabs(sum)>1e-6);
			for(size_t alpha=0;alpha<p.size();alpha++) p[alpha] /= sum;
		}

		bool checkSites(size_t site) const
		{
			for (size_t i=1;i<site+1;i++) {
				bool seen = (std::find(sitesSeen_.begin(),sitesSeen_.end(),i) != sitesSeen_.end());
				if (!seen) return false;
			}
			return true;
		}

		void setCollapseBasis()
		{
			size_t nk = collapseBasis_.n_row();

			for (size_t i=0;i<nk;i++)
				for (size_t j=0;j<nk;j++)
					collapseBasis_(i,j) = (i==j) ? 1.0 : 0.0;
			if (nk!=4) return;

			rotation4d(collapseBasis_);

		}

		void rotation4d(MatrixType& m) const
		{
			assert(m.n_row()==4 && m.n_col()==4);

			MatrixType aux1(m.n_row(),m.n_col());
			RealType theta = M_PI*rng_();
			rotation2d(aux1,0,1,theta);

			theta = M_PI*rng_();
			MatrixType aux2(m.n_row(),m.n_col());
			rotation2d(aux2,1,2,theta);

			MatrixType aux3(m.n_row(),m.n_col());
			theta = M_PI*rng_();
			rotation2d(aux3,2,3,theta);

			MatrixType aux4(m.n_row(),m.n_col());
			theta = M_PI*rng_();
			rotation2d(aux4,3,0,theta);

			m = (aux1 * aux2)*(aux3*aux4);
		}

		void rotation2d(MatrixType& m,size_t x,size_t y,const RealType& theta) const
		{
			for (size_t i=0;i<m.n_row();i++) m(i,i) = 1.0;
			m(x,x) = m(y,y) = cos(theta);
			m(x,y) = sin(theta);
			m(y,x) = -sin(theta);
		}

		const MettsStochasticsType& mettsStochastics_;
		const LeftRightSuperType& lrs_;
		RngType& rng_;
		PsimagLite::ProgressIndicator progress_;
		size_t prevDirection_;
		MatrixType collapseBasis_;
		std::vector<size_t> sitesSeen_;
	};  //class MettsCollapse
} // namespace Dmrg
/*@}*/
#endif //METTS_COLLAPSE_H
