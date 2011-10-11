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
	template<typename VectorWithOffsetType,typename MettsStochasticsType>
	class MettsCollapse  {
		typedef typename VectorWithOffsetType::VectorType VectorType;
		typedef typename MettsStochasticsType::PairType PairType;
		typedef typename MettsStochasticsType::LeftRightSuperType LeftRightSuperType;
		typedef typename MettsStochasticsType::RealType RealType;
		typedef PsimagLite::Matrix<RealType> MatrixType;
		enum {EXPAND_ENVIRON=ProgramGlobals::EXPAND_ENVIRON,
		      EXPAND_SYSTEM=ProgramGlobals::EXPAND_SYSTEM};
	public:
		typedef PsimagLite::PackIndices PackIndicesType;

		MettsCollapse(const MettsStochasticsType& mettsStochastics,
		              const LeftRightSuperType& lrs)
		: mettsStochastics_(mettsStochastics),
		  lrs_(lrs),
		  hilbertSize_(mettsStochastics.hilbertSizePerSite())
		{}

		void operator()(VectorWithOffsetType& dest2,
		                const VectorWithOffsetType& src2,
		                size_t site,
		                size_t direction)
		{
			if (dest2.size()==0) {
				dest2 =  src2;
			}
			std::vector<RealType> p(hilbertSize_,0);
			probability(p,dest2,direction);
			RealType sum = 0;
			for (size_t i=0;i<p.size();i++)
				sum += p[i];
			assert(fabs(sum-1.0)<1e-6);
			VectorWithOffsetType dest;
			size_t indexFixed = mettsStochastics_.chooseRandomState(p);
			collapseVector(dest,dest2,direction,indexFixed);
			RealType x = std::norm(dest);
			assert(x>1e-6);
			dest2 = (1.0/x) * dest;
		}

	private:
		void collapseVector(VectorWithOffsetType& dest,
		                    const VectorWithOffsetType& src,
		                    size_t direction,
		                    size_t indexFixed) const
		{
			assert(src.sectors()==1);
			dest = src;
			for (size_t ii=0;ii<src.sectors();ii++) {
				size_t i0 = src.sector(ii);
				VectorType vdest,vsrc;
				src.extract(vsrc,i0);
				collapseVector(vdest,vsrc,direction,i0,indexFixed);
				dest.setDataInSector(vdest,i0);
			}
			//assert(std::norm(dest)>1e-6);
		}

		void collapseVector(VectorType& w,
		                    const VectorType& v,
		                    size_t direction,
		                    size_t m,
		                    size_t indexFixed) const
		{
			int offset = lrs_.super().partition(m);
			int total = lrs_.super().partition(m+1) - offset;
			
			size_t nk = hilbertSize_;
			size_t ns = lrs_.left().size();
			PackIndicesType packSuper(ns);
			PackIndicesType packLeft(ns/nk);
			PackIndicesType packRight(nk);
			w.resize(total);
			for (size_t i=0;i<size_t(total);i++) {
				w[i] = 0;
				size_t alpha,beta;
				packSuper.unpack(alpha,beta,lrs_.super().permutation(i+offset));
				size_t alpha0,alpha1;
				packLeft.unpack(alpha0,alpha1,alpha);
				size_t beta0,beta1;
				packRight.unpack(beta0,beta1,beta);
				if (direction==EXPAND_SYSTEM && alpha1!=indexFixed) continue;
				if (direction==EXPAND_ENVIRON && beta0 != indexFixed) continue;
				w[i] = v[i];
			}
		}

		void probability(std::vector<RealType>& p,
		                 const VectorWithOffsetType& src,
						 size_t direction) const
		{
			RealType sum = 0;
			for(size_t alpha=0;alpha<hilbertSize_;alpha++) {
				VectorWithOffsetType dest;
				collapseVector(dest,src,direction,alpha);
				RealType x = std::norm(dest);
				sum += x*x;
				p[alpha] = x*x;
			}
			assert(fabs(sum)>1e-6);
			for(size_t alpha=0;alpha<p.size();alpha++) p[alpha] /= sum;
		}

		const MettsStochasticsType& mettsStochastics_;
		const LeftRightSuperType& lrs_;
		size_t hilbertSize_;
	};  //class MettsCollapse
} // namespace Dmrg
/*@}*/
#endif //METTS_COLLAPSE_H
