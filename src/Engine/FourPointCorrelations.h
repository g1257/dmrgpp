// BEGIN LICENSE BLOCK
/*
Copyright © 2008 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 1.0.0]
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

/*! \file FourPointCorrelations.h
 *
 *  A class to perform post-processing calculation of 
 *  4-point correlations of the form A_i B_j C_k D_l
 *  with i\<j\<k \< l
 *
 */
#ifndef FOURPOINT_C_H
#define FOURPOINT_C_H

#include "CrsMatrix.h"
#include "CorrelationsSkeleton.h"

namespace Dmrg {
	template<typename FieldType,typename RealType,typename IoType,typename MatrixType,template<typename> class VectorTemplate>
	class FourPointCorrelations {
		//typedef typename MatrixType::value_type FieldType;
		typedef size_t IndexType;
		typedef Precomputed<FieldType,RealType,IoType,MatrixType,VectorTemplate> PrecomputedType;
		typedef CorrelationsSkeleton<FieldType,RealType,IoType,MatrixType,VectorTemplate> CorrelationsSkeletonType;
		static size_t const GROW_RIGHT = CorrelationsSkeletonType::GROW_RIGHT;
		
	public:
		FourPointCorrelations(PrecomputedType& precomp,CorrelationsSkeletonType& skeleton,bool verbose=false)
			: precomp_(precomp),skeleton_(skeleton),verbose_(verbose)
		{
		}
				
		//! Four-point: these are expensive and uncached!!!
		//! requires i<=j<=k<=l
		FieldType operator()(
				char mod1,size_t i1,const MatrixType& O1,
				char mod2,size_t i2,const MatrixType& O2,
				char mod3,size_t i3,const MatrixType& O3,
				char mod4,size_t i4,const MatrixType& O4,
				int fermionicSign)
		{
			if (i1>i2 || i3>i4)
				throw std::runtime_error("calcCorrelation: FourPoint needs ordered points\n");
			if (i1==i2 || i3==i4)
				throw std::runtime_error("calcCorrelation: FourPoint needs distinct points\n");
			
			// Take care of modifiers
			MatrixType O1m, O2m,O3m,O4m;
			skeleton_.createWithModification(O1m,O1,mod1);
			skeleton_.createWithModification(O2m,O2,mod2);
			skeleton_.createWithModification(O3m,O3,mod3);
			skeleton_.createWithModification(O4m,O4,mod4);
			if (verbose_) {
				std::cerr<<"O1m, mod="<<mod1<<"\n";
				std::cerr<<O1m;
				std::cerr<<"O2m, mod="<<mod2<<"\n";
				std::cerr<<O2m;
			}
			
			// Multiply and grow ("snowball")
			MatrixType O1g,O2g,O3g,O4g;

			int ns = i2-1;
			if (ns<0) ns = 0;
			skeleton_.growDirectly(O1g,O1m,i1,fermionicSign,ns);
			skeleton_.dmrgMultiply(O2g,O1g,O2m,fermionicSign,ns);
			
			//std::cerr<<"Result="<<bracket(O2g)<<"\n";
			
			//ns++;
			precomp_.setPointer(ns);
			size_t trunc = precomp_.transform().n_col();
			MatrixType O2gt(trunc,trunc);
			precomp_.transform(O2gt,O2g);
			if (verbose_) {
				std::cerr<<"O2gt\n";
				std::cerr<<O2gt;
			}
			
			ns = i3-1;
			if (ns<0) ns = 0;
			precomp_.setPointer(ns);
			MatrixType Otmp;
			growDirectly4p(Otmp,O2gt,i2+1,fermionicSign,ns);
			if (verbose_) {
				std::cerr<<"Otmp\n";
				std::cerr<<Otmp;
			}
			skeleton_.dmrgMultiply(O3g,Otmp,O3m,fermionicSign,ns);
			if (verbose_) {
				std::cerr<<"O3g\n";
				std::cerr<<O3g;
			}
			precomp_.setPointer(ns);
			trunc = precomp_.transform().n_col();
			MatrixType O3gt(trunc,trunc);
			precomp_.transform(O3gt,O3g);
			if (verbose_) {
				std::cerr<<"O3gt\n";
				std::cerr<<O3gt;
			}
			ns = i4-1;
			if (ns<0) ns = 0;
			precomp_.setPointer(ns);
			growDirectly4p(Otmp,O3gt,i3+1,fermionicSign,ns);
			if (verbose_) {
				std::cerr<<"Otmp\n";
				std::cerr<<Otmp;
			}
			skeleton_.dmrgMultiply(O4g,Otmp,O4m,fermionicSign,ns);
			
			return skeleton_.bracket(O4g);
		}
			
	private:
		//! i can be zero here!!
		void growDirectly4p(MatrixType& Odest,const MatrixType& Osrc,size_t i,int fermionicSign,size_t ns)
		{
			Odest =Osrc;
			std::vector<int> signs;
			// from 0 --> i
			int nt=i-1;
			if (nt<0) nt=0;
			
			for (size_t s=nt;s<ns;s++) {
				// set appropriate privates which are:
				// SpermutationInverse_(s) and transform_(s)
				/*io_.rewind();
				calcSpermutation(s);
				//std::cerr<<"*****************======="<<transform_.n_row()<<"\n";
				io_.readMatrix(transform_,"#TRANSFORM_sites",s);*/
				precomp_.setPointer(s);
				//std::cerr<<"%%%%%%%%%%%%%%%%%======="<<transform_.n_row()<<"\n";
				int growOption = GROW_RIGHT;
				//if (i==1 && s==0) growOption = GROW_LEFT;// <-- not needed since nt>0
				
				/*io_.rewind();
				io_.read(electrons_,"#ELECTRONS_sites=",s);*/
				skeleton_.createSigns(signs,fermionicSign);
				MatrixType Onew(precomp_.transform().n_col(),precomp_.transform().n_col());
				skeleton_.fluffUp(Onew,Odest,signs,growOption);
				Odest = Onew;
				
			}
		}
				
// 		void getCombined(MatrixType& ret,MatrixType& O1,MatrixType& O2,MatrixType& O3)
// 		{
// 			int nBig = O1.n_row();
// 			int nSmall = nBig;
// 			MatrixType fmTmp(nSmall,nBig);
// 			FieldType alpha=1.0,beta=0.0;
// 			ret.resize(nSmall,nSmall);
// 			
// 			if (isZero(O2)) throw std::runtime_error("O2 is zero\n");
// 			psimag::BLAS::GEMM('N','N',nBig,nSmall,nBig,alpha,
// 					   &(O1(0,0)),nBig,&(O2(0,0)),nBig,beta,&(fmTmp(0,0)),nBig);
// 			if (isZero(fmTmp)) throw std::runtime_error("fmTmp is zero\n");
// 			psimag::BLAS::GEMM('N','N',nSmall,nSmall,nBig,alpha,
// 					   &(fmTmp(0,0)),nBig,&(O3(0,0)),nBig,beta,&(ret(0,0)),nSmall);
// 			if (isZero(O3)) throw std::runtime_error("O3 is zero\n");
// 			if (isZero(ret)) throw std::runtime_error("ret is zero\n");
// 			
// 		}
		
		
		
		
		
		PrecomputedType& precomp_; // <-- NB: not the owner
		CorrelationsSkeletonType& skeleton_; // <-- NB: not the owner
		bool verbose_;
	};  //class FourPointCorrelations
} // namespace Dmrg

/*@}*/
#endif
