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

/*! \file DmrgObserve.h
 *
 *  A class to perform post-processing calculation of observables
 *
 */
#ifndef DMRG_OBSERVE_H
#define DMRG_OBSERVE_H
#include "Precomputed.h"
#include "CrsMatrix.h"
#include "CorrelationsSkeleton.h"
#include "FourPointCorrelations.h"

namespace Dmrg {
	
	template<typename RealType,typename FieldType,typename IoType,typename ConcurrencyType>
	class Observe {
		typedef size_t IndexType;
		typedef SparseVector<FieldType> VectorType;
		typedef psimag::Matrix<FieldType> MatrixType;
		typedef Precomputed<RealType,FieldType,IoType,MatrixType,SparseVector> PrecomputedType;
		typedef CorrelationsSkeleton<RealType,FieldType,IoType,MatrixType,SparseVector> CorrelationsSkeletonType;
		typedef FourPointCorrelations<RealType,FieldType,IoType,MatrixType,SparseVector>  FourPointCorrelationsType;
		
		static size_t const GROW_RIGHT = CorrelationsSkeletonType::GROW_RIGHT;
		static size_t const GROW_LEFT = CorrelationsSkeletonType::GROW_LEFT;
		static size_t const DIAGONAL = CorrelationsSkeletonType::DIAGONAL;
		static size_t const NON_DIAGONAL = CorrelationsSkeletonType::NON_DIAGONAL;
		
	public:
		Observe(const std::string& filename,size_t n,size_t n1,size_t stepTimes,ConcurrencyType& concurrency,bool verbose=false)
		: precomp_(filename,2*n,stepTimes,verbose),halfLatticeSize_(n),
			    oneSiteHilbertSize_(n1),skeleton_(precomp_),fourpoint_(precomp_,skeleton_),concurrency_(concurrency),
				verbose_(verbose)
		{}
		
		psimag::Matrix<FieldType> correlations(size_t n,const MatrixType& O1,const MatrixType& O2,int fermionicSign,
						      size_t n1=0,size_t nf=0)
		{
			if (nf==0) nf = 2*n;
			if (n1==0) n1 = n;
			/*clearCache(n1, nf);
			precomputeGrowth(O1,fermionicSign,n1,nf);*/
			initCache(O1,n1,nf,fermionicSign);
			psimag::Matrix<FieldType> w(n,nf);
			for (size_t i=0;i<n1;i++) {
				concurrency_.loopCreate(nf-1);
				std::vector<FieldType> v(nf-1);
				size_t j = i;
				while(concurrency_.loop(j)) {
				//for (size_t j=i;j<nf-1;j++) {
					//std::cerr<<"About to do i="<<i<<" and j="<<j<<"\n";
					//try {
						v[j]  = calcCorrelation(i,j,O1,O2,fermionicSign);
						if (verbose_) std::cerr<<"Result for i="<<i<<" and j="<<j<<" is "<<v[j]<<"\n";
					//} catch (std::exception& e) {
					//	std::cerr<<"Result for i="<<i<<" and j="<<j<<" exception caught: "<<e.what()<<"\n";
					//}
				}
				concurrency_.gather(v);
				for (j=i;j<nf-1;j++) w(i,j) = v[j];
			}
			return w;
		}

		void initCache(const MatrixType& O1,size_t n1, size_t nf,int fermionicSign)
		{
			clearCache(n1, nf);
                        precomputeGrowth(O1,fermionicSign,n1,nf);
		}
	
		// Return the vector: O1 * O2 |psi>
		// where |psi> is the g.s. 
		// Note1: O1 is applied to site i and O2 is applied to site j
		// Note2: O1 and O2 operators must commute or anti-commute (set fermionicSign accordingly)
		FieldType calcCorrelation(
					size_t i,
					size_t j,
					const MatrixType& O1,
					const MatrixType& O2,
					int fermionicSign)
		{
			FieldType c = 0;
			if (i==j) {
				c=calcDiagonalCorrelation(i,O1,O2,fermionicSign);
			} else if (i>j) {
				c= -calcCorrelation_(j,i,O2,O1,fermionicSign);
			} else {
				c=calcCorrelation_(i,j,O1,O2,fermionicSign);
			}
			return c;
		}
		
		MatrixType multiplyTranspose(const MatrixType& O1,const MatrixType& O2)
		{
			size_t n=O1.n_row();
			MatrixType ret(n,n);
			for (size_t s=0;s<n;s++) 
				for (size_t t=0;t<n;t++) 
					for (size_t w=0;w<n;w++) 
						ret(s,t) += std::conj(O1(s,w))*O2(w,t);
			return ret;
		}
		
		MatrixType add(const MatrixType& O1,const MatrixType& O2)
		{
			size_t n=O1.n_row();
			MatrixType ret(n,n);
			for (size_t s=0;s<n;s++) for (size_t t=0;t<n;t++)  ret(s,t) += O1(s,t)+O2(s,t);
			return ret;
		}

		FieldType fourPoint(
				char mod1,size_t i1,const MatrixType& O1,
				char mod2,size_t i2,const MatrixType& O2,
				char mod3,size_t i3,const MatrixType& O3,
				char mod4,size_t i4,const MatrixType& O4,
				int fermionicSign)
		{
			return fourpoint_(mod1,i1,O1,mod2,i2,O2,mod3,i3,O3,mod4,i4,O4,fermionicSign);
		}

		FieldType onePoint(
					size_t i,
					const MatrixType& O1,
					int fermionicSign)
		{
			
			size_t n = O1.n_row();
			MatrixType Oid=identity(n);
			if (i==0) return calcCorrelation_(0,1,O1,Oid,1,NON_DIAGONAL,PrecomputedType::USETIMEVECTOR);
			return calcCorrelation_(i-1,i,Oid,O1,1,DIAGONAL,PrecomputedType::USETIMEVECTOR);
		}
	
	private:
		
		FieldType calcDiagonalCorrelation(
					size_t i,
					const MatrixType& O1,
					const MatrixType& O2,
					int fermionicSign)
		{
			
			size_t n = O1.n_row();
			MatrixType O1new=identity(n);
			//for (size_t s=0;s<n;s++)  O1new(s,s)=static_cast<FieldType>(1.0);
			MatrixType O2new=multiplyTranspose(O1,O2);
			//for (size_t s=0;s<n;s++) for (size_t t=0;t<n;t++) for (size_t w=0;w<n;w++) O2new(s,t) += conj(O1(w,s))*O2(w,t);
			if (i==0) return calcCorrelation_(0,1,O2new,O1new,1);
			return calcCorrelation_(i-1,i,O1new,O2new,1,DIAGONAL);
		}

		//! FIXME : make sure that i less than j, which is assumed here
		FieldType calcCorrelation_(
					size_t i,
					size_t j,
					const MatrixType& O1,
					const MatrixType& O2,
					int fermionicSign,
					size_t isDiagonal=NON_DIAGONAL,
					size_t useTimeVector=PrecomputedType::NOTIMEVECTOR)
		{
			MatrixType O1g,O2g,O1m,O2m;
			skeleton_.createWithModification(O1m,O1,'n');
			skeleton_.createWithModification(O2m,O2,'n');
			
			int ns = j-1;
			if (ns<0) ns = 0;
			skeleton_.growDirectly(O1g,O1m,i,fermionicSign,ns);
			skeleton_.dmrgMultiply(O2g,O1g,O2m,fermionicSign,ns);
			
			return skeleton_.bracket(O2g,useTimeVector);
		}

		MatrixType identity(size_t n)
		{
			MatrixType ret(n,n);
			for (size_t s=0;s<n;s++)  ret(s,s)=static_cast<FieldType>(1.0);	
			return ret;
		}

		void multiply(MatrixType& O1,FieldType x)
		{
			for (size_t i=0;i<O1.n_row();i++) for (size_t j=0;j<O1.n_col();j++) O1(i,j) *= x;
		}

		MatrixType sustract(const MatrixType& O1,const MatrixType& O2)
		{
			size_t n=O1.n_row();
			MatrixType ret(n,n);
			for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) ret(i,j) = O1(i,j) - O2(i,j);
			return ret;
		}
		
		void clearCache(size_t  ns,size_t nf)
		{
			growCached_.clear();
			grownOperators_.clear();
			std::vector<MatrixType> v;
			for (size_t i=0;i<nf;i++) {
				MatrixType tmp(1,1);
				v.push_back(tmp);
			}
			for (size_t i=0;i<ns;i++) 
				grownOperators_.push_back(v);
		}

		//! i can be zero here!!
		const MatrixType* grow(const MatrixType& Osrc,size_t i,int fermionicSign,size_t ns,size_t isDiagonal)
		{
			if (isDiagonal==DIAGONAL) {
				static MatrixType Ox;
				skeleton_.growDirectly(Ox,Osrc,i,fermionicSign,ns);
				return &Ox;
			}	
			int nt=i-1;
			if (nt<0) nt=0;
			return &(grownOperators_[i][ns-nt-1]);
			//std::cerr<<"done grow for i="<<i<<" and ns="<<ns<<"\n";
		}

		//! i can be zero here!!
		void precomputeGrowth(const MatrixType& Osrc,int fermionicSign,size_t ns,size_t nfinal)	 
		{
			//MatrixType tmp;
			for (size_t i=0;i<ns;i++) {
				int nt=i-1;
				if (nt<0) nt=0;
				MatrixType Oinc = Osrc;
				if (verbose_) std::cerr<<"Precomputing "<<i<<" out of "<<(ns-1)<<"\n";
				for (size_t s=nt+1;s<nfinal;s++) {
					if (verbose_) std::cerr<<"\tPrecomputing "<<s<<" out of "<<(nfinal-1)<<"\n";
					growRecursive(grownOperators_[i][s-nt-1],Oinc,i,fermionicSign,s-1);
					Oinc = grownOperators_[i][s-nt-1];
					//growDirectly(grownOperators_[i][s-nt-1],Osrc,i,fermionicSign,s);
					//if (tmp!=grownOperators_[i][s-nt-1]) throw std::runtime_error("Not equal\n");
				}
			}
			std::cerr<<"precomputeGrowth done\n";
		}

		//! i can be zero here!!
		void growRecursive(MatrixType& Odest,const MatrixType& Osrc,size_t i,int fermionicSign,size_t s)
		{
			std::vector<int> signs;
			// from 0 --> i
			int nt=i-1;
			if (nt<0) nt=0;
			
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
			if (s==size_t(nt)) {
				growOption = GROW_LEFT;
				if (i==0) growOption = GROW_RIGHT;
			}
			/* io_.rewind();
			io_.read(electrons_,"#ELECTRONS_sites=",s);*/
			skeleton_.createSigns(signs,fermionicSign);
			MatrixType Onew(precomp_.transform().n_col(),precomp_.transform().n_col());
			Odest = Onew;
			skeleton_.fluffUp(Odest,Osrc,signs,growOption);
		}

		PrecomputedType precomp_;
		size_t halfLatticeSize_;
		size_t oneSiteHilbertSize_;
		CorrelationsSkeletonType skeleton_;
		FourPointCorrelationsType fourpoint_;
		ConcurrencyType& concurrency_;
		bool verbose_;
		std::vector<size_t> growCached_;
		std::vector<std::vector<MatrixType> > grownOperators_;
		
	};  //class Observe
} // namespace Dmrg

/*@}*/
#endif
