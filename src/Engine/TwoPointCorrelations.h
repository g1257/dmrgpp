// BEGIN LICENSE BLOCK
/*
Copyright (c) 2008 , UT-Battelle, LLC
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

/*! \file TwoPointCorrelations.h
 *
 *  A class to perform post-processing calculation of TwoPointCorrelations
 *  <state1 | A_i B_j |state2>
 *
 */
#ifndef TWO_POINT_H
#define TWO_POINT_H
#include "CrsMatrix.h"
#include "VectorWithOffsets.h" // for operator*
#include "VectorWithOffset.h" // for operator*
#include "Profiling.h"

namespace Dmrg {
	
	template<typename CorrelationsSkeletonType,
		typename ConcurrencyType>
	class TwoPointCorrelations {
		typedef typename CorrelationsSkeletonType::ObserverHelperType
			ObserverHelperType;
		typedef typename ObserverHelperType::MatrixType MatrixType;
		typedef typename ObserverHelperType::VectorType VectorType ;
		typedef typename ObserverHelperType::VectorWithOffsetType
			VectorWithOffsetType;
		typedef typename ObserverHelperType::BasisWithOperatorsType
			BasisWithOperatorsType ;
		typedef typename VectorType::value_type FieldType;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef PsimagLite::Profiling ProfilingType;

		static size_t const GROW_RIGHT = CorrelationsSkeletonType::GROW_RIGHT;
		static size_t const GROW_LEFT = CorrelationsSkeletonType::GROW_LEFT;
		static size_t const DIAGONAL = CorrelationsSkeletonType::DIAGONAL;
		static size_t const NON_DIAGONAL = CorrelationsSkeletonType::NON_DIAGONAL;
		static const size_t EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
		static const size_t EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;

	public:
		TwoPointCorrelations(
				ObserverHelperType& helper,
				CorrelationsSkeletonType& skeleton,
				ConcurrencyType& concurrency,
				bool verbose=false)
		: helper_(helper),
		  skeleton_(skeleton),
		  concurrency_(concurrency),
		  verbose_(verbose)
		{}

		PsimagLite::Matrix<FieldType> operator()(
				const MatrixType& O1,
				const MatrixType& O2,
				int fermionicSign,
				size_t rows,
				size_t cols)
		{
			PsimagLite::Matrix<FieldType> w(rows,cols);

			initCache(O1,rows,cols,fermionicSign);

			for (size_t i=0;i<rows;i++) {
				concurrency_.loopCreate(cols);
				std::vector<FieldType> v(cols);
				size_t j = i;
				//ProfilingType profile("correlations loop i=" + utils::ttos(i));
				while(concurrency_.loop(j)) {
					v[j]  = calcCorrelation(i,j,O1,O2,fermionicSign);
					if (verbose_) {
						std::cerr<<"Result for i="<<i;
						std::cerr<<" and j="<<j<<" is "<<v[j]<<"\n";
					}
				}
				concurrency_.gather(v);
				for (j=i;j<v.size();j++) w(i,j) = v[j];
			}

			return w;
		}

	private:

		void initCache(const MatrixType& O1,size_t n1, size_t nf,int fermionicSign)
		{
			clearCache(n1, nf);
            precomputeGrowth(O1,fermionicSign,n1,nf-1);
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

		MatrixType multiplyTranspose(
				const MatrixType& O1,
				const MatrixType& O2)
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
			for (size_t s=0;s<n;s++) for (size_t t=0;t<n;t++)
				ret(s,t) += O1(s,t)+O2(s,t);
			return ret;
		}
			
		FieldType calcDiagonalCorrelation(
					size_t i,
					const MatrixType& O1,
					const MatrixType& O2,
					int fermionicSign)
		{
			size_t n = O1.n_row();
			MatrixType O1new=identity(n);

			MatrixType O2new=multiplyTranspose(O1,O2);
			if (i==0) return calcCorrelation_(0,1,O2new,O1new,1);
			return calcCorrelation_(i-1,i,O1new,O2new,1);
		}

		FieldType calcCorrelation_(
					size_t i,
					size_t j,
					const MatrixType& O1,
					const MatrixType& O2,
					int fermionicSign)
		{
			
			if (i>=j) throw std::runtime_error(
					"Observer::calcCorrelation_(...): i must be smaller than j\n");
			MatrixType O1m,O2m;
			skeleton_.createWithModification(O1m,O1,'n');
			skeleton_.createWithModification(O2m,O2,'n');

			if (j==skeleton_.numberOfSites()-1) {
				if (i==j-1) {
					helper_.setPointer(j-2);
					size_t ni = helper_.basisS().size()/helper_.basisE().size();
					MatrixType O1g(ni,ni);
					for (size_t x=0;x<O1g.n_row();x++) O1g(x,x) = 1.0;
					return skeleton_.bracketRightCorner(O1g,O1m,O2m,fermionicSign);
				}
				MatrixType O1g,O2g;
				skeleton_.growDirectly(O1g,O1m,i,fermionicSign,j-2);
				helper_.setPointer(j-2);
				return skeleton_.bracketRightCorner(O1g,O2m,fermionicSign);
			}
			MatrixType O1g,O2g;
			size_t ns = j-1;
			skeleton_.growDirectly(O1g,O1m,i,fermionicSign,ns);
			skeleton_.dmrgMultiply(O2g,O1g,O2m,fermionicSign,ns);

			return skeleton_.bracket(O2g,fermionicSign);
		}

		MatrixType identity(size_t n)
		{
			MatrixType ret(n,n);
			for (size_t s=0;s<n;s++)  ret(s,s)=static_cast<RealType>(1.0);
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
		const MatrixType* grow(
				const MatrixType& Osrc,
				size_t i,
				int fermionicSign,
				size_t ns,
				size_t isDiagonal)
		{
			if (isDiagonal==DIAGONAL) {
				static MatrixType Ox;
				skeleton_.growDirectly(Ox,Osrc,i,fermionicSign,ns);
				return &Ox;
			}	
			int nt=i-1;
			if (nt<0) nt=0;
			return &(grownOperators_[i][ns-nt-1]);
		}

		//! i can be zero here!!
		void precomputeGrowth(
				const MatrixType& Osrc,
				int fermionicSign,
				size_t ns,
				size_t nfinal)
		{
			for (size_t i=0;i<ns;i++) {
				int nt=i-1;
				if (nt<0) nt=0;
				MatrixType Oinc = Osrc;
				if (verbose_)
					std::cerr<<"Precomputing "<<i<<" out of "<<(ns-1)<<"\n";
				for (size_t s=nt+1;s<nfinal;s++) {
					if (verbose_) {
						std::cerr<<"\tPrecomputing "<<s;
						std::cerr<<" out of "<<(nfinal-1)<<"\n";
					}
					growRecursive(grownOperators_[i][s-nt-1],Oinc,i,
							fermionicSign,s-1);
					Oinc = grownOperators_[i][s-nt-1];
				}
			}
			if (verbose_) std::cerr<<"precomputeGrowth done\n";
		}

		//! i can be zero here!!
		void growRecursive(MatrixType& Odest,
				const MatrixType& Osrc,
				size_t i,
				int fermionicSign,
				size_t s)
		{
			std::vector<int> signs;
			// from 0 --> i
			int nt=i-1;
			if (nt<0) nt=0;
			
			helper_.setPointer(s);
			size_t growOption = skeleton_.growthDirection(s,nt,i);

			MatrixType Onew(helper_.columns(),helper_.columns());
			Odest = Onew;
			skeleton_.fluffUp(Odest,Osrc,fermionicSign,growOption);
		}

		ObserverHelperType& helper_;
		CorrelationsSkeletonType& skeleton_;
		ConcurrencyType& concurrency_;
		bool verbose_;
		std::vector<size_t> growCached_;
		std::vector<std::vector<MatrixType> > grownOperators_;
		
	};  //class TwoPointCorrelations
} // namespace Dmrg

/*@}*/
#endif // TWO_POINT_H
