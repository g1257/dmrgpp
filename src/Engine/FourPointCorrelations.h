/*
Copyright (c) 2008-2012, UT-Battelle, LLC
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
#include "Braket.h"
#include "AnsiColors.h"

namespace Dmrg {
template<typename CorrelationsSkeletonType>
class FourPointCorrelations {

	typedef typename CorrelationsSkeletonType::ObserverHelperType
	ObserverHelperType;
	typedef typename ObserverHelperType::VectorType VectorType ;
	typedef typename ObserverHelperType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ObserverHelperType::BasisWithOperatorsType BasisWithOperatorsType ;

	typedef SizeType IndexType;
	static SizeType const GROW_RIGHT = CorrelationsSkeletonType::GROW_RIGHT;
	typedef typename VectorType::value_type FieldType;
	typedef typename BasisWithOperatorsType::RealType RealType;

public:

	typedef typename CorrelationsSkeletonType::SparseMatrixType SparseMatrixType;
	typedef typename CorrelationsSkeletonType::OperatorType OperatorType;
	typedef typename ObserverHelperType::MatrixType MatrixType;
	typedef PsimagLite::Vector<char>::Type VectorCharType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<SparseMatrixType>::Type VectorSparseMatrixType;
	typedef typename CorrelationsSkeletonType::BraketType BraketType;

	FourPointCorrelations(ObserverHelperType& precomp,
	                      CorrelationsSkeletonType& skeleton,
	                      bool verbose=false)
	    : helper_(precomp),skeleton_(skeleton),verbose_(verbose)
	{
	}

	//! Four-point: these are expensive and uncached!!!
	//! requires i1<i2<i3<i4
	FieldType operator()(SizeType i1,
	                     SizeType i2,
	                     SizeType i3,
	                     SizeType i4,
	                     const BraketType& braket,
	                     SizeType threadId) const
	{
		if (i1>i2 || i3>i4 || i2>i3)
			throw PsimagLite::RuntimeError("calcCorrelation: FourPoint needs ordered points\n");
		if (i1==i2 || i3==i4)
			throw PsimagLite::RuntimeError("calcCorrelation: FourPoint needs distinct points\n");
		if (i1==i3 || i2==i4)
			throw PsimagLite::RuntimeError("calcCorrelation: FourPoint needs distinct points\n");
		if (i2==i3 || i1==i4)
			throw PsimagLite::RuntimeError("calcCorrelation: FourPoint needs distinct points\n");

		SparseMatrixType O2gt;

		firstStage(O2gt,'N',i1,'N',i2,braket,0,1,threadId);

		return secondStage(O2gt,i2,'C',i3,'C',i4,braket,2,3,threadId);
	}

	//! 3-point: these are expensive and uncached!!!
	//! requires i1<i2<i3
	FieldType threePoint(SizeType i1,
	                     SizeType i2,
	                     SizeType i3,
	                     const BraketType& braket,
	                     SizeType threadId) const
	{
		if (i1>i2 || i2>i3)
			throw PsimagLite::RuntimeError("calcCorrelation: FourPoint needs ordered points\n");
		if (i1==i2 || i1==i3 || i2==i3)
			throw PsimagLite::RuntimeError("calcCorrelation: FourPoint needs distinct points\n");

		SparseMatrixType O2gt;

		firstStage(O2gt,'N',i1,'N',i2,braket,0,1,threadId);

		return secondStage(O2gt,i2,'N',i3,braket,2,threadId);
	}

	//! 4-points or more: these are expensive and uncached!!!
	//! requires i0<i1<i2<i3<...<i_{n-1}
	FieldType anyPoint(const BraketType& braket, SizeType threadId) const
	{
		SizeType n = braket.points();
		if (n < 4)
			err("anyPoint: FATAL ERROR: 4 or more points expected\n");

		checkIndicesForStrictOrdering(braket);

		SparseMatrixType O2gt;
		firstStage(O2gt,
		           'N',
		           braket.site(0),
		           'N',
		           braket.site(1),
		           braket,
		           0,
		           1,
		           threadId);

		assert(n > 3);
		SizeType end = n - 2;
		// do the middle
		for (SizeType i = 2; i < end; ++i) {
			SizeType i2 = braket.site(i - 1);
			SparseMatrixType OsoFar;
			middleStage(OsoFar,
			            O2gt,
			            i2,
			            'N',
			            braket.site(i),
			            braket.op(i),
			            braket.op(i - 1).fermionSign,
			            threadId);
			O2gt = OsoFar;
		}

		// handle the last 2 sites
		SizeType i2 = braket.site(n - 3);
		SizeType i3 = braket.site(n - 2);
		SizeType i4 = braket.site(n - 1);

		return secondStage(O2gt,
		                   i2,
		                   'N',
		                   i3,
		                   'N',
		                   i4,
		                   braket,
		                   n - 2,
		                   n - 1,
		                   threadId);
	}

	//! requires i1<i2
	void firstStage(SparseMatrixType& O2gt,
	                char mod1,
	                SizeType i1,
	                char mod2,
	                SizeType i2,
	                const BraketType& braket,
	                SizeType index0,
	                SizeType index1,
	                SizeType threadId) const
	{

		// Take care of modifiers
		SparseMatrixType O1m, O2m;
		skeleton_.createWithModification(O1m,braket.op(index0).data,mod1);
		skeleton_.createWithModification(O2m,braket.op(index1).data,mod2);
		if (verbose_) {
			std::cerr<<"O1m, mod="<<mod1<<"\n";
			std::cerr<<O1m;
			std::cerr<<"O2m, mod="<<mod2<<"\n";
			std::cerr<<O2m;
		}

		// Multiply and grow ("snowball")
		SparseMatrixType O1g,O2g;

		int ns = i2-1;
		if (ns<0) ns = 0;
		skeleton_.growDirectly(O1g,O1m,i1,braket.op(index0).fermionSign,ns,true,threadId);
		skeleton_.dmrgMultiply(O2g,O1g,O2m,braket.op(index1).fermionSign,ns,threadId);

		helper_.setPointer(threadId,ns);
		helper_.transform(O2gt,O2g,threadId);
		if (verbose_) {
			std::cerr<<"O2gt\n";
			std::cerr<<O2gt;
		}
	}

	//! requires i2<i3<i4
	FieldType secondStage(const SparseMatrixType& O2gt,
	                      SizeType i2,
	                      char mod3,
	                      SizeType i3,
	                      char mod4,
	                      SizeType i4,
	                      const BraketType& braket,
	                      SizeType index0,
	                      SizeType index1,
	                      SizeType threadId) const
	{
		// Take care of modifiers
		SparseMatrixType O3m,O4m;
		skeleton_.createWithModification(O3m,braket.op(index0).data,mod3);
		skeleton_.createWithModification(O4m,braket.op(index1).data,mod4);

		int ns = i3-1;
		if (ns<0) ns = 0;
		helper_.setPointer(threadId,ns);
		SparseMatrixType Otmp;
		if (index0 == 0) err("secondStage\n");

		int fermionS = braket.op(index0 - 1).fermionSign;
		growDirectly4p(Otmp,O2gt,i2+1,fermionS,ns,threadId);
		if (verbose_) {
			std::cerr<<"Otmp\n";
			std::cerr<<Otmp;
		}

		SparseMatrixType O3g,O4g;
		if (i4==skeleton_.numberOfSites(threadId)-1) {
			if (i3<i4-1) { // not tested
				std::cerr<<PsimagLite::AnsiColor::red;
				std::cerr<<"WARNING: Correlation code patch might be WRONG!";
				std::cerr<<PsimagLite::AnsiColor::reset<<"\n";
				skeleton_.dmrgMultiply(O3g,Otmp,O3m,braket.op(index0).fermionSign,ns,threadId);
				skeleton_.growDirectly(Otmp,
				                       O3g,
				                       i3,
				                       braket.op(index0).fermionSign,
				                       i4 - 2,
				                       true,
				                       threadId);
				helper_.setPointer(threadId,i4-2);
				return skeleton_.bracketRightCorner(Otmp,
				                                    O4m,
				                                    braket.op(index1).fermionSign,
				                                    threadId);
			}

			helper_.setPointer(threadId,i4-2);
			return skeleton_.bracketRightCorner(Otmp,
			                                    O3m,
			                                    O4m,
			                                    braket.op(index1).fermionSign,
			                                    threadId);
		}

		skeleton_.dmrgMultiply(O3g,Otmp,O3m,braket.op(index0).fermionSign,ns,threadId);
		if (verbose_) {
			std::cerr<<"O3g\n";
			std::cerr<<O3g;
		}

		helper_.setPointer(threadId,ns);

		SparseMatrixType O3gt;
		helper_.transform(O3gt,O3g,threadId);
		if (verbose_) {
			std::cerr<<"O3gt\n";
			std::cerr<<O3gt;
		}

		ns = i4-1;
		if (ns<0) ns = 0;
		helper_.setPointer(threadId,ns);
		growDirectly4p(Otmp,O3gt,i3+1,braket.op(index0).fermionSign,ns,threadId);
		if (verbose_) {
			std::cerr<<"Otmp\n";
			std::cerr<<Otmp;
		}

		skeleton_.dmrgMultiply(O4g,Otmp,O4m,braket.op(index1).fermionSign,ns,threadId);
		return skeleton_.bracket(O4g,braket.op(index1).fermionSign,threadId);
	}

	//! requires i2<i3<i4
	void middleStage(SparseMatrixType& dest,
	                 const SparseMatrixType& OsoFar,
	                 SizeType i2,
	                 char mod3,
	                 SizeType i3,
	                 const OperatorType& Op3,
	                 int fermionS,
	                 SizeType threadId) const
	{
		// Take care of modifiers
		if (i2 > i3)
			throw PsimagLite::RuntimeError("calcCorrelation: FourPoint needs ordered points\n");
		if (i2 == i3)
			throw PsimagLite::RuntimeError("calcCorrelation: FourPoint needs distinct points\n");

		SparseMatrixType O3m;
		skeleton_.createWithModification(O3m,Op3.data,mod3);

		int ns = i3-1;
		if (ns < 0) ns = 0;
		helper_.setPointer(threadId,ns);
		SparseMatrixType Otmp;
		growDirectly4p(Otmp,OsoFar,i2+1,fermionS,ns,threadId);
		if (verbose_) {
			std::cerr<<"Otmp\n";
			std::cerr<<Otmp;
		}

		SparseMatrixType O3g;
		skeleton_.dmrgMultiply(O3g,Otmp,O3m,Op3.fermionSign,ns,threadId);
		if (verbose_) {
			std::cerr<<"O3g\n";
			std::cerr<<O3g;
		}

		helper_.setPointer(threadId,ns);

		helper_.transform(dest,O3g,threadId);
		if (verbose_) {
			std::cerr<<"dest\n";
			std::cerr<<dest;
		}
	}


private:

	//! requires i2<i3
	FieldType secondStage(const SparseMatrixType& O2gt,
	                      SizeType i2,
	                      char mod3,
	                      SizeType i3,
	                      const BraketType& braket,
	                      SizeType index,
	                      SizeType threadId) const
	{
		// Take care of modifiers
		SparseMatrixType O3m;
		skeleton_.createWithModification(O3m,braket.op(index).data,mod3);

		int ns = i3-1;
		if (ns<0) ns = 0;
		SparseMatrixType Otmp;
		if (index == 0) err("secondStage\n");
		int fermionS = braket.op(index - 1).fermionSign;
		growDirectly4p(Otmp,O2gt,i2+1,fermionS,ns,threadId);
		if (verbose_) {
			std::cerr<<"Otmp\n";
			std::cerr<<Otmp;
		}

		if (i3 == skeleton_.numberOfSites(threadId)-1) {
			helper_.setPointer(threadId,i3-2);
			return skeleton_.bracketRightCorner(Otmp,
			                                    O3m,
			                                    braket.op(index).fermionSign,
			                                    threadId);
		}

		helper_.setPointer(threadId,ns);
		SparseMatrixType O3g;
		skeleton_.dmrgMultiply(O3g,Otmp,O3m,braket.op(index).fermionSign,ns,threadId);
		if (verbose_) {
			std::cerr<<"O3g\n";
			std::cerr<<O3g;
		}

		helper_.setPointer(threadId,ns);
		return skeleton_.bracket(O3g,braket.op(index).fermionSign,threadId);
	}

	//! i can be zero here!!
	void growDirectly4p(SparseMatrixType& Odest,
	                    const SparseMatrixType& Osrc,
	                    SizeType i,
	                    int fermionicSign,
	                    SizeType ns,
	                    SizeType threadId) const
	{
		Odest =Osrc;

		// from 0 --> i
		int nt=i-1;
		if (nt<0) nt=0;

		for (SizeType s=nt;s<ns;s++) {
			helper_.setPointer(threadId,s);
			int growOption = GROW_RIGHT;

			SparseMatrixType Onew(helper_.columns(threadId),helper_.columns(threadId));
			skeleton_.fluffUp(Onew,Odest,fermionicSign,growOption,true,threadId);
			Odest = Onew;

		}
	}

	void checkIndicesForStrictOrdering(const BraketType& braket) const
	{
		if (braket.points() < 2) return;

		bool flag = true;
		SizeType prev = braket.site(0);
		for (SizeType i = 1; i < braket.points(); ++i) {
			if (braket.site(i) <= prev) {
				flag = false;
				break;
			}
		}

		if (flag) return;
		throw PsimagLite::RuntimeError("AnyPoint: Point must be strictly ordered\n");
	}

	ObserverHelperType& helper_; // <-- NB: not the owner
	CorrelationsSkeletonType& skeleton_; // <-- NB: not the owner
	bool verbose_;
};  //class FourPointCorrelations
} // namespace Dmrg

/*@}*/
#endif

