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
#include "Parallel2PointCorrelations.h"
#include "Concurrency.h"
#include "Parallelizer.h"

namespace Dmrg {

template<typename CorrelationsSkeletonType>
class TwoPointCorrelations {
	typedef typename CorrelationsSkeletonType::ObserverHelperType
	ObserverHelperType;
	typedef typename ObserverHelperType::VectorType VectorType ;
	typedef typename ObserverHelperType::VectorWithOffsetType
	VectorWithOffsetType;
	typedef typename ObserverHelperType::BasisWithOperatorsType
	BasisWithOperatorsType ;
	typedef typename VectorType::value_type FieldType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef PsimagLite::Profiling ProfilingType;
	typedef TwoPointCorrelations<CorrelationsSkeletonType> ThisType;

	static SizeType const GROW_RIGHT = CorrelationsSkeletonType::GROW_RIGHT;
	static SizeType const GROW_LEFT = CorrelationsSkeletonType::GROW_LEFT;
	static SizeType const DIAGONAL = CorrelationsSkeletonType::DIAGONAL;
	static SizeType const NON_DIAGONAL = CorrelationsSkeletonType::NON_DIAGONAL;
	static const SizeType EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
	static const SizeType EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;

public:

	typedef typename CorrelationsSkeletonType::SparseMatrixType SparseMatrixType;
	typedef typename ObserverHelperType::MatrixType MatrixType;
	typedef Parallel2PointCorrelations<ThisType> Parallel2PointCorrelationsType;
	typedef typename Parallel2PointCorrelationsType::PairType PairType;

	TwoPointCorrelations(ObserverHelperType& helper,
	                     CorrelationsSkeletonType& skeleton,
	                     bool verbose=false)
	    : helper_(helper),
	      skeleton_(skeleton),
	      verbose_(verbose)
	{}

	void operator()(PsimagLite::Matrix<FieldType>& w,
	                const SparseMatrixType& O1,
	                const SparseMatrixType& O2,
	                int fermionicSign)
	{
		SizeType rows = w.n_row();
		SizeType cols = w.n_col();

		typename PsimagLite::Vector<PairType>::Type pairs;
		for (SizeType i=0;i<rows;i++) {
			for (SizeType j=i;j<cols;j++) {
				if (i>j) continue;
				pairs.push_back(PairType(i,j));
			}
		}

		typedef PsimagLite::Parallelizer<Parallel2PointCorrelationsType> ParallelizerType;
		ParallelizerType threaded2Points(PsimagLite::Concurrency::npthreads,
		                                 PsimagLite::MPI::COMM_WORLD);

		Parallel2PointCorrelationsType helper2Points(w,*this,pairs,O1,O2,fermionicSign);

		threaded2Points.loopCreate(helper2Points);
	}

	// Return the vector: O1 * O2 |psi>
	// where |psi> is the g.s.
	// Note1: O1 is applied to site i and O2 is applied to site j
	// Note2: O1 and O2 operators must commute or anti-commute (set fermionicSign accordingly)
	FieldType calcCorrelation(
	        SizeType i,
	        SizeType j,
	        const SparseMatrixType& O1,
	        const SparseMatrixType& O2,
	        int fermionicSign,
	        SizeType threadId)
	{
		FieldType c = 0;
		if (i==j) {
			c=calcDiagonalCorrelation(i,O1,O2,fermionicSign,threadId);
		} else if (i>j) {
			c= -calcCorrelation_(j,i,O2,O1,fermionicSign,threadId);
		} else {
			c=calcCorrelation_(i,j,O1,O2,fermionicSign,threadId);
		}
		return c;
	}

private:

	SparseMatrixType add(const SparseMatrixType& O1,const SparseMatrixType& O2)
	{
		SizeType n=O1.n_row();
		SparseMatrixType ret(n,n);
		for (SizeType s=0;s<n;s++) for (SizeType t=0;t<n;t++)
			ret(s,t) += O1(s,t)+O2(s,t);
		return ret;
	}

	FieldType calcDiagonalCorrelation(
	        SizeType i,
	        const SparseMatrixType& O1,
	        const SparseMatrixType& O2,
	        int,
	        SizeType threadId)
	{
		SizeType n = O1.row();
		SparseMatrixType O1new=identity(n);

		SparseMatrixType O2new = O1 * O2;
		if (i==0) return calcCorrelation_(0,1,O2new,O1new,1,threadId);
		return calcCorrelation_(i-1,i,O1new,O2new,1,threadId);
	}

	FieldType calcCorrelation_(SizeType i,
	                           SizeType j,
	                           const SparseMatrixType& O1,
	                           const SparseMatrixType& O2,
	                           int fermionicSign,
	                           SizeType threadId)
	{

		if (i>=j) throw PsimagLite::RuntimeError(
		            "Observer::calcCorrelation_(...): i must be smaller than j\n");
		SparseMatrixType O1m,O2m;
		skeleton_.createWithModification(O1m,O1,'n');
		skeleton_.createWithModification(O2m,O2,'n');

		if (j==skeleton_.numberOfSites(threadId)-1) {
			if (i==j-1) {
				helper_.setPointer(threadId,j-2);
				SizeType ni = helper_.leftRightSuper(threadId).left().size()/
				        helper_.leftRightSuper(threadId).right().size();

				SparseMatrixType O1g;
				O1g.makeDiagonal(ni,1.0);

				return skeleton_.bracketRightCorner(O1g,O1m,O2m,fermionicSign,threadId);
			}
			SparseMatrixType O1g;
			skeleton_.growDirectly(O1g,O1m,i,fermionicSign,j-2,true,threadId);
			helper_.setPointer(threadId,j-2);
			return skeleton_.bracketRightCorner(O1g,O2m,fermionicSign,threadId);
		}

		SparseMatrixType O1g,O2g;
		SizeType ns = j-1;

		skeleton_.growDirectly(O1g,O1m,i,fermionicSign,ns,true,threadId);
		skeleton_.dmrgMultiply(O2g,O1g,O2m,fermionicSign,ns,threadId);

		return skeleton_.bracket(O2g,fermionicSign,threadId);
	}

	SparseMatrixType identity(SizeType n)
	{
		SparseMatrixType ret(n,n);
		ret.makeDiagonal(n,1.0);
		return ret;
	}

	//! i can be zero here!!
	void growRecursive(SparseMatrixType& Odest,
	                   const SparseMatrixType& Osrc,
	                   size_t i,
	                   int fermionicSign,
	                   size_t s,
	                   size_t threadId)
	{
		// from 0 --> i
		int nt=i-1;
		if (nt<0) nt=0;

		helper_.setPointer(threadId,s);
		SizeType growOption = skeleton_.growthDirection(s,nt,i,threadId);

		SparseMatrixType Onew(helper_.columns(threadId),helper_.columns(threadId));
		Odest = Onew;
		skeleton_.fluffUp(Odest,Osrc,fermionicSign,growOption,true,threadId);
	}

	ObserverHelperType& helper_;
	CorrelationsSkeletonType& skeleton_;
	bool verbose_;
};  //class TwoPointCorrelations
} // namespace Dmrg

/*@}*/
#endif // TWO_POINT_H
