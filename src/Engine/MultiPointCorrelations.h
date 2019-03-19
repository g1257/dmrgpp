/*
Copyright (c) 2013, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5..0]
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

/*! \file MultiPointCorrelations.h
 *
 *
 *
 */
#ifndef MULTI_POINT_CORRELATIONS_H
#define MULTI_POINT_CORRELATIONS_H
#include "CrsMatrix.h"
#include "VectorWithOffsets.h" // for operator*
#include "VectorWithOffset.h" // for operator*

namespace Dmrg {

template<typename CorrelationsSkeletonType>
class MultiPointCorrelations {
	typedef typename CorrelationsSkeletonType::ObserverHelperType
	ObserverHelperType;
	typedef typename ObserverHelperType::VectorType VectorType ;
	typedef typename ObserverHelperType::VectorWithOffsetType
	VectorWithOffsetType;
	typedef typename ObserverHelperType::BasisWithOperatorsType
	BasisWithOperatorsType ;
	typedef typename VectorType::value_type FieldType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef MultiPointCorrelations<CorrelationsSkeletonType> ThisType;
	typedef typename CorrelationsSkeletonType::SparseMatrixType SparseMatrixType;

public:

	typedef typename ObserverHelperType::MatrixType MatrixType;

	MultiPointCorrelations(SizeType nthreads,
	                       ObserverHelperType& helper,
	                       CorrelationsSkeletonType& skeleton,
	                       bool verbose=false)
	    : nthreads_(nthreads),
	      helper_(helper),
	      skeleton_(skeleton),
	      verbose_(verbose)
	{}

	template<typename VectorLikeType>
	typename PsimagLite::EnableIf
	<PsimagLite::IsVectorLike<VectorLikeType>::True,void>::Type
	operator()(VectorLikeType& result,
	           const SparseMatrixType& O,
	           SizeType rows,
	           SizeType cols)
	{
		assert(rows == cols);
		size_t threadId = 0;
		result.resize(rows);
		SparseMatrixType Og;
		SparseMatrixType identity(O.rows(),O.cols());

		identity.makeDiagonal(O.rows(),1.0);

		size_t rowsOver2 = static_cast<size_t>(rows/2);

		for (SizeType i=0;i<rowsOver2;i++) {
			result[i] = calcCorrelation_(Og,i,O,identity,threadId);
		}
		for (SizeType i=rowsOver2; i<rows; i++) result[i]=0;
	}

private:

	// from i to i+1
	FieldType calcCorrelation_(SparseMatrixType& O2gt,
	                           SizeType i,
	                           const SparseMatrixType& O,
	                           const SparseMatrixType& identity,
	                           SizeType threadId)
	{

		if (i>=skeleton_.numberOfSites(threadId)-1)
			throw PsimagLite::RuntimeError("calcCorrelation: i must be < sites-1\n");
		ProgramGlobals::FermionOrBosonEnum fermionicSign =
		        ProgramGlobals::FermionOrBosonEnum::BOSON;

		SizeType ns = i;
		SparseMatrixType O2g;
		if (i==0) {
			skeleton_.growDirectly(O2gt,O,i,fermionicSign,ns,true,threadId);
			skeleton_.dmrgMultiply(O2g,O2gt,identity,fermionicSign,ns,threadId);
			FieldType ret = skeleton_.bracket(O2g,fermionicSign,threadId);
			return ret;
		}

		//			if (i==5) {
		skeleton_.dmrgMultiply(O2g,O2gt,O,fermionicSign,ns-1,threadId);
		//			} else {
		//				skeleton_.dmrgMultiply(O2g,O2gt,identity,fermionicSign,ns-1,threadId);
		//			}
		O2gt.clear();
		FieldType ret = skeleton_.bracket(O2g,fermionicSign,threadId);
		helper_.setPointer(threadId,ns-1);
		helper_.transform(O2gt,O2g,threadId);
		return ret;
	}

	SizeType nthreads_;
	ObserverHelperType& helper_;
	CorrelationsSkeletonType& skeleton_;
	bool verbose_;
};  //class MultiPointCorrelations
} // namespace Dmrg

/*@}*/
#endif // MULTI_POINT_CORRELATIONS_H
