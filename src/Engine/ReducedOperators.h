/*
Copyright (c) 2009-2015, 2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file ReducedOperators.h
 *
 *  FIXME
 *
 */
#ifndef REDUCEDOP_IMPL_H
#define REDUCEDOP_IMPL_H

#include "Su2SymmetryGlobals.h"
#include "Operator.h"
#include "ChangeOfBasis.h"
#include "BlockOffDiagMatrix.h"
#include "../KronUtil/MatrixDenseOrSparse.h"

namespace Dmrg {
template<typename BasisType>
class ReducedOperators {

	typedef typename BasisType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef Operator<SparseElementType> OperatorType_;
	typedef typename OperatorType_::StorageType OperatorStorageType;
	typedef typename BasisType::RealType RealType;
	typedef typename OperatorType_::PairType PairType;
	typedef ClebschGordanCached<RealType> ClebschGordanType;
	typedef PsimagLite::Matrix<SparseElementType> DenseMatrixType;
	typedef Su2SymmetryGlobals<RealType> Su2SymmetryGlobalsType;
	typedef typename PsimagLite::Vector<OperatorType_>::Type VectorOperatorType;
	typedef typename PsimagLite::Vector<const OperatorType_*>::Type
	VectorPointerOperatorType;
	typedef typename PsimagLite::Vector<SparseElementType>::Type VectorType;
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef ChangeOfBasis<OperatorStorageType, DenseMatrixType> ChangeOfBasisType;

public:

	typedef typename ChangeOfBasisType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef OperatorType_ OperatorType;


	void clear()
	{
		momentumOfOperators_.clear();
		basisrinverse_.clear();
		reducedOperators_.clear();
		reducedHamiltonian_.clear();
		lfactorLeft_.clear();
		lfactorRight_.clear();
		lfactorHamLeft_.clear();
		lfactorHamRight_.clear();
		reducedMapping_.clear();
		fastBasisLeft_.clear();
		fastBasisRight_.clear();
		flavorIndexCached_.clear();
		changeOfBasis_.clear();
		su2Transform_.clear();
		su2TransformT_.clear();
	}

private:

	PsimagLite::Vector<SizeType>::Type momentumOfOperators_;
	//PsimagLite::Vector<SizeType>::Type basisrinverse_;
	//typename PsimagLite::Vector<OperatorType>::Type reducedOperators_;
	//OperatorStorageType reducedHamiltonian_;
	ChangeOfBasisType changeOfBasis_;
}; // ReducedOperators
}// namespace Dmrg

/*@}*/
#endif

