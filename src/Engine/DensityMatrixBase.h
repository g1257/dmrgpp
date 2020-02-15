/*
Copyright (c) 2009-2017, UT-Battelle, LLC
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

/*! \file DensityMatrixBase.h
 *
 *
 *
 */
#ifndef DENSITY_MATRIX_BASE_H
#define DENSITY_MATRIX_BASE_H

#include "BlockDiagonalMatrix.h"
#include "PsimagLite.h"

namespace Dmrg {
template<typename TargetingType>
class DensityMatrixBase {

public:

	typedef typename TargetingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename TargetingType::TargetVectorType::value_type DensityMatrixElementType;
	typedef typename PsimagLite::Real<DensityMatrixElementType>::Type RealType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Matrix<DensityMatrixElementType> MatrixType;
	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;

	struct Params {

		Params(bool u,
		       ProgramGlobals::DirectionEnum d,
		       bool de,
		       bool enablePersistentSvd_,
		       bool serialSvd_)
		    : useSvd(u),
		      direction(d),
		      debug(de),
		      enablePersistentSvd(enablePersistentSvd_),
		      serialSvd(serialSvd_)
		{}

		bool useSvd;
		ProgramGlobals::DirectionEnum direction;
		bool debug;
		bool enablePersistentSvd;
		bool serialSvd;
	};

	typedef typename BlockDiagonalMatrixType::BuildingBlockType BuildingBlockType;

	virtual ~DensityMatrixBase()
	{}

	virtual const BlockDiagonalMatrixType& operator()()=0;

	virtual void diag(typename PsimagLite::Vector<RealType>::Type&, char) = 0;

	virtual const typename PsimagLite::Vector<MatrixType>::Type& vts() const
	{
		return vtsEmpty_;
	}

	virtual const typename PsimagLite::Vector<VectorRealType>::Type& s() const
	{
		return sEmpty_;
	}

	virtual const typename BasisWithOperatorsType::VectorQnType& qns() const
	{
		return qnsEmpty_;
	}

private:

	typename PsimagLite::Vector<MatrixType>::Type vtsEmpty_;
	typename PsimagLite::Vector<VectorRealType>::Type sEmpty_;
	typename BasisWithOperatorsType::VectorQnType qnsEmpty_;

}; // class DensityMatrixBase
} // namespace Dmrg

/*@}*/
#endif

