/*
Copyright (c) 2009,-2014 UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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
/** \file ParallelTriDiag.h
*/

#ifndef PARALLEL_TRIDIAG_H
#define PARALLEL_TRIDIAG_H

#include "Mpi.h"
#include "Concurrency.h"

namespace Dmrg {

template<typename ModelType,typename LanczosSolverType, typename VectorWithOffsetType>
class ParallelTriDiag {

	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Concurrency ConcurrencyType;

public:

	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type TargetVectorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixComplexOrRealType;
	typedef typename PsimagLite::Vector<MatrixComplexOrRealType>::Type VectorMatrixFieldType;

	ParallelTriDiag(const VectorWithOffsetType& phi,
	                VectorMatrixFieldType& T,
	                VectorMatrixFieldType& V,
	                typename PsimagLite::Vector<SizeType>::Type& steps,
	                const LeftRightSuperType& lrs,
	                RealType currentTime,
	                const ModelType& model,
	                InputValidatorType& io)
	    : phi_(phi),
	      T_(T),
	      V_(V),
	      steps_(steps),
	      lrs_(lrs),
	      currentTime_(currentTime),
	      model_(model),
	      io_(io)
	{}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      ConcurrencyType::MutexType*)
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = PsimagLite::Concurrency::npthreads;

		ConcurrencyType::mpiDisableIfNeeded(mpiRank,blockSize,"ParallelTriDiag",total);

		for (SizeType p=0;p<blockSize;p++) {
			SizeType ii = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (ii >= total) continue;

			SizeType i = phi_.sector(ii);
			steps_[ii] = triDiag(phi_,T_[ii],V_[ii],i,threadNum);
		}
	}

private:

	SizeType triDiag(const VectorWithOffsetType& phi,
	                 MatrixComplexOrRealType& T,
	                 MatrixComplexOrRealType& V,
	                 SizeType i0,
	                 SizeType threadNum)
	{
		SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
		typename ModelType::ModelHelperType modelHelper(p,lrs_,currentTime_,threadNum);
		typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);

		typename LanczosSolverType::ParametersSolverType params(io_,"Tridiag");
		params.lotaMemory = true;
		params.threadId = threadNum;

		LanczosSolverType lanczosSolver(lanczosHelper,params,&V);

		TridiagonalMatrixType ab;
		SizeType total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2,i0);
		lanczosSolver.decomposition(phi2,ab);
		lanczosSolver.buildDenseMatrix(T,ab);
		return lanczosSolver.steps();
	}

	const VectorWithOffsetType& phi_;
	VectorMatrixFieldType& T_;
	VectorMatrixFieldType& V_;
	typename PsimagLite::Vector<SizeType>::Type& steps_;
	const LeftRightSuperType& lrs_;
	RealType currentTime_;
	const ModelType& model_;
	InputValidatorType& io_;
}; // class ParallelTriDiag
} // namespace Dmrg

/*@}*/
#endif // PARALLEL_TRIDIAG_H

