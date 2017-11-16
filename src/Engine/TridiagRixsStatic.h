/*
Copyright (c) 2009,2017 UT-Battelle, LLC
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
/** \file TridiagRixsStatic.h
 * This class performs the tridiagonalization of a modified Hamiltonian, which is
 * H0 + A_i for i the current site and A a local operator
 * This is needed to implement RixsStatic
*/
#ifndef TRIDIAGRIXSSTATIC_H
#define TRIDIAGRIXSSTATIC_H
#include "ApplyOperatorLocal.h"

namespace Dmrg {

template<typename ModelType,typename LanczosSolverType, typename VectorWithOffsetType>
class TridiagRixsStatic {

	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef PsimagLite::Matrix<SparseElementType> MatrixType;
	typedef typename LanczosSolverType::ParametersSolverType ParametersSolverType;
	typedef ApplyOperatorLocal<LeftRightSuperType, VectorWithOffset<ComplexOrRealType> >
	ApplyOperatorLocalType;
	typedef typename VectorWithOffset<ComplexOrRealType>::VectorSizeType VectorSizeType;

	class MyMatrixVector : public LanczosSolverType::LanczosMatrixType {

		typedef typename LanczosSolverType::LanczosMatrixType BasisType;

	public:

		MyMatrixVector(ModelType const *model,
		               ModelHelperType const *modelHelper,
		               const OperatorType& A,
		               ProgramGlobals::DirectionEnum dir,
		               typename ApplyOperatorLocalType::BorderEnum corner,
		               const VectorSizeType& weights)
		    : BasisType(model, modelHelper),
		      applyOperatorLocal_(modelHelper->leftRightSuper()),
		      A_(A),
		      dir_(dir),
		      corner_(corner),
		      fs_(modelHelper->leftRightSuper().left().electronsVector()), // FIXME CHECK
		      x2_(weights, modelHelper->leftRightSuper().super()),
		      y2_(weights, modelHelper->leftRightSuper().super())
		{}

		template<typename SomeVectorType>
		void matrixVectorProduct(SomeVectorType &x,SomeVectorType const &y) const
		{
			BasisType::matrixVectorProduct(x, y);
			MatrixType fullA;
			crsMatrixToFullMatrix(fullA,A_.data);
			if(isZero(fullA))
				return;
			// add here x += Ay
			x2_.setDataInSector(x,0);
			y2_.setDataInSector(y,0);
			applyOperatorLocal_(x2_, y2_, A_, fs_, dir_, corner_);
			x2_.extract(x, 0);
		}

	private:

		ApplyOperatorLocalType applyOperatorLocal_;
		const OperatorType& A_;
		ProgramGlobals::DirectionEnum dir_;
		typename ApplyOperatorLocalType::BorderEnum corner_;
		FermionSign fs_;
		mutable VectorWithOffset<ComplexOrRealType> x2_;
		mutable VectorWithOffset<ComplexOrRealType> y2_;
	}; // class MyMatrixVector

	typedef MyMatrixVector MyMatrixVectorType;
	typedef PsimagLite::LanczosSolver<ParametersSolverType,
	MyMatrixVectorType,
	typename MyMatrixVectorType::VectorType> MyLanczosSolverType;
	typedef typename MyLanczosSolverType::TridiagonalMatrixType MyTridiagonalMatrixType;

public:

	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type TargetVectorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixComplexOrRealType;
	typedef typename PsimagLite::Vector<MatrixComplexOrRealType>::Type VectorMatrixFieldType;

	TridiagRixsStatic(const LeftRightSuperType& lrs,
	                  const ModelType& model,
	                  InputValidatorType& io,
	                  SizeType site,
	                  ProgramGlobals::DirectionEnum direction)
	    : lrs_(lrs),
	      model_(model),
	      io_(io),
	      A_(io, model, false, "RS:"),
	      direction_(direction)
	{
		SizeType numberOfSites = model.geometry().numberOfSites();

		int site2 = ProgramGlobals::findBorderSiteFrom(site, direction, numberOfSites);
		corner_ = (site2 >= 0) ? ApplyOperatorLocalType::BORDER_YES :
		                         ApplyOperatorLocalType::BORDER_NO;
	}

	void operator()(const VectorWithOffsetType& phi,
	                VectorMatrixFieldType& T,
	                VectorMatrixFieldType& V,
	                VectorSizeType& steps)
	{
		for (SizeType ii = 0; ii < phi.sectors(); ++ii) {
			SizeType i = phi.sector(ii);
			steps[ii] = triDiag(phi, T[ii], V[ii], i);
		}
	}

private:

	SizeType triDiag(const VectorWithOffsetType& phi,
	                 MatrixComplexOrRealType& T,
	                 MatrixComplexOrRealType& V,
	                 SizeType i0)
	{
		VectorSizeType weights(lrs_.super().partition(), 0);
		weights[i0] = phi.effectiveSize(i0);
		SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
		SizeType threadNum = 0;
		SizeType currentTime = 0;
		typename ModelType::ModelHelperType modelHelper(p,lrs_,currentTime,threadNum);
		typename MyLanczosSolverType::LanczosMatrixType lanczosHelper(&model_,
		                                                              &modelHelper,
		                                                              A_,
		                                                              direction_,
		                                                              corner_,
		                                                              weights);

		ParametersSolverType params(io_,"Tridiag");
		params.lotaMemory = true;
		params.threadId = threadNum;

		MyLanczosSolverType lanczosSolver(lanczosHelper,params,&V);

		MyTridiagonalMatrixType ab;
		SizeType total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2,i0);
		lanczosSolver.decomposition(phi2,ab);
		lanczosSolver.buildDenseMatrix(T,ab);
		return lanczosSolver.steps();
	}

	const LeftRightSuperType& lrs_;
	const ModelType& model_;
	InputValidatorType& io_;
	OperatorType A_;
	typename ApplyOperatorLocalType::BorderEnum corner_;
	ProgramGlobals::DirectionEnum direction_;
}; // class TridiagRixsStatic
} // namespace Dmrg
#endif // TRIDIAGRIXSSTATIC_H
