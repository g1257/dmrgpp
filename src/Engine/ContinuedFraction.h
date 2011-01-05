
/*
Copyright (c) 2009, UT-Battelle, LLC
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

#ifndef CONTINUED_FRACTION_H 
#define CONTINUED_FRACTION_H
#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "LanczosSolver.h"

namespace Dmrg {
	template<
    	typename ModelType_,
	 	typename ConcurrencyType_>
	class ContinuedFraction  {
	public:
		
		typedef ModelType_ ModelType;
		typedef ConcurrencyType_ ConcurrencyType;
		typedef typename ModelType::RealType RealType;
		typedef typename ModelType::VectorType VectorType;
		typedef typename ModelType::SparseMatrixType SparseMatrixType;
		typedef typename VectorType::value_type FieldType;
		typedef typename ModelType::BasisType BasisType;

		typedef LanczosSolver<RealType,SparseMatrixType,VectorType> LanczosSolverType;
		typedef psimag::Matrix<FieldType> MatrixType;
		typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;

		static const size_t parallelRank_ = 0; // ContF needs to support concurrency FIXME 
		
		ContinuedFraction(const ModelType& model)
			: model_(model),progress_("ContinuedFraction",0)
		{
			// printHeader();
			// task 1: Compute Hamiltonian and
			// task 2: Compute ground state |phi>
			computeGroundState();
			
		} 
		

		RealType gsEnergy() const
		{
			return gsEnergy_;
		} 

		void getGreenFunction(TridiagonalMatrixType& ab,RealType& norma,
				size_t i,size_t j) const
		{
			// task 3: compute |initVector> =\sum_x c_x|phi>, where
			// c_x are some operator
			VectorType initVector;
			computeInitVector(initVector,i,j);

			// task 4: tridiag H starting with |initVector>
			triDiagonalize(ab,initVector);

			norma = initVector*initVector;
		}
		 
	
	private:
		

		void computeGroundState()
		{
			model_.setupHamiltonian(hamiltonian_);
			MatrixType fm;
			crsMatrixToFullMatrix(fm,hamiltonian_);
			std::cerr<<fm;
			if (!isHermitian(fm)) throw std::runtime_error("Hamiltonian non Hermitian\n");
			//std::cerr<<hamiltonian_;
			std::cerr<<"Done setting up Hamiltonian\n";

			RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;
			size_t parallelRank = 0;

			LanczosSolverType lanczosSolver(hamiltonian_,iter,eps,parallelRank);
			gsVector_.resize(hamiltonian_.rank());
			lanczosSolver.computeGroundState(gsEnergy_,gsVector_);
		} 

		void computeInitVector(VectorType& initVector,size_t i,size_t j) const
		{
			initVector.resize(model_.size());
			VectorType tmpVector(initVector.size());
			size_t spin = ModelType::SPIN_UP;
			size_t destruction = ModelType::DESTRUCTOR;
			SparseMatrixType ci;
			model_.getOperator(ci,destruction,i,spin);
			SparseMatrixType cj;
			model_.getOperator(cj,destruction,j,spin);
			ci.matrixVectorProduct(tmpVector,gsVector_);
			cj.matrixVectorProduct(initVector,gsVector_);
			initVector += tmpVector;
		} 

		void triDiagonalize(TridiagonalMatrixType& ab,const VectorType& initVector) const
		{
			// tridiagonalize starting with tmpVector = c^\dagger_i|gsVector>
			MatrixType V;

			RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;
			size_t parallelRank = 0;

			LanczosSolverType lanczosSolver(hamiltonian_,iter,eps,parallelRank);

			lanczosSolver.tridiagonalDecomposition(initVector,ab,V);
			//ab.buildDenseMatrix(T);
			//return lanczosSolver.steps();
		} 
		
		
		const ModelType& model_;
		ProgressIndicator progress_;
		SparseMatrixType hamiltonian_;
		RealType gsEnergy_;
		VectorType gsVector_; 
	}; // class ContinuedFraction
} // namespace Dmrg

#endif 
