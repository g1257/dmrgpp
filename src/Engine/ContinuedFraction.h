
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
#include "ApplyOperatorLocal.h"

namespace Dmrg {
	template<
		template<typename,typename> class LanczosSolverTemplate,
    	typename ModelType_,
	 	typename ConcurrencyType_,
    	typename IoType_,
       	template<typename> class VectorType>
	class ContinuedFraction  {
	public:
		
		typedef ModelType_ ModelType;
		typedef ConcurrencyType_ ConcurrencyType;
		typedef IoType_ IoType;
		typedef typename ModelType::RealType RealType;
		typedef std::complex<RealType> ComplexType;
		typedef typename ModelType::OperatorsType OperatorsType;
		typedef std::vector<ComplexType> ComplexVectorType;
		typedef LanczosSolverTemplate<InternalProductType,ComplexVectorType> LanczosSolverType;
		typedef std::vector<RealType> VectorType;
		typedef psimag::Matrix<ComplexType> ComplexMatrixType;
		typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
		typedef ComplexVectorType TargetVectorType;

		static const size_t parallelRank_ = 0; // ContF needs to support concurrency FIXME
		
		
		ContinuedFraction(const ModelType& model)
			: model_(model),progress_("ContinuedFraction",0)
		{
			printHeader();
			// task 1: Compute Hamiltonian and c operators
			SparseMatrixType hamiltonian;
			std::vector<OperatorType> creationOps;
			computeHamiltonian(hamiltonian,creationOps);
						
			// task 2: Compute ground state |phi>
			computeGroundState(hamiltonian);
			
			// task 3: compute |initVector> =\sum_x c_x|phi>, where 
			// c_x are some operator
			VectorType initVector;
			computeInitVector(initVector);
			MatrixType T;
			
			// task 4: tridiag H starting with |initVector>
			tridiagonalize(T,initVector);
			
			// task 5: diag. T and store the result
			// to be able to produce GreenFuction(i,j)
			// analytically
			diagonalizeAndStore(T,initVector);
		}
		
		
					RealType greenFunction(size_t i,size_t j,RealType t) const
					{
						throw std::runtime_error("Unimplemented\n");
					}
		
	
	private:
		

		void computeHamiltonian(SparseMatrixType& hamiltonian,
				std::vector<OperatorType>& creationOps)
		{
			BasisDataType q;
			Block block(geometry.numberOfSites());
			for (size_t i=0;i<block.size();i++) block[i] = i;
			model.setNaturalBasis(creationMatrix,hamiltonian,q,block);
		}
		

		void computeGroundState(const SparseMatrixType& hamiltonian)
		{
			RealType eps= 0.01*ProgramGlobals::LanczosTolerance;
			size_t iter= ProgramGlobals::LanczosSteps;
			size_t parallelRank = 0;

			LanczosSolverType lanczosSolver(hamiltonian,iter,eps,parallelRank);

			lanczosSolver.computeGroundState(gsEnergy_,gsVector_);
		}
		

		void computeInitVector(VectorType& initVector)
		{
			VectorType tmpVector;
			creationMatrix[6].multiply(tmpVector,gsVector);
			creationMatrix[8].multiply(initVector,gsVector);
			initVector += tmpVector;
		}
		

		void triDiagonalize(MatrixType& T,const VectorType& initVector)
		{
			// tridiagonalize starting with tmpVector = c^\dagger_i|gsVector>
			TridiagonalMatrixType ab;
			MatrixType V;
			lanczosSolver.tridiagonalDecomposition(initVector,ab,V);
			ab.buildDenseMatrix(T);
			//return lanczosSolver.steps();
		}
		

		void diagonalizeAndStore(MatrixType& T,const VectorType& initVector)
		{
			std::vector<>RealType> eig(T.n_row());
			utils::diag(T,eigs,'V');
			RealType norma = norm2(initVector);
		}
		
		
		
		const ModelType& model_;
		ProgressIndicator progress_;
		RealType gsEnergy_;
		VectorType gsVector_;
		
	};

} // namespace Dmrg

#endif
