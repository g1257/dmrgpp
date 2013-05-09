// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009-2011, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
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
/** \ingroup PsimagLite */
/*@{*/

/*! \file ChebyshevSolver.h
 *
 *  A class to represent a generic Lanczos Solver
 *
 */

#ifndef CHEBYSHEV_SOLVER_H_
#define CHEBYSHEV_SOLVER_H_
#include "ProgressIndicator.h"
#include "TridiagonalMatrix.h"
#include "Vector.h"
#include "Matrix.h"
#include "Random48.h"
#include "TypeToString.h"
#include "ChebyshevSerializer.h"
#include "LanczosSolver.h"
#include "String.h"

namespace PsimagLite {

	/** MatrixType must have the following interface:
	 * RealType type to indicate the matrix type
	 * rank() member function to indicate the rank of the matrix
	 * matrixVectorProduct(typename Vector< RealType>::Type& x,const 
	 *        typename Vector< RealType>::Type& const y) 
	 *    member function that implements the operation x += Hy
	 * 
	 * SolverParametersType is just a structure with a few things
	 * like eMax, eMin, the spectrum bounds (needed to scale the Hamiltonian)
	 * steps, the number of moments to compute, you can use 
	 * ParametersForSolver class, if you want.
	 * 
	 */
	template<typename SolverParametersType,typename MatrixType,typename VectorType>
	class ChebyshevSolver {

		typedef typename SolverParametersType::RealType RealType;
		typedef LanczosVectors<RealType,MatrixType,VectorType> LanczosVectorsType;

	public:

		typedef SolverParametersType ParametersSolverType;
		typedef MatrixType LanczosMatrixType;
		typedef typename Vector<RealType>::Type TridiagonalMatrixType;
		typedef typename VectorType::value_type VectorElementType;
// 		typedef Matrix<VectorElementType> DenseMatrixType;
		typedef ChebyshevSerializer<RealType,TridiagonalMatrixType> PostProcType;
		typedef PsimagLite::Random48<RealType> RngType;

		enum {WITH_INFO=1,DEBUG=2,ALLOWS_ZERO=4};

		ChebyshevSolver(MatrixType const &mat,
		                SolverParametersType& params,
		                Matrix<VectorElementType>* storageForLanczosVectors=0)
		: progress_("ChebyshevSolver",0),
		  mat_(mat),
		  params_(params),
		  mode_(WITH_INFO),
		  rng_(343311),
		  lanczosVectors_(mat,params.lotaMemory,storageForLanczosVectors)
		{
			params.steps=400;
			setMode(params.options);
			computeAandB();
			PsimagLite::OstringStream msg;
			msg<<"Constructing... mat.rank="<<mat_.rank()<<" steps="<<params.steps;
			progress_.printline(msg,std::cout);
		}

		void computeGroundState(RealType& gsEnergy,VectorType& z)
		{
			unimplemented("computeGroundState");
		}

		void computeGroundState(RealType &gsEnergy,
		                        VectorType &z,
	                            const VectorType& initialVector)
		{
			if (mode_ & DEBUG) {
				computeGroundStateTest(gsEnergy,z,initialVector);
				return;
			}
			unimplemented("computeGroundState");
		}

		void buildDenseMatrix( Matrix<VectorElementType>& T,const TridiagonalMatrixType& ab) const
		{
			unimplemented("buildDenseMatrix");
		}

		void push(TridiagonalMatrixType& ab,const RealType& a,const RealType& b) const
		{
			ab.push_back(a);
			ab.push_back(b);
		}

		//! ab.a contains the even moments
		//! ab.b contains the odd moments
		void decomposition(const VectorType& initVector,
    	                   TridiagonalMatrixType& ab)
		{
			VectorType x(initVector.size(),0.0);
			VectorType y = initVector;

			lanczosVectors_.reset(y.size(),params_.steps);
			ab.resize(2*params_.steps,0);
			for (size_t j=0; j < lanczosVectors_.n_col(); j++) {
				for (size_t i = 0; i < mat_.rank(); i++)
					lanczosVectors_(i,j) = y[i];
				RealType atmp = 0;
				RealType btmp = 0;
				oneStepDecomposition(x,y,atmp,btmp,j==0);
				ab[2*j] = 2*atmp-ab[0];
				ab[2*j+1] = 2*btmp-ab[1];
			}
		}

		//! atmp = < phi_n | phi_n>
		//! btmp = < phi_n | phi_{n+1}>
		void oneStepDecomposition(VectorType& x,
		                          VectorType& y,
		                          RealType& atmp,
		                          RealType& btmp,
		                          bool isFirst) const
		{
			VectorType z(x.size(),0.0);
			mat_.matrixVectorProduct (z, y); // z+= Hy
			// scale matrix:
			z -= params_.b*y;
			z *= params_.oneOverA;

			RealType val = (isFirst) ? 1.0 : 2.0;

			atmp = 0.0;
			for (size_t i = 0; i < mat_.rank(); i++) 
				atmp += std::real(y[i]*std::conj(y[i]));

			for (size_t i = 0; i < mat_.rank(); i++) {
				VectorElementType tmp = val*z[i] - x[i];
				x[i] = y[i];
				y[i] = tmp;
			}

			btmp = 0.0;
			for (size_t i = 0; i < mat_.rank(); i++)
				btmp += std::real(y[i]*std::conj(x[i]));
		}

		size_t steps() const {return params_.steps; }

	private:

		void unimplemented(const String& s) const
		{
			String s2("Hmmm...this ain't looking good...");
			s2 += String(__FILE__) + " " + ttos(__LINE__) + " ";
			s2 += s;
			throw std::runtime_error(s);
		}

		void setMode(const String& options)
		{
			if (options.find("lanczosdebug")!=String::npos)
				mode_ |=  DEBUG;

			if (options.find("lanczosAllowsZero")!=String::npos)
				mode_ |= ALLOWS_ZERO;
		}

		void info(RealType energyTmp,const VectorType& x,std::ostream& os)
		{
			RealType norma=norm(x);

			if (norma<1e-5 || norma>100) {
				std::cerr<<"norma="<<norma<<"\n";
				//throw std::runtime_error("Norm\n");
			}

			PsimagLite::OstringStream msg;
			msg.precision(8);
			msg<<"Found Energy="<<energyTmp<<" after "<<params_.steps;
			msg<<" iterations, "<<" orig. norm="<<norma;
			progress_.printline(msg,os);
		}

		//! only for debugging:
		void computeGroundStateTest(RealType &gsEnergy,
		                            VectorType& z,
		                            const VectorType& initialVector)
		{
			unimplemented("computeGroundStateTest");
		}

		class InternalMatrix {
		public:
			InternalMatrix(const MatrixType& mat)
			: matx_(mat),y_(matx_.rank())
			{}

			size_t rank() const { return matx_.rank(); }

			void matrixVectorProduct (VectorType &x,const VectorType &y) const
			{
				for (size_t i=0;i<y_.size();i++) y_[i] = -y[i];
				matx_.matrixVectorProduct(x,y_);
			}
			
			VectorElementType operator()(size_t i,size_t j) const
			{
				return matx_(i,j);
			}

		private:
			const MatrixType& matx_;
			mutable VectorType y_;
		}; // class InternalMatrix

		void computeAandB()
		{
			PsimagLite::OstringStream msg;
			msg<<"Asking LanczosSolver to compute spectrum bounds...";
			progress_.printline(msg,std::cout);

			SolverParametersType params;
			params.steps = 200;
			params.tolerance = 1e-10;
			params.stepsForEnergyConvergence=10000;

			InternalMatrix mat2(mat_);
			RealType eMax = 0;
			LanczosSolver<SolverParametersType,InternalMatrix,VectorType> 
			                                         lanczosSolver2(mat2,params);

			VectorType z2(mat_.rank(),0);
			lanczosSolver2.computeGroundState(eMax,z2);
			
			VectorType z(mat_.rank(),0);
			LanczosSolver<SolverParametersType,MatrixType,VectorType> 
			                                         lanczosSolver(mat_,params);
			RealType eMin = 0;
			lanczosSolver.computeGroundState(eMin,z);

			eMax = -eMax;
				
			eMax *= 3;
			eMin *= 3;
// 			eMax = 20.1;
// 			eMin = -20.1;
// 			eMax = -eMin;
			assert(eMax-eMin>1e-2);

			params_.oneOverA=2.0/(eMax-eMin);
			params_.b=(eMax+eMin)/2;
			
			PsimagLite::OstringStream msg2;
			msg2<<"Spectrum bounds computed, eMax="<<eMax<<" eMin="<<eMin;
			progress_.printline(msg2,std::cout);
		}

		ProgressIndicator progress_;
		MatrixType const& mat_;
		SolverParametersType& params_;
		size_t mode_;
		RngType rng_;
		LanczosVectorsType lanczosVectors_;
		//! Scaling factors for the Chebyshev expansion
	}; // class ChebyshevSolver
} // namespace PsimagLite
/*@}*/
#endif //CHEBYSHEV_SOLVER_H_
