/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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
#include "LanczosOrDavidsonBase.h"

namespace PsimagLite {

/** MatrixType must have the following interface:
	 * RealType type to indicate the matrix type
	 * rows() member function to indicate the rank of the matrix
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
template<typename SolverParametersType,typename MatrixType_,typename VectorType>
class ChebyshevSolver  {

	typedef LanczosOrDavidsonBase<SolverParametersType,MatrixType_,VectorType> NotBaseType;
	typedef typename SolverParametersType::RealType RealType;
	typedef LanczosVectors<MatrixType_,VectorType> LanczosVectorsType;
	typedef typename LanczosVectorsType::DenseMatrixType DenseMatrixType;
	typedef typename LanczosVectorsType::DenseMatrixRealType DenseMatrixRealType;
	typedef typename LanczosVectorsType::VectorVectorType VectorVectorType;

public:

	typedef SolverParametersType ParametersSolverType;
	typedef MatrixType_ MatrixType;
	typedef TridiagonalMatrix<RealType> TridiagonalMatrixType;
	typedef typename VectorType::value_type VectorElementType;
	typedef ChebyshevSerializer<TridiagonalMatrixType> PostProcType;
	typedef PsimagLite::Random48<RealType> RngType;

	enum {WITH_INFO=1,DEBUG=2,ALLOWS_ZERO=4};

	ChebyshevSolver(MatrixType const &mat,
	                SolverParametersType& params)
	    : progress_("ChebyshevSolver"),
	      mat_(mat),
	      params_(params),
	      mode_(WITH_INFO),
	      rng_(343311),
	      lanczosVectors_(mat,
	                      params.lotaMemory,
	                      params.steps,
	                      NotBaseType::isReorthoEnabled(params))
	{
		params.steps=400;
		setMode(params.options);
		computeAandB();
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Constructing... mat.rank="<<mat_.rows()<<" steps="<<params.steps;
		progress_.printline(msgg, std::cout);
	}

	void computeGroundState(RealType&, VectorType&)
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

	void buildDenseMatrix(DenseMatrixType&,
	                      const TridiagonalMatrixType&) const
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
		VectorType x(initVector.size(), 0.0);
		VectorType y = initVector;

		lanczosVectors_.prepareMemory(y.size(), lanczosVectors_.cols());
		ab.resize(2*params_.steps, 0);
		SizeType cols = lanczosVectors_.cols();
		for (SizeType j = 0; j < cols; ++j) {
			if (lanczosVectors_.lotaMemory())
				lanczosVectors_.saveVector(y, j);

			RealType atmp = 0;
			RealType btmp = 0;
			oneStepDec(x, y, atmp, btmp, j);
			ab.a(j) =     2*atmp - ab.a(0);
			ab.b(j) = 2*btmp - ab.b(0);
		}

		// lanczosVectors_.resize(cols); <--- not needed because all steps are performed
		//                                    and there is no early exit here
	}

	//! atmp = < phi_n | phi_n>
	//! btmp = < phi_n | phi_{n+1}>
	void oneStepDec(VectorType& x,
	                VectorType& y,
	                RealType& atmp,
	                RealType& btmp,
	                SizeType jind) const
	{
		bool isFirst = (jind == 0);
		VectorType z(x.size(),0.0);
		mat_.matrixVectorProduct (z, y); // z+= Hy
		// scale matrix:
		z -= params_.b*y;
		z *= params_.oneOverA;

		RealType val = (isFirst) ? 1.0 : 2.0;

		atmp = 0.0;
		for (SizeType i = 0; i < mat_.rows(); i++)
			atmp += PsimagLite::real(y[i]*PsimagLite::conj(y[i]));

		for (SizeType i = 0; i < mat_.rows(); i++) {
			VectorElementType tmp = val*z[i] - x[i];
			x[i] = y[i];
			y[i] = tmp;
		}

		btmp = 0.0;
		for (SizeType i = 0; i < mat_.rows(); i++)
			btmp += PsimagLite::real(y[i]*PsimagLite::conj(x[i]));
	}

	SizeType steps() const {return params_.steps; }

	const DenseMatrixType& lanczosVectors() const
	{
		const DenseMatrixType* ptr = lanczosVectors_.data();
		if (!ptr)
			err("LanczosSolver::lanczosVectors() called but no data stored\n");
		return *(ptr);
	}

private:

	void unimplemented(const String& s) const
	{
		String s2("Hmmm...this ain't looking good...");
		s2 += String(__FILE__) + " " + ttos(__LINE__) + " ";
		s2 += s;
		throw RuntimeError(s);
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
			//throw RuntimeError("Norm\n");
		}

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Found Energy="<<energyTmp<<" after "<<params_.steps;
		msg<<" iterations, "<<" orig. norm="<<norma;
		progress_.printline(msgg, os);
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
		    : matx_(mat),y_(matx_.rows())
		{}

		SizeType rows() const { return matx_.rows(); }

		void matrixVectorProduct (VectorType &x,const VectorType &y) const
		{
			for (SizeType i=0;i<y_.size();i++) y_[i] = -y[i];
			matx_.matrixVectorProduct(x,y_);
		}

		VectorElementType operator()(SizeType i,SizeType j) const
		{
			return matx_(i,j);
		}

	private:
		const MatrixType& matx_;
		mutable VectorType y_;
	}; // class InternalMatrix

	void computeAandB()
	{
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Asking LanczosSolver to compute spectrum bounds...";
		progress_.printline(msgg, std::cout);

		SolverParametersType params;
		InternalMatrix mat2(mat_);
		RealType eMax = 0;
		LanczosSolver<SolverParametersType,InternalMatrix,VectorType>
		        lanczosSolver2(mat2,params);

		VectorType z2(mat_.rows(),0);
		VectorType init(z2.size());
		PsimagLite::fillRandom(init);
		lanczosSolver2.computeOneState(eMax, z2, init, 0);

		VectorType z(mat_.rows(),0);
		LanczosSolver<SolverParametersType,MatrixType,VectorType>
		        lanczosSolver(mat_,params);
		RealType eMin = 0;
		lanczosSolver.computeOneState(eMin, z, init, 0);

		eMax = -eMax;
		eMax *= 3;
		eMin *= 3;
		assert(eMax-eMin>1e-2);

		params_.oneOverA=2.0/(eMax-eMin);
		params_.b=(eMax+eMin)/2;

		PsimagLite::OstringStream msgg2(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
		msg2<<"Spectrum bounds computed, eMax="<<eMax<<" eMin="<<eMin;
		progress_.printline(msgg2, std::cout);
	}

	ProgressIndicator progress_;
	MatrixType const& mat_;
	SolverParametersType& params_;
	SizeType mode_;
	RngType rng_;
	LanczosVectorsType lanczosVectors_;
	//! Scaling factors for the Chebyshev expansion
}; // class ChebyshevSolver
} // namespace PsimagLite
/*@}*/
#endif //CHEBYSHEV_SOLVER_H_

