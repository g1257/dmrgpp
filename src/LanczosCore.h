/*
Copyright (c) 2009-2011-2017, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 2.]
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

/*! \file LanczosCore.h
 *
 *  A class to represent a generic Lanczos Solver
 *
 */

#ifndef PSI_LANCZOS_CORE_H
#define PSI_LANCZOS_CORE_H
#include "LanczosVectors.h"
#include "ProgressIndicator.h"
#include "TridiagonalMatrix.h"
#include <cassert>
#include "Vector.h"
#include "Matrix.h"
#include "Random48.h"
#include "ContinuedFraction.h"

namespace PsimagLite {

//! MatrixType must have the following interface:
//! 	RealType type to indicate the matrix type
//! 	rows() member function to indicate the rank of the matrix
//! 	matrixVectorProduct(typename Vector< RealType>::Type& x,const
//!     typename Vector< RealType>::Type& const y)
//!    	   member function that implements the operation x += Hy

template<typename SolverParametersType, typename MatrixType, typename VectorType_>
class LanczosCore {

	typedef LanczosVectors<MatrixType,VectorType_> LanczosVectorsType;

public:

	typedef VectorType_ VectorType;
	typedef typename LanczosVectorsType::DenseMatrixType DenseMatrixType;
	typedef typename LanczosVectorsType::VectorVectorType VectorVectorType;
	typedef typename SolverParametersType::RealType RealType;
	typedef typename LanczosVectorsType::DenseMatrixRealType DenseMatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef SolverParametersType ParametersSolverType;
	typedef MatrixType LanczosMatrixType;
	typedef typename LanczosVectorsType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename VectorType::value_type VectorElementType;
	typedef ContinuedFraction<TridiagonalMatrixType> PostProcType;

	LanczosCore(const MatrixType& mat,
	            const SolverParametersType& params,
	            bool isReorthoEnabled)
	    : progress_("LanczosCore"),
	      mat_(mat),
	      params_(params),
	      steps_(params.steps),
	      lanczosVectors_(mat_,
	                      params.lotaMemory,
	                      params.steps,
	                      isReorthoEnabled)
	{
		OstringStream msg(std::cout.precision());
		msg()<<"Constructing... mat.rank="<<mat_.rows();
		msg()<<" maximum steps="<<steps_<<" maximum eps="<<params_.tolerance<<" requested";
		progress_.printline(msg, std::cout);
	}

	/**
	  *     In each step of the Lanczos algorithm the values of a[]
	  *     and b[] are computed.
	  *     then a tridiagonal matrix T[j] is formed from the matrix
	  *     T[j-1] as
	  *
	  *            | a[0]  b[0]                            |
	  *            | b[0]  a[1]  b[1]                      |
	  *     T(j) = |       b[1]    .     .                 |
	  *            |               .     .  a[j-2]  b[j-2] |
	  *            |               .     .  b[j-2]  a[j-1] |
	  */
	void decomposition(const VectorType& initVector,
	                   TridiagonalMatrixType& ab,
	                   SizeType excitedForStop)
	{
		SizeType& max_nstep = steps_;
		SizeType matsize = mat_.rows();

		if (initVector.size()!=matsize) {
			String msg("decomposition: vector size ");
			msg += ttos(initVector.size()) + " but matrix size ";
			msg += ttos(matsize) + "\n";
			throw RuntimeError(msg);
		}

		VectorType V2(matsize, 0);  // v2
		VectorType V1(matsize, 0);  // v1
		VectorType V0 = initVector;  // v0

		RealType atmp = 0;
		for (SizeType i = 0; i < matsize; ++i)
			atmp += PsimagLite::real(V0[i]*PsimagLite::conj(V0[i]));

		atmp = 1.0/sqrt(atmp);
		for (SizeType i = 0; i < matsize; ++i) V0[i] *= atmp;

		if (max_nstep > matsize) max_nstep = matsize;
		ab.resize(max_nstep,0);

		bool exitFlag = false;
		SizeType j = 0;
		lanczosVectors_.saveInitialVector(V0);
		lanczosVectors_.prepareMemory(matsize, max_nstep);

		// -- 1st step --
		ab.b(0) = 0.0;							// beta_0 = 0 always
		if (lanczosVectors_.lotaMemory()) lanczosVectors_.saveVector(V0, 0);
		for (SizeType i = 0; i < matsize; ++i) V1[i] = 0.0;
		mat_.matrixVectorProduct(V1, V0); // V1 = H|V0>
		atmp = 0.0;
		for (SizeType i = 0; i < matsize; ++i)
			atmp += PsimagLite::real(V1[i]*PsimagLite::conj(V0[i])); // alpha = <V0|V1>
		ab.a(0) = atmp;

		RealType btmp = 0.0;
		for (SizeType i = 0; i < matsize; ++i) {
			V1[i] -= atmp*V0[i];								// V1 = V1 - alpha*V0
			btmp += PsimagLite::real(V1[i]*PsimagLite::conj(V1[i]));
		}

		btmp = sqrt(btmp);
		ab.b(1) = btmp;			// beta = sqrt(V1*V1)

		if (btmp > 0) {
			btmp = 1.0/btmp;
			for (SizeType i = 0; i < matsize; ++i) V1[i] *= btmp;		// normalize V1
		}

		if (lanczosVectors_.lotaMemory() && 1 < max_nstep)
			lanczosVectors_.saveVector(V1, 1);

		VectorRealType tmpEigs(ab.size(), 0);
		VectorRealType eold(ab.size(), 0);
		RealType deltaMax = 0;

		// ---- Starting the loop -------
		for (j = 1; j < max_nstep; ++j) {

			lanczosVectors_.oneStepDecomposition(V0,V1,V2,ab,j);
			ab.diag(tmpEigs, j + 1);
			const SizeType eigsForError = std::min(excitedForStop,
			                                       static_cast<SizeType>(tmpEigs.size()));

			deltaMax = 0;
			for (SizeType ii = 0; ii < eigsForError; ++ii) {
				const RealType delta = fabs(tmpEigs[ii] - eold[ii]);
				eold[ii] = tmpEigs[ii];
				if (delta > deltaMax) deltaMax = delta;
			}

			if (deltaMax < params_.tolerance) exitFlag=true;
			if (j == max_nstep-1) exitFlag=true;
			if (exitFlag && mat_.rows()<=4) break;
			if (exitFlag && j >= params_.minSteps) break;

			if (lanczosVectors_.lotaMemory())
				lanczosVectors_.saveVector(V2, j + 1);

			for (SizeType i = 0; i < matsize; ++i) {
				V0[i] = V1[i];
				V1[i] = V2[i];
			}
		}

		if (j < max_nstep) {
			max_nstep = j + 1;
			ab.resize(max_nstep);
		}

		lanczosVectors_.resize(max_nstep);

		OstringStream msg(std::cout.precision());
		msg()<<"Decomposition done for mat.rank="<<mat_.rows();
		msg()<<" after "<<j<<" steps";
		if (params_.tolerance>0) msg()<<", actual eps="<<deltaMax;

		progress_.printline(msg, std::cout);

		if (j == max_nstep && j != mat_.rows()) {
			OstringStream msg2(std::cout.precision());
			msg2()<<"WARNING: Maximum number of steps used. ";
			msg2()<<"Increasing this maximum is recommended.";
			progress_.printline(msg2, std::cout);
		}
	}

	void excitedVector(VectorType& z, const DenseMatrixType& ritz, SizeType excited) const
	{
		lanczosVectors_.excitedVector(z, ritz, excited);
	}

	SizeType steps() const {return steps_; }

	const DenseMatrixType& lanczosVectors() const
	{
		const DenseMatrixType* ptr = lanczosVectors_.data();
		if (!ptr)
			throw RuntimeError("LanczosCore::lanczosVectors() called but no data stored\n");

		return *(ptr);
	}

	const MatrixType& matrix() const { return mat_; }

private:

	ProgressIndicator progress_;
	const MatrixType& mat_;
	const SolverParametersType& params_;
	SizeType steps_;
	LanczosVectorsType lanczosVectors_;
}; // class LanczosCore

} // namespace PsimagLite

/*@}*/
#endif

