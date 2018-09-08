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

/*! \file LanczosSolver.h
 *
 *  A class to represent a generic Lanczos Solver
 *
 */

#ifndef LANCZOSSOLVER_HEADER_H
#define LANCZOSSOLVER_HEADER_H
#include "LanczosVectors.h"
#include "ProgressIndicator.h"
#include "TridiagonalMatrix.h"
#include <cassert>
#include "Vector.h"
#include "Matrix.h"
#include "Random48.h"
#include "ContinuedFraction.h"
#include "LanczosOrDavidsonBase.h"

namespace PsimagLite {

//! MatrixType must have the following interface:
//! 	RealType type to indicate the matrix type
//! 	rows() member function to indicate the rank of the matrix
//! 	matrixVectorProduct(typename Vector< RealType>::Type& x,const
//!     typename Vector< RealType>::Type& const y)
//!    	   member function that implements the operation x += Hy

template<typename SolverParametersType,typename MatrixType,typename VectorType>
class LanczosSolver : public LanczosOrDavidsonBase<SolverParametersType,MatrixType,VectorType> {

	typedef LanczosOrDavidsonBase<SolverParametersType,MatrixType,VectorType> BaseType;
	typedef typename SolverParametersType::RealType RealType;
	typedef LanczosVectors<MatrixType,VectorType> LanczosVectorsType;
	typedef typename LanczosVectorsType::DenseMatrixType DenseMatrixType;
	typedef typename LanczosVectorsType::VectorVectorType VectorVectorType;

public:

	typedef typename LanczosVectorsType::DenseMatrixRealType DenseMatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef SolverParametersType ParametersSolverType;
	typedef MatrixType LanczosMatrixType;
	typedef typename LanczosVectorsType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename VectorType::value_type VectorElementType;
	typedef ContinuedFraction<TridiagonalMatrixType> PostProcType;

	enum {WITH_INFO=1,DEBUG=2,ALLOWS_ZERO=4};

	LanczosSolver(MatrixType const &mat,
	              const SolverParametersType& params)
	    : progress_("LanczosSolver",params.threadId),
	      mat_(mat),
	      params_(params),
	      steps_(params.steps),
	      mode_(WITH_INFO),
	      rng_(343311),
	      lanczosVectors_(mat_,
	                      params.lotaMemory,
	                      params.steps,
	                      BaseType::isReorthoEnabled(params))
	{
		setMode(params.options);
		OstringStream msg;
		msg<<"Constructing... mat.rank="<<mat_.rows();
		msg<<" maximum steps="<<steps_<<" maximum eps="<<params_.tolerance<<" requested";
		progress_.printline(msg,std::cout);
	}

	// FIXME : Deprecate this function
	virtual void computeGroundState(RealType& gsEnergy,
	                                VectorType& z)
	{
		SizeType n =mat_.rows();
		RealType atmp=0.0;
		VectorType y(n);

		for (SizeType i=0;i<n;i++) {
			y[i]=rng_()-0.5;
			atmp += PsimagLite::real(y[i]*PsimagLite::conj(y[i]));
		}
		if (mode_ & DEBUG) {
			computeExcitedStateTest(gsEnergy,z,y,0);
			return;
		}
		atmp = 1.0 / sqrt (atmp);
		for (SizeType i = 0; i < n; i++) y[i] *= atmp;
		computeGroundState(gsEnergy,z,y);
	}

	// FIXME : Deprecate this function
	virtual void computeGroundState(RealType &gsEnergy,
	                                VectorType &z,
	                                const VectorType& initialVector)
	{
		if (mode_ & DEBUG) {
			computeExcitedStateTest(gsEnergy,z,initialVector,0);
			return;
		}

		SizeType n=mat_.rows();
		TridiagonalMatrixType ab;

		decomposition(initialVector, ab);
		typename Vector<RealType>::Type c(steps_);
		groundAllocations(steps_ + 2,c.size() > 0);
		try {
			ground(gsEnergy,steps_, ab, c);
		} catch (std::exception &e) {
			std::cerr<<"MatrixRank="<<n<<"\n";
			throw e;
		}

		lanczosVectors_.hookForZ(z,c,ab);

		String str = "LanczosSolver: computeGroundState: ";
		if (norm(z)<1e-6)
			throw RuntimeError(str + " norm is zero\n");

		if (mode_ & WITH_INFO) info(gsEnergy,initialVector,0,std::cout);
	}

	virtual void computeExcitedState(RealType& gsEnergy,VectorType& z,SizeType excited)
	{
		if (excited == 0) return computeGroundState(gsEnergy,z);

		SizeType n =mat_.rows();
		RealType atmp=0.0;
		VectorType y(n);

		for (SizeType i=0;i<n;i++) {
			y[i]=rng_()-0.5;
			atmp += PsimagLite::real(y[i]*PsimagLite::conj(y[i]));
		}

		if (mode_ & DEBUG) {
			computeExcitedStateTest(gsEnergy,z,y,0);
			return;
		}

		atmp = 1.0 / sqrt (atmp);
		for (SizeType i = 0; i < n; i++) y[i] *= atmp;
		computeExcitedState(gsEnergy,z,y,excited);
	}

	/* calculate states in the original basis given
	 * (1) the eigen-vectors (triEigenVec_) of the tridiagonal Matrix
	 * (2) basis states of Lanczos algorithm
	 * \Psi_0(h) = lanczosVectors_(h,m) triEigenVec_(m,0)
	*/
	virtual void computeAnyState(RealType& eval, VectorType& z, SizeType& excited)
	{
		OstringStream msg1;
		msg1 << " WARNING: Must perform decomposition first!!!! \n";
		SizeType nn=mat_.rows();
		SizeType max_iter = triEigenVec_.rows();
		assert(max_iter==triEigenVec_.cols());

		if (max_iter<2) {
			String msg("WARNING: Must perform decomposition first!!!! \n");
			msg += " I need eigen-vectors of the tridiagonal matrix \n ";
			throw RuntimeError(msg);
		}

		z.resize(nn);
		for (SizeType h=0; h<nn; h++) {
			for (SizeType m=0; m<max_iter; m++) {
				VectorElementType tmp = lanczosVectors_.data(h,m);
				z[h] += tmp*triEigenVec_(m,excited);
			}
		}
		eval = eigs_[excited];
	}

	RealType ExpectationH(VectorType& z1) {
		SizeType nn = z1.size();
		VectorType z2(nn,0.0);

		mat_.matrixVectorProduct(z2, z1); // z2 = H |z1>

		VectorElementType out=0.0;
		for (SizeType h = 0; h < nn; h++) {
				out += z2[h]*PsimagLite::conj(z1[h]); // <z1|z2>
		}

		return PsimagLite::real(out);
	}

	virtual void computeExcitedState(RealType &gsEnergy,
	                                 VectorType &z,
	                                 const VectorType& initialVector,
	                                 SizeType excited)
	{
		if (excited == 0) return computeGroundState(gsEnergy,z,initialVector);

		if (mode_ & DEBUG) {
			computeExcitedStateTest(gsEnergy,z,initialVector,excited);
			return;
		}

		TridiagonalMatrixType ab;

		decomposition(initialVector, ab);
		gsEnergy = ab.excited(z,excited);

		if (mode_ & WITH_INFO) info(gsEnergy,initialVector,excited,std::cout);
	}

	void buildDenseMatrix(DenseMatrixType& T,const TridiagonalMatrixType& ab) const
	{
		ab.buildDenseMatrix(T);
	}

	RealType eigs(int i)
	{
		return eigs_[i];
	}

	void push(TridiagonalMatrixType& ab,const RealType& a,const RealType& b) const
	{
		ab.push(a,b);
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
	                   TridiagonalMatrixType& ab)
	{
		SizeType& max_nstep = steps_;
		SizeType matsize=mat_.rows();
		assert(matsize==mat_.cols());

		if (initVector.size()!=matsize) {
			String msg("decomposition: vector size ");
			msg += ttos(initVector.size()) + " but matrix size ";
			msg += ttos(matsize) + "\n";
			throw RuntimeError(msg);
		}

		VectorType V2(matsize);  // v2
		VectorType V1(matsize);  // v1
		VectorType V0 = initVector;  // v0

		RealType atmp = 0;
		for (SizeType i = 0; i < matsize; i++) {
			V2[i] = 0;
			V1[i] = 0;
			atmp += PsimagLite::real(V0[i]*PsimagLite::conj(V0[i]));
		}

		for (SizeType i = 0; i < matsize; i++) V0[i] /= sqrt(atmp);

		if (max_nstep > matsize) max_nstep = matsize;
		ab.resize(max_nstep,0);

		if (mode_ & ALLOWS_ZERO && lanczosVectors_.isHyZero(V0,ab)) return;

		RealType eold = 100.;
		bool exitFlag=false;
		SizeType j = 0;
		RealType enew = 0;
		lanczosVectors_.saveInitialVector(V0);
		typename Vector<RealType>::Type nullVector;
		groundAllocations(max_nstep + 2,false);
		lanczosVectors_.prepareMemory(matsize, max_nstep);

		// -- 1st step --
		ab.b(0) = 0.0;							// beta_0 = 0 always
		if (lanczosVectors_.lotaMemory()) lanczosVectors_.saveVector(V0, 0);
		for (SizeType i = 0; i < matsize; i++) V1[i] = 0.0;
		mat_.matrixVectorProduct(V1, V0); // V1 = H|V0>
		atmp = 0.0;
		for (SizeType i = 0; i < matsize; i++)
			atmp += PsimagLite::real(V1[i]*PsimagLite::conj(V0[i])); // alpha = <V0|V1>
		ab.a(0) = atmp;

		RealType btmp = 0.0;
		for (SizeType i = 0; i < matsize; i++) {
			V1[i] = V1[i] - atmp*V0[i];								// V1 = V1 - alpha*V0
			btmp += PsimagLite::real(V1[i]*PsimagLite::conj(V1[i]));
		}
		btmp = sqrt(btmp);
		ab.b(1) = btmp;			// beta = sqrt(V1*V1)

		for (SizeType i = 0; i < matsize; i++) V1[i] = V1[i]/btmp;		// normalize V1
		if (lanczosVectors_.lotaMemory()) lanczosVectors_.saveVector(V1, 1);

		// ---- Starting the loop -------
		for (j=1; j < max_nstep; ++j) {

			lanczosVectors_.oneStepDecomposition(V0,V1,V2,ab,j);
			enew = TriDiag('N', ab, j+1); // only need eigen-values here
			if (fabs (enew - eold) < params_.tolerance) exitFlag=true;
			if (j == max_nstep-1) exitFlag=true;
			if (exitFlag && mat_.rows()<=4) break;
			if (exitFlag && j >= params_.minSteps) break;

			if (lanczosVectors_.lotaMemory())
				lanczosVectors_.saveVector(V2, j+1);

			for (SizeType i = 0; i < matsize; i++) {
				V0[i] = V1[i];
				V1[i] = V2[i];
			}

			eold = enew;
		}

		// -- calculate the eigen-vectors of Tridiag Matrix
		enew = TriDiag('V', ab, j+1);

		if (j < max_nstep) {
			max_nstep = j + 1;
			ab.resize(max_nstep);
		}

		lanczosVectors_.resize(max_nstep);

		OstringStream msg;
		msg<<"Decomposition done for mat.rank="<<mat_.rows();
		msg<<" after "<<j<<" steps";
		if (params_.tolerance>0) msg<<", actual eps="<<fabs(enew - eold);

		progress_.printline(msg,std::cout);

		if (j == max_nstep && j != mat_.rows()) {
			OstringStream msg2;
			msg2<<"WARNING: Maximum number of steps used. ";
			msg2<<"Increasing this maximum is recommended.";
			progress_.printline(msg2,std::cout);
		}
	}

	double TriDiag(char jobz, TridiagonalMatrixType& ab, int n){

		int lda=n;
		eigs_.resize(n);

		triEigenVec_.resize(n,lda,0.0);
		for (int i=0;i<n;i++){
			for (int j=0;j<n;j++){
				triEigenVec_(i,j) = 0.0;
			}
		}

		for (int i=0;i<n-1;i++){
			triEigenVec_(i,i) = ab.a(i);
			triEigenVec_(i,i+1) = ab.b(i+1);
			triEigenVec_(i+1,i) = PsimagLite::conj(ab.b(i+1));
		}
		triEigenVec_(n-1,n-1) = ab.a(n-1);

		diag(triEigenVec_,eigs_,jobz);

		return eigs_[0];
	}

	void oneStepDec(VectorType& x,
	                VectorType& y,
	                RealType& atmp,
	                RealType& btmp,
	                RealType& b0,
	                SizeType it) const
	{
		lanczosVectors_.oneStepDecomposition(x,y,atmp,btmp,b0,it);
	}

	SizeType steps() const {return steps_; }

	const DenseMatrixType& lanczosVectors() const
	{
		const DenseMatrixType* ptr = lanczosVectors_.data();
		if (!ptr)
			throw RuntimeError("LanczosSolver::lanczosVectors() called but no data stored\n");

		return *(ptr);
	}

private:

	void setMode(const String& options)
	{
		if (options.find("lanczosdebug")!=String::npos) mode_ |=  DEBUG;
		if (options.find("lanczosAllowsZero")!=String::npos) mode_ |= ALLOWS_ZERO;
	}

	void info(RealType energyTmp,
	          const VectorType& x,
	          SizeType excited,
	          std::ostream& os)
	{
		RealType norma=norm(x);
		SizeType& iter = steps_;

		if (norma<1e-5 || norma>100) {
			std::cerr<<"norma="<<norma<<"\n";
		}

		OstringStream msg;
		msg.precision(os.precision());
		String what = "lowest";
		if (excited > 0) what = ttos(excited) + " excited";
		msg<<"Found "<<what<<" eigenvalue= "<<energyTmp<<" after "<<iter;
		msg<<" iterations, "<<" orig. norm="<<norma;
		if (excited > 0)
			msg<<"\nLanczosSolver: EXPERIMENTAL feature excited > 0 is in use";
		progress_.printline(msg,os);
	}

	void groundAllocations(SizeType n, bool extra)
	{
		if (groundD_.size() != n) {
			groundD_.clear();
			groundD_.resize(n);
		}

		if (groundE_.size() != n) {
			groundE_.clear();
			groundE_.resize(n);
		}

		if (!extra) return;

		if (groundV_.size() != n*n) {
			groundV_.clear();
			groundV_.resize(n*n);
		}
	}

	void ground(RealType &s,
	            int n,
	            const TridiagonalMatrixType& ab,
	            typename Vector<RealType>::Type& gs)
	{
		RealType* vki = 0;
		long int maxCounter = params_.stepsForEnergyConvergence;

		if (gs.size() > 0) {
			std::fill(groundV_.begin(),groundV_.end(),0.0);
			vki = &(groundV_[0]);
			for (int k = 0; k < n; k++, vki += (n + 1)) (*vki) = 1.0;
		}

		for (int i = 0; i < n; i++) {
			groundD_[i] = ab.a(i);
			groundE_[i] = ab.b(i);
		}

		long int intCounter=0;
		int m = 0;
		int l = 0;
		for (; l < n; l++) {
			do {
				intCounter++;
				if (intCounter > maxCounter) {
					std::cerr<<"lanczos: ground: premature exit ";
					std::cerr<<"(may indicate an internal error)\n";
					break;
				}

				for (m = l; m < n - 1; m++) {
					RealType dd = fabs(groundD_[m]) + fabs(groundD_[m + 1]);
					if ((fabs(groundE_[m]) + dd) == dd) break;
				}

				if (m != l) {
					RealType g = (groundD_[l + 1] - groundD_[l])/(2.0*groundE_[l]);
					RealType r = sqrt(g*g + 1.0);
					g = groundD_[m] - groundD_[l] + groundE_[l]/
					        (g + (g >= 0 ? fabs(r) : -fabs(r)));
					RealType p = 0.0;
					RealType c = 1.0;
					int i = m -1;
					for (s = 1.0; i >= l; i--) {
						RealType f = s * groundE_[i];
						RealType h = c * groundE_[i];
						groundE_[i + 1] = (r = sqrt(f * f + g * g));
						if (r == 0.0) {
							groundD_[i + 1] -= p;
							groundE_[m] = 0.0;
							break;
						}

						s = f / r;
						c = g / r;
						g = groundD_[i + 1] - p;
						r = (groundD_[i] - g) * s + 2.0 * c * h;
						groundD_[i + 1] = g + (p = s * r);
						g = c * r - h;
						if (gs.size() > 0) {
							vki = &(groundV_[0]) + i;
							for (int k = 0; k < n; k++, vki += n) {
								f = vki[1];
								vki[1] = s * vki[0] + c * f;
								vki[0] = c * vki[0] - s * f;
							}
						}
					}

					if (r == 0.0 && i >= l) continue;
					groundD_[l] -= p;
					groundE_[l] = g;
					groundE_[m] = 0.0;
				}
			} while (m != l);
		}

		s = groundD_[l = 0];
		        for (int i = 1; i < n; i++)
		        if (groundD_[i] < s) s = groundD_[l = i];

		        if (gs.size() > 0) {
			vki = &(groundV_[0]) + l;
			for (int k = 0; k < n; k++, vki += n) gs[k] = (*vki);
		}

		if (intCounter > maxCounter)
		throw RuntimeError("LanczosSolver::ground(): internal error\n");
	}

	void getColumn(MatrixType const &mat,VectorType& x,SizeType col)
	{
		SizeType n =x.size();
		VectorType y(n);
		for (SizeType i=0;i<n;i++) {
			x[i]=0;
			y[i]=0;
			if (i==col) y[i]=1.0;
		}
		mat.matrixVectorProduct (x, y);
	}

	//! only for debugging:
	void computeExcitedStateTest(RealType&,
	                             VectorType&,
	                             const VectorType&,
	                             SizeType excited)
	{
		SizeType n =mat_.rows();
		Matrix<VectorElementType> a(n,n);
		for (SizeType i=0;i<n;i++) {
			VectorType x(n);
			getColumn(mat_,x,i);
			for (SizeType j=0;j<n;j++) a(i,j)=x[j];
		}
		bool ih  = isHermitian(a,true);
		if (!ih) throw RuntimeError("computeGroundState: Matrix not hermitian\n");

		std::cerr<<"Matrix hermitian="<<ih<<"\n";
		std::cerr<<a;
		std::cerr<<"----------------\n";

		typename Vector<RealType>::Type eigs(a.rows());
		diag(a,eigs,'V');
		for (SizeType i=0;i<a.rows();i++) std::cerr<<a(i,0)<<" ";
		std::cerr<<"\n";
		std::cerr<<"--------------------------------\n";
		std::cerr<<"eigs["<<excited<<"]="<<eigs[excited]<<"\n";

		throw RuntimeError("testing lanczos solver\n");
	}

	ProgressIndicator progress_;
	DenseMatrixRealType triEigenVec_;
	MatrixType const& mat_;
	const SolverParametersType& params_;
	SizeType steps_;
	SizeType mode_;
	Random48<RealType> rng_;
	LanczosVectorsType lanczosVectors_;
	typename Vector<RealType>::Type eigs_;
	VectorRealType groundD_;
	VectorRealType groundE_;
	VectorRealType groundV_;
}; // class LanczosSolver

} // namespace PsimagLite

/*@}*/
#endif

