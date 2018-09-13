/*
Copyright (c) 2009-2013-2018, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 5.]
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

/*! \file LanczosVectors.h
 *
 *  to store or not to store lanczos vectors
 *
 */

#ifndef LANCZOS_VECTORS_HEADER_H
#define LANCZOS_VECTORS_HEADER_H
#include "ProgressIndicator.h"
#include <cassert>
#include "Vector.h"
#include "Matrix.h"
#include "Random48.h"
#include "ContinuedFraction.h"

namespace PsimagLite {

template<typename MatrixType_,typename VectorType_>
class LanczosVectors {

	typedef typename VectorType_::value_type ComplexOrRealType;
	typedef typename Real<ComplexOrRealType>::Type RealType;
	typedef LanczosVectors<MatrixType_, VectorType_> ThisType;

public:

	typedef MatrixType_ MatrixType;
	typedef VectorType_ VectorType;
	typedef TridiagonalMatrix<RealType> TridiagonalMatrixType;
	typedef typename VectorType::value_type VectorElementType;
	typedef Matrix<VectorElementType> DenseMatrixType;
	typedef Matrix<RealType> DenseMatrixRealType;
	typedef ContinuedFraction<TridiagonalMatrixType> PostProcType;
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;

	enum {WITH_INFO=1,DEBUG=2,ALLOWS_ZERO=4};

	LanczosVectors(const MatrixType& mat,
	               bool lotaMemory,
	               SizeType steps,
	               bool isReorthoEnabled)
	    : progress_("LanczosVectors"),
	      mat_(mat),
	      lotaMemory_(lotaMemory),
	      isReorthoEnabled_(isReorthoEnabled),
	      dummy_(0),
	      needsDelete_(false),
	      ysaved_(0),
	      data_(0),
	      overlap_(0)
	{
		if (!lotaMemory)
			throw RuntimeError("LanczosVectors: support for lotaMemory=false has been removed\n");

		dealWithOverlapStorage(steps);
	}

	~LanczosVectors()
	{
		delete overlap_;
		overlap_ = 0;

		if (!needsDelete_) return;

		delete data_;
		data_ = 0;
	}

	void saveVector(const VectorType& y, SizeType j)
	{
		if (!lotaMemory_) return;

		SizeType rows = data_->rows();
		for (SizeType i = 0; i < rows; ++i)
			data_->operator()(i, j) = y[i];
	}

	void prepareMemory(SizeType rows, SizeType steps)
	{
		if (!lotaMemory_) return;

		dealWithStorageOfV(rows, steps);

		if (overlap_)
			overlap_->resize(steps, 0);
		else
			overlap_ = new VectorType(steps, 0);
	}

	void resize(SizeType x)
	{
		SizeType rows = data_->rows();
		data_->resize(rows, x);
	}

	const DenseMatrixType* data() const
	{
		return data_;
	}

	DenseMatrixType* data()
	{
		return data_;
	}

	VectorElementType data(SizeType i, SizeType j)
	{
		return data_->operator()(i, j);
	}

	SizeType cols() const { return data_->cols(); }

	SizeType rows() const { return data_->rows(); }

	bool lotaMemory() const { return lotaMemory_; }

	void saveInitialVector(const VectorType& y)
	{
		ysaved_ = y;
	}

	void excitedVector(VectorType& z, const DenseMatrixType& ritz, SizeType excited) const
	{
		SizeType small = data_->cols();
		SizeType big = data_->rows();
		for (SizeType j = 0; j <small; j++) {
			ComplexOrRealType ctmp = ritz(j, excited);
			for (SizeType i = 0; i < big; i++)
				z[i] += ctmp * data_->operator()(i, j);
		}
	}

	void oneStepDecomposition(VectorType& V0,
	                          VectorType& V1,
	                          VectorType& V2,
	                          TridiagonalMatrixType& ab,
	                          SizeType iter) const
	{
		SizeType nn = V1.size();
		for (SizeType h = 0; h < nn; h++) V2[h] = 0.0;
		mat_.matrixVectorProduct(V2, V1); // V2 = H|V1>

		RealType atmp = 0.0;
		for (SizeType h = 0; h < nn; h++)
			atmp += PsimagLite::real(V2[h]*PsimagLite::conj(V1[h])); // <V1|V2>
		ab.a(iter) = atmp;

		RealType btmp = 0.0;
		for (SizeType h = 0; h < nn; h++) {
			V2[h] = V2[h] - ab.a(iter)*V1[h] - ab.b(iter)*V0[h];  // V2 = V2 - alpha*V1 - beta*V0;
			btmp += PsimagLite::real(V2[h]*PsimagLite::conj(V2[h]));
		}

		btmp = sqrt(btmp);
		if (iter + 1 < ab.size())
			ab.b(iter+1) = btmp;	// beta = sqrt(V2*V2)

		for (SizeType i = 0; i < nn; i++) V2[i] = V2[i]/btmp;		// normalize V2

		reorthoIfNecessary(V2, iter);
	}

private:

	void reorthoIfNecessary(VectorType& V2, SizeType iter) const
	{
		if (!isReorthoEnabled_) return;

		if (!overlap_ || overlap_->size() <= iter)
			throw RuntimeError("reorthoIfNecessary failed\n");

		SizeType nn = V2.size();
		//cout << "  Re-orthogonalization " << endl;
		for(SizeType i=0; i<iter+1; i++) {

			VectorElementType rij=0.0;
			for(SizeType h=0; h<nn; h++)
				rij += PsimagLite::conj(data_->operator()(h,i))*V2[h];

			// V2 = V2 - <Vi|V2> Vi -- gram-schmid
			for (SizeType h=0; h<nn; h++)
				V2[h] = V2[h] - data_->operator()(h,i)*rij;
		}

		RealType ntmp = 0.0;
		for (SizeType h=0;h<nn;h++)
			ntmp += PsimagLite::real(V2[h]*PsimagLite::conj(V2[h]));

		ntmp = 1.0/sqrt(ntmp);
		for (SizeType h=0; h<nn; h++) V2[h] = V2[h]*ntmp;
	}

	void dealWithStorageOfV(SizeType rows, SizeType cols)
	{
		if (!lotaMemory_) return;

		if (data_) {
			if (rows != data_->rows() || cols != data_->cols())
				throw RuntimeError("LanczosVectors: data has already been set!\n");
			return;
		}

		data_ = new DenseMatrixType(rows, cols);
		needsDelete_ = true;
		OstringStream msg;
		msg<<"lotaMemory_=true";
		progress_.printline(msg,std::cout);
	}

	void dealWithOverlapStorage(SizeType steps)
	{
		if (!isReorthoEnabled_) return;

		OstringStream msg;
		msg<<"Reortho enabled";
		progress_.printline(msg,std::cout);
		return;

		SizeType maxNstep = std::min(steps , mat_.rows());
		overlap_  = new VectorType(maxNstep, 0);
	}

	//! copy ctor and assigment operator are invalid
	//! because this class contains a pointer:
	ThisType& operator=(const ThisType& other);

	LanczosVectors(const ThisType& copy);

	ProgressIndicator progress_;
	const MatrixType& mat_;
	bool lotaMemory_;
	bool isReorthoEnabled_;
	VectorElementType dummy_;
	bool needsDelete_;
	VectorType ysaved_;
	DenseMatrixType* data_;
	VectorType* overlap_;
}; // class LanczosVectors

} // namespace PsimagLite

/*@}*/
#endif // LANCZOS_VECTORS_HEADER_H

