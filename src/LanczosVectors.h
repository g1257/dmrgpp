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

	SizeType cols() const { return data_->cols(); }

	SizeType rows() const { return data_->rows(); }

	bool lotaMemory() const { return lotaMemory_; }

	void saveInitialVector(const VectorType& y)
	{
		ysaved_ = y;
	}

	void hookForZ(VectorType& z,
	              const typename Vector<RealType>::Type& c,
	              const TridiagonalMatrixType&)
	{
		if (!lotaMemory_) {
			VectorType x(z.size(),0.0);
			VectorType y = ysaved_;
			for (SizeType i = 0; i < z.size(); i++)
				z[i] = 0.0;
			RealType atmp = 0.0;
			for (SizeType j=0; j < c.size(); j++) {
				RealType ctmp = c[j];
				for (SizeType i = 0; i < y.size(); i++)
					z[i] += ctmp * y[i];
				RealType btmp = 0;
				oneStepDecomposition(x,y,atmp,btmp, j);
			}

			return;
		}

		for (SizeType j = 0; j < data_->cols(); j++) {
			RealType ctmp = c[j];
			for (SizeType i = 0; i < data_->rows(); i++) {
				z[i] += ctmp * data_->operator()(i,j);
			}
		}

	}

	// provides a gracious way to exit if Ay == 0 (we assume that then A=0)
	bool isHyZero(const VectorType& y,
	              TridiagonalMatrixType& ab)
	{
		if (!lotaMemory_) return false;

		OstringStream msg;
		msg<<"Testing whether matrix is zero...";
		progress_.printline(msg,std::cout);

		VectorType x(mat_.rows());

		for (SizeType i = 0; i < x.size(); i++) x[i] = 0.0;

		mat_.matrixVectorProduct (x, y); // x+= Hy

		for (SizeType i = 0; i < x.size(); i++)
			if (PsimagLite::real(x[i]*PsimagLite::conj(x[i]))!=0) return false;

		for (SizeType j=0; j < data_->cols(); j++) {
			for (SizeType i = 0; i < mat_.rows(); i++) {
				data_->operator()(i,j) = (i==j) ? 0.0 : 1.1;
			}
			ab.a(j) = 0.0;
			ab.b(j) = 0.0;
		}
		return true;
	}

	void oneStepDecomposition(VectorType& x,
	                          VectorType& y,
	                          RealType& atmp,
	                          RealType& btmp,
	                          SizeType it) const
	{
		mat_.matrixVectorProduct (x, y); // x+= Hy

		atmp = 0.0;
		for (SizeType i = 0; i < mat_.rows(); i++)
			atmp += PsimagLite::real(y[i]*PsimagLite::conj(x[i]));
		btmp = 0.0;
		for (SizeType i = 0; i < mat_.rows(); i++) {
			x[i] -= atmp * y[i];
			btmp += PsimagLite::real(x[i]*PsimagLite::conj(x[i]));
		}

		btmp = sqrt (btmp);

		if (fabs(btmp)<1e-10) {
			for (SizeType i = 0; i < mat_.rows(); i++) {
				VectorElementType tmp = y[i];
				y[i] = x[i];
				x[i] = -btmp * tmp;
			}
			return;
		}

		assert(btmp >= 1e-12);
		RealType inverseBtmp = 1.0/btmp;
		for (SizeType i = 0; i < mat_.rows(); i++) {
			//lanczosVectors(i,j) = y[i];
			VectorElementType tmp = y[i];
			y[i] = x[i] * inverseBtmp;
			x[i] = -btmp * tmp;
		}

		reorthoIfNecessary(x, it);
	}

private:

	void reorthoIfNecessary(VectorType& x, SizeType it) const
	{
		if (!isReorthoEnabled_) return;

		if (!overlap_ || overlap_->size() <= it)
			throw RuntimeError("reorthoIfNecessary failed\n");

		for (SizeType j = 0; j < it; ++j) {
			overlap_->operator[](j) = 0;
			for (SizeType i = 0; i < x.size(); ++i)
				overlap_->operator[](j) += PsimagLite::conj(data_->operator()(i, j))*x[i];
		}

		RealType oldNorm = 0;
		RealType newNorm = 0;
		for (SizeType i = 0; i < x.size(); ++i) {
			ComplexOrRealType sum = 0.0;
			for (SizeType j = 0; j < it; ++j)
				sum -= overlap_->operator[](j)*data_->operator()(i, j);
			ComplexOrRealType tmp = x[i];
			oldNorm += PsimagLite::real(tmp*PsimagLite::conj(tmp));
			tmp += sum;
			x[i] = tmp;
			newNorm += PsimagLite::real(tmp*PsimagLite::conj(tmp));
		}

		RealType factor = oldNorm/newNorm;
		for (SizeType i = 0; i < x.size(); ++i) x[i] *= factor;
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

