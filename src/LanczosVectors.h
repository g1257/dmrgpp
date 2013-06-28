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
#include "String.h"

namespace PsimagLite {

template<typename MatrixType,typename VectorType>
class LanczosVectors {

	typedef typename VectorType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef LanczosVectors<MatrixType,VectorType> ThisType;

public:

	typedef TridiagonalMatrix<RealType> TridiagonalMatrixType;
	typedef typename VectorType::value_type VectorElementType;
	typedef PsimagLite::Matrix<VectorElementType> DenseMatrixType;
	typedef PsimagLite::ContinuedFraction<TridiagonalMatrixType> PostProcType;

	enum {WITH_INFO=1,DEBUG=2,ALLOWS_ZERO=4};

	LanczosVectors(const MatrixType& mat,
	               bool lotaMemory,
	               DenseMatrixType* storage)
	    : progress_("LanczosVectors"),
	      mat_(mat),
	      lotaMemory_(lotaMemory),
	      dummy_(0),
	      needsDelete_(false),
	      ysaved_(0)
	{
		if (storage) {
			data_ = storage;
			return;
		}
		data_ = new DenseMatrixType();
		needsDelete_ = true;
	}

	~LanczosVectors()
	{
		if (needsDelete_) delete data_;
	}

	void resize(SizeType matrixRank,SizeType steps)
	{
		steps_= steps;
		if (!lotaMemory_) return;
		data_->reset(matrixRank,steps);
	}

	void reset(SizeType matrixRank,SizeType steps)
	{
		if (!lotaMemory_) return;
		data_->reset(matrixRank,steps);
	}

	VectorElementType& operator()(SizeType i,SizeType j)
	{
		if (!lotaMemory_) return dummy_;
		return data_->operator()(i,j);
	}

	const VectorElementType& operator()(SizeType i,SizeType j) const
	{
		if (!lotaMemory_) return dummy_;
		return data_->operator()(i,j);
	}

	SizeType n_col() const { return data_->n_col(); }

	SizeType n_row() const { return data_->n_row(); }

	bool lotaMemory() const { return lotaMemory_; }

	void saveInitialVector(const VectorType& y)
	{
		ysaved_ = y;
	}

	void hookForZ(VectorType& z,
	              const typename Vector<RealType>::Type& c,
	              const TridiagonalMatrixType& ab)
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
				oneStepDecomposition(x,y,atmp,btmp);
			}
			return;
		}

		for (SizeType j = 0; j < data_->n_col(); j++) {
			RealType ctmp = c[j];
			for (SizeType i = 0; i < data_->n_row(); i++) {
				z[i] += ctmp * data_->operator()(i,j);
			}
		}

	}

	// provides a gracious way to exit if Ay == 0 (we assume that then A=0)
	bool isHyZero(const VectorType& y,
	              TridiagonalMatrixType& ab)
	{
		if (!lotaMemory_) return false;

		PsimagLite::OstringStream msg;
		msg<<"Testing whether matrix is zero...";
		progress_.printline(msg,std::cout);

		VectorType x(mat_.rank());

		for (SizeType i = 0; i < x.size(); i++) x[i] = 0.0;

		mat_.matrixVectorProduct (x, y); // x+= Hy

		for (SizeType i = 0; i < x.size(); i++)
			if (std::real(x[i]*std::conj(x[i]))!=0) return false;

		for (SizeType j=0; j < data_->n_col(); j++) {
			for (SizeType i = 0; i < mat_.rank(); i++) {
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
	                          RealType& btmp) const
	{
		mat_.matrixVectorProduct (x, y); // x+= Hy

		atmp = 0.0;
		for (SizeType i = 0; i < mat_.rank(); i++)
			atmp += std::real(y[i]*std::conj(x[i]));
		btmp = 0.0;
		for (SizeType i = 0; i < mat_.rank(); i++) {
			x[i] -= atmp * y[i];
			btmp += std::real(x[i]*std::conj(x[i]));
		}

		btmp = sqrt (btmp);

		if (fabs(btmp)<1e-10) {
			for (SizeType i = 0; i < mat_.rank(); i++) {
				VectorElementType tmp = y[i];
				y[i] = x[i];
				x[i] = -btmp * tmp;
			}
			return;
		}

		for (SizeType i = 0; i < mat_.rank(); i++) {
			//lanczosVectors(i,j) = y[i];
			VectorElementType tmp = y[i];
			y[i] = x[i] / btmp;
			x[i] = -btmp * tmp;
		}
	}

private:

	//! copy ctor and assigment operator are invalid
	//! because this class contains a pointer:
	ThisType& operator=(const ThisType& other);
	LanczosVectors(const ThisType& copy);

	ProgressIndicator progress_;
	const MatrixType& mat_;
	bool lotaMemory_;
	VectorElementType dummy_;
	bool needsDelete_;
	SizeType steps_;
	VectorType ysaved_;
	DenseMatrixType* data_;
}; // class LanczosVectors

} // namespace PsimagLite

/*@}*/
#endif // LANCZOS_VECTORS_HEADER_H

