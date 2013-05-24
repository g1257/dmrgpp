// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
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

/*! \file LinearPrediction.h
 *
 *  Extrapolating a "time" series
 *  see extrapolation.tex for more details
 */
  
#ifndef LINEAR_PREDICTION_H
#define LINEAR_PREDICTION_H
#include "Matrix.h"
#include "LAPACK.h"
#include "BLAS.h"

namespace PsimagLite {

	template<typename FieldType>
	class LinearPrediction {
		typedef Matrix<FieldType> MatrixType;
	public:
		LinearPrediction(const typename Vector<FieldType>::Type& y)
		: y_(y)
		{
			SizeType ysize = y.size();
			if (ysize&1) throw RuntimeError(
				"LinearPrediction::ctor(...): data set must contain an even number of points\n");
			SizeType n = ysize/2;
			MatrixType A(n,n);
			typename Vector<FieldType>::Type B(n);
			computeA(A);
			computeB(B);
			computeD(A,B);
		}

		const FieldType& operator()(SizeType i) const
		{
			return y_[i];
		}

		void predict(SizeType p)
		{
			SizeType n = y_.size();
			for (SizeType i=n;i<n+p;i++) {
				FieldType sum = 0;
				for (SizeType j=0;j<d_.size();j++)
					sum += d_[j]*y_[i-j-1];
				y_.push_back(sum);
			}
		}

	private:
		//! Note: A and B cannot be const. here due to the ultimate
		//! call to BLAS::GEMV
		void computeD(MatrixType& A,typename Vector<FieldType>::Type& B)
		{
			SizeType n = B.size();
			typename Vector<int>::Type ipiv(n); // use signed integers here!!
			int info = 0;
			psimag::LAPACK::GETRF(n, n, &(A(0,0)), n, &(ipiv[0]), info);

			typename Vector<FieldType>::Type work(2);
			int lwork = -1; // query mode
			psimag::LAPACK::GETRI(n, &(A(0,0)), n,  &(ipiv[0]),
					&(work[0]), lwork,info );
			lwork = static_cast<int>(work[0]);
			if (lwork<=0) throw
				RuntimeError("LinearPrediction:: internal error\n");
			work.resize(lwork);
			// actual work:
			psimag::LAPACK::GETRI(n, &(A(0,0)), n,  &(ipiv[0]),
								&(work[0]), lwork,info );

			d_.resize(n);
			psimag::BLAS::GEMV('N',n,n,1.0,&(A(0,0)),n,&(B[0]),1,0.0,&(d_[0]),1);
		}

		void computeA(MatrixType& A) const
		{
			SizeType n = A.n_row();
			for (SizeType l=0;l<n;l++) {
				for (SizeType j=0;j<n;j++) {
					A(l,j) = 0;
					for (SizeType i=n;i<2*n;i++)
						A(l,j) += y_[i-l-1] * y_[i-j-1];
				}
			}
		}

		void computeB(typename Vector<FieldType>::Type& B) const
		{
			SizeType n = B.size();
			for (SizeType l=0;l<n;l++) {
				B[l] = 0;
				for (SizeType i=n;i<2*n;i++)
					B[l] += y_[i-l-1] * y_[i];
			}
		}

		typename Vector<FieldType>::Type y_;
		typename Vector<FieldType>::Type d_;
	}; // class LinearPrediction
} // namespace PsimagLite 

/*@}*/	
#endif // LINEAR_PREDICTION_H

