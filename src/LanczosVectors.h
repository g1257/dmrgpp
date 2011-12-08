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

	template<typename RealType,typename MatrixType,typename VectorType>
	class LanczosVectors {

		typedef LanczosVectors<RealType,MatrixType,VectorType> ThisType;

	public:

		typedef TridiagonalMatrix<RealType> TridiagonalMatrixType;
		typedef typename VectorType::value_type VectorElementType;
		typedef PsimagLite::Matrix<VectorElementType> DenseMatrixType;
		typedef PsimagLite::ContinuedFraction<RealType,TridiagonalMatrixType>
		                    PostProcType;

		enum {WITH_INFO=1,DEBUG=2,ALLOWS_ZERO=4};

		LanczosVectors(const MatrixType& mat,
		               bool lotaMemory,
		               DenseMatrixType* storage)
		: progress_("LanczosVectors",0),
		  mat_(mat),
		  lotaMemory_(lotaMemory),
		  dummy_(0),
		  needsDelete_(false)
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

		size_t nullSize() const
		{
			return (lotaMemory_) ? 0 : steps_;
		}

		void resize(size_t matrixRank,size_t steps)
		{
			z_.clear();
			z_.resize(matrixRank);
			for (size_t i=0;i<z_.size();i++) z_[i] = 0; 
			steps_= steps;
			if (!lotaMemory_) return;
			data_->reset(matrixRank,steps);
		}
		
		void reset(size_t matrixRank,size_t steps)
		{
			if (!lotaMemory_) return;
			data_->reset(matrixRank,steps);
		}

		VectorElementType& operator()(size_t i,size_t j)
		{
			if (!lotaMemory_) return dummy_;
			return data_->operator()(i,j);
		}
		
		const VectorElementType& operator()(size_t i,size_t j) const
		{
			if (!lotaMemory_) return dummy_;
			return data_->operator()(i,j);
		}

		size_t n_col() const { return data_->n_col(); }

		size_t n_row() const { return data_->n_row(); }

		bool lotaMemory() const { return lotaMemory_; }
		
		void hookForZ(const VectorType& y,
		              const RealType& ctmp)
		{
			if (lotaMemory_) return;
			for (size_t i = 0; i < y.size(); i++) z_[i] += ctmp * y[i];
		}
		
		void hookForZ(VectorType& z,const std::vector<RealType>& c)
		{
			if (!lotaMemory_) {
				z=z_;
				return;
			}

			for (size_t j = 0; j < data_->n_col(); j++) {
 				//mat_.matrixVectorProduct (x, y);
 				//atmp = ab.a(j);
				//RealType btmp = ab.b(j);
				RealType ctmp = c[j];
 				for (size_t i = 0; i < data_->n_row(); i++) {
					z_[i] += ctmp * data_->operator()(i,j);
					
					//x[i] -= atmp * y[i];
					//VectorElementType tmp = lanczosVectors(i,j);
					//if (fabs(x[i] - lanczosVectors(i,j))>1e-6) throw std::runtime_error("Different\n");
					//VectorElementType tmp = y[i];
					//y[i] = tmp / btmp;
					//x[i] = -btmp * tmp;
				}
			}
			z=z_;
		}

		// provides a gracious way to exit if Ay == 0 (we assume that then A=0)
		bool isHyZero(const VectorType& y,
		              TridiagonalMatrixType& ab)
		{
			if (!lotaMemory_) return false;

			std::ostringstream msg;
			msg<<"Testing whether matrix is zero...";
			progress_.printline(msg,std::cout);

			VectorType x(mat_.rank());

			for (size_t i = 0; i < x.size(); i++) x[i] = 0.0;

			mat_.matrixVectorProduct (x, y); // x+= Hy

			for (size_t i = 0; i < x.size(); i++)
				if (std::real(x[i]*std::conj(x[i]))!=0) return false;

			for (size_t j=0; j < data_->n_col(); j++) {
				for (size_t i = 0; i < mat_.rank(); i++) {
						data_->operator()(i,j) = (i==j) ? 0.0 : 1.1;
				}
				ab.a(j) = 0.0;
				ab.b(j) = 0.0;
			}
			return true;
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
		VectorType z_;
		size_t steps_;
		DenseMatrixType* data_;
	}; // class LanczosVectors
} // namespace PsimagLite

/*@}*/
#endif // LANCZOS_VECTORS_HEADER_H
