/*
Copyright (c) 2012-2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.]
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
/** \ingroup DMRG */
/*@{*/

/*! \file KronConnections.h
 *
 *
 */

#ifndef KRON_CONNECTIONS_H
#define KRON_CONNECTIONS_H

#include "InitKron.h"
#include "Matrix.h"
#include "Concurrency.h"

namespace Dmrg {

template<typename ModelType,typename ModelHelperType_>
class KronConnections {

	typedef InitKron<ModelType,ModelHelperType_> InitKronType;
	typedef typename InitKronType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename InitKronType::GenGroupType GenGroupType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;
	typedef typename InitKronType::RealType RealType;

	KronConnections(const InitKronType& initKron,VectorType& y,const VectorType& x)
	    : initKron_(initKron),y_(y),x_(x),maxVsize_(0)
	{
		init();
	}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      ConcurrencyType::MutexType*)
	{
		bool useSymmetry = initKron_.useSymmetry();
		VectorType yi(maxVsize_,0.0);
		VectorType xi(maxVsize_, 0.0);
		VectorType xj(maxVsize_, 0.0);
		VectorType yij(maxVsize_,0.0);

		for (SizeType p = 0; p < blockSize; ++p) {
			SizeType outPatch = threadNum*blockSize + p;
			if (outPatch >= total) break;

			SizeType i1 = vstart_[outPatch];
			SizeType i2 = i1 + vsize_[outPatch];

			if (i1 == i2) continue;

			assert(i1 < i2);

			std::fill(yi.begin(), yi.begin() + vsize_[outPatch], 0.0);

			/**
			* Symmetry : Access only upper triangular matrix
			**/
			SizeType jpatchStart = 0;
			SizeType jpatchEnd = outPatch;

			if (!useSymmetry) {
				jpatchEnd = total;
			} else {				
				for (SizeType i = i1; i < i2; i++)
					xi[i - i1] = x_[i];
			}

			for (SizeType inPatch = jpatchStart; inPatch < jpatchEnd; ++inPatch) {
				SizeType j1 = vstart_[inPatch];
				SizeType j2 = j1 + vsize_[inPatch];

				if (j1 == j2) continue;

				assert(j1 < j2);

				for (SizeType j = j1; j < j2; j++)
					xj[j - j1] = x_[j];

				std::fill(yij.begin(), yij.begin() + vsize_[outPatch], 0.0);

				SizeType sizeListK = initKron_.connections();

				for (SizeType k = 0; k < sizeListK; ++k) {
					const MatrixDenseOrSparseType& Ak = initKron_.xc(k)(outPatch, inPatch);
					const MatrixDenseOrSparseType& Bk = initKron_.yc(k)(outPatch, inPatch);

					if (Ak.isZero() || Bk.isZero() == 0)
						continue;

					bool diagonal = (outPatch == inPatch);

					// FIXME: Check that kronMult overwrites yi insead of accumulating
					if (useSymmetry && !diagonal) {
						kronMult(yi, xj, 'n', 'n', Ak, Bk);
						for (SizeType i = 0; i < vsize_[outPatch]; ++i)
							yij[i] += yi[i];

						kronMult(yi, xi, 't', 't', Ak, Bk);
						for (SizeType j = j1; j < j2; j++)
							y_[j] += yi[j-j1];
					} else {
						kronMult(yi, xj, 'n','n', Ak, Bk);
						for (SizeType i = 0; i < vsize_[outPatch]; ++i)
							yij[i] += yi[i];
					}
				} // end loop over k

				for (SizeType i = 0; i< vsize_[outPatch]; ++i)
					yi[i] += yij[i];
			} // end inside patch

			for (SizeType i = i1; i < i2; ++i)
				y_[i] = yi[i - i1];
		} // end outside patch
	}

	void sync()
	{}

private:

	void init()
	{
		SizeType npatches = initKron_.connections();
		vsize_.resize(npatches, 0);
		vstart_.resize(npatches, 0);

		const GenGroupType& istartLeft = initKron_.istartLeft();
		const GenGroupType& istartRight = initKron_.istartRight();

		/**
		  *  Calculating the size of the patches
		  *  and row and column indeces for X and Y
		  *  matrices.
		**/
		SizeType ip = 0;
		for (SizeType ipatch = 0; ipatch < npatches; ipatch++){
			//  No of rows for the Lindex
			SizeType nrowL = istartLeft(ipatch+1) - istartLeft(ipatch);
			// No of rows for the Rindex
			SizeType nrowR  = istartRight(ipatch+1) - istartRight(ipatch);
			vsize_[ipatch] = nrowL*nrowR;
			vstart_[ipatch] = ip;
			ip += vsize_[ipatch]; // ip: Points to start of each patch

			if (maxVsize_ < vsize_[ipatch])
				maxVsize_ = vsize_[ipatch];
		}
	}

	const InitKronType& initKron_;
	VectorType& y_;
	const VectorType& x_;
	SizeType maxVsize_;
	VectorSizeType vstart_;
	VectorSizeType vsize_;
}; //class KronConnections

} // namespace PsimagLite

/*@}*/

#endif // KRON_CONNECTIONS_H
